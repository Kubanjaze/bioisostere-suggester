import sys
import os

if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors

RDLogger.DisableLog("rdApp.*")


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_compounds(path: str) -> pd.DataFrame:
    """Load CSV with compound_name and smiles; skip invalid/empty SMILES."""
    df = pd.read_csv(path)
    required = {"compound_name", "smiles"}
    if not required.issubset(df.columns):
        raise ValueError(f"CSV must have columns: {required}")

    rows = []
    n_invalid = 0
    for _, row in df.iterrows():
        smi = str(row["smiles"]).strip() if pd.notna(row["smiles"]) else ""
        if not smi:
            n_invalid += 1
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            n_invalid += 1
            print(f"  WARNING: invalid SMILES skipped: {row['compound_name']} ({smi})")
            continue
        rows.append({"compound_name": row["compound_name"], "smiles": smi, "mol": mol})

    total = len(df)
    valid = len(rows)
    print(f"Loaded {total} rows, {valid} valid, {n_invalid} invalid/empty SMILES skipped")
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Bioisostere library
# ---------------------------------------------------------------------------

def build_bioisostere_library() -> "dict[str, dict]":
    """Return dict of 6 bioisostere transforms with compiled RDKit reactions."""
    entries = [
        {
            "name": "aryl_F_to_Cl",
            "smirks": "[c:1][F]>>[c:1]Cl",
            "description": "Aromatic F -> Cl (increases lipophilicity)",
        },
        {
            "name": "aryl_Cl_to_F",
            "smirks": "[c:1][Cl]>>[c:1]F",
            "description": "Aromatic Cl -> F (decreases lipophilicity, reduces metabolic soft-spot)",
        },
        {
            "name": "aryl_Cl_to_Br",
            "smirks": "[c:1][Cl]>>[c:1]Br",
            "description": "Aromatic Cl -> Br (heavier halogen, increases lipophilicity)",
        },
        {
            "name": "aryl_OMe_to_OH",
            # [OX2] (not [OX2H0]) to match "OC" written without explicit H in SMILES
            "smirks": "[c:1][OX2][CH3]>>[c:1][OH]",
            "description": "Aryl OMe -> phenol (increases TPSA, improves solubility)",
        },
        {
            "name": "aryl_OH_to_F",
            "smirks": "[c:1][OX2H1]>>[c:1]F",
            "description": "Phenol -> fluorobenzene (removes H-bond donor, reduces TPSA)",
        },
        {
            "name": "amide_N_methyl",
            # [NH1:1] targets secondary amide NH only, avoids ring NHs and primary amines
            "smirks": "[NH1:1][C:2](=[O:3])>>[N:1](C)[C:2](=[O:3])",
            "description": "N-methylation of amide NH (blocks H-bond donor, may improve BBB)",
        },
    ]

    library = {}
    for entry in entries:
        rxn = AllChem.ReactionFromSmarts(entry["smirks"])
        if rxn is None:
            raise ValueError(f"Failed to compile SMIRKS for {entry['name']}: {entry['smirks']}")
        library[entry["name"]] = {
            "smirks": entry["smirks"],
            "description": entry["description"],
            "rxn": rxn,
        }
    return library


# ---------------------------------------------------------------------------
# Transform application
# ---------------------------------------------------------------------------

def apply_transforms(df: pd.DataFrame, library: "dict[str, dict]") -> pd.DataFrame:
    """Apply all transforms to all parent compounds; return variants DataFrame."""
    rows = []
    for _, parent in df.iterrows():
        mol = parent["mol"]
        parent_smi = parent["smiles"]
        name = parent["compound_name"]
        parent_canonical = Chem.MolToSmiles(mol, canonical=True)

        for transform_name, entry in library.items():
            rxn = entry["rxn"]
            try:
                product_sets = rxn.RunReactants((mol,))
            except Exception:
                continue

            seen_smiles: set = set()
            for product_tuple in product_sets:
                for p in product_tuple:
                    try:
                        Chem.SanitizeMol(p)
                    except Exception:
                        continue
                    smi = Chem.MolToSmiles(p, canonical=True)
                    if smi is None or smi == "" or smi == parent_canonical:
                        continue
                    if smi in seen_smiles:
                        continue
                    seen_smiles.add(smi)
                    rows.append({
                        "compound_name": name,
                        "parent_smiles": parent_smi,
                        "transform_name": transform_name,
                        "transform_description": entry["description"],
                        "variant_smiles": smi,
                    })

    return pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=["compound_name", "parent_smiles", "transform_name",
                 "transform_description", "variant_smiles"]
    )


# ---------------------------------------------------------------------------
# Descriptor computation
# ---------------------------------------------------------------------------

def _descriptors(mol) -> dict:
    return {
        "logp": round(Descriptors.MolLogP(mol), 2),
        "tpsa": round(Descriptors.TPSA(mol), 1),
        "mw": round(Descriptors.ExactMolWt(mol), 1),
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
    }


def compute_parent_descriptors(df: pd.DataFrame) -> pd.DataFrame:
    """Compute LogP + TPSA for each parent compound."""
    rows = []
    for _, row in df.iterrows():
        d = _descriptors(row["mol"])
        rows.append({
            "compound_name": row["compound_name"],
            "parent_smiles": row["smiles"],
            "parent_logp": d["logp"],
            "parent_tpsa": d["tpsa"],
        })
    return pd.DataFrame(rows)


def score_variants(variants_df: pd.DataFrame, parent_desc_df: pd.DataFrame) -> pd.DataFrame:
    """Compute variant descriptors, merge parent descriptors, compute deltas + RO5."""
    if variants_df.empty:
        return variants_df.assign(
            parent_logp=pd.Series(dtype=float),
            parent_tpsa=pd.Series(dtype=float),
            variant_logp=pd.Series(dtype=float),
            variant_tpsa=pd.Series(dtype=float),
            delta_logp=pd.Series(dtype=float),
            delta_tpsa=pd.Series(dtype=float),
            mw=pd.Series(dtype=float),
            hbd=pd.Series(dtype=int),
            hba=pd.Series(dtype=int),
            logp=pd.Series(dtype=float),
            tpsa=pd.Series(dtype=float),
            ro5_pass=pd.Series(dtype=bool),
        )

    variant_rows = []
    for _, row in variants_df.iterrows():
        mol = Chem.MolFromSmiles(row["variant_smiles"])
        if mol is None:
            continue
        d = _descriptors(mol)
        variant_rows.append({
            "compound_name": row["compound_name"],
            "parent_smiles": row["parent_smiles"],
            "transform_name": row["transform_name"],
            "transform_description": row["transform_description"],
            "variant_smiles": row["variant_smiles"],
            "variant_logp": d["logp"],
            "variant_tpsa": d["tpsa"],
            "mw": d["mw"],
            "hbd": d["hbd"],
            "hba": d["hba"],
            "logp": d["logp"],
            "tpsa": d["tpsa"],
        })

    if not variant_rows:
        return pd.DataFrame()

    scored = pd.DataFrame(variant_rows)
    scored = scored.merge(parent_desc_df, on=["compound_name", "parent_smiles"], how="left")
    scored["delta_logp"] = (scored["variant_logp"] - scored["parent_logp"]).round(2)
    scored["delta_tpsa"] = (scored["variant_tpsa"] - scored["parent_tpsa"]).round(1)
    scored["ro5_pass"] = (
        (scored["mw"] <= 500)
        & (scored["logp"] <= 5)
        & (scored["hbd"] <= 5)
        & (scored["hba"] <= 10)
    )

    col_order = [
        "compound_name", "parent_smiles", "transform_name", "transform_description",
        "variant_smiles", "parent_logp", "parent_tpsa", "variant_logp", "variant_tpsa",
        "delta_logp", "delta_tpsa", "mw", "hbd", "hba", "logp", "tpsa", "ro5_pass",
    ]
    return scored[[c for c in col_order if c in scored.columns]]


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_delta_scatter(scored_df: pd.DataFrame, output_path: str) -> None:
    """Scatter ΔLogP vs ΔTPSA colored by transform_name."""
    if scored_df.empty:
        print("  No variants to plot (delta_scatter skipped)")
        return

    transforms = sorted(scored_df["transform_name"].unique())
    cmap = matplotlib.colormaps["tab10"]
    color_map = {t: cmap(i / max(len(transforms) - 1, 1)) for i, t in enumerate(transforms)}

    fig, ax = plt.subplots(figsize=(9, 6))
    for transform in transforms:
        sub = scored_df[scored_df["transform_name"] == transform]
        ax.scatter(
            sub["delta_logp"],
            sub["delta_tpsa"],
            label=transform,
            color=color_map[transform],
            s=70,
            edgecolors="k",
            linewidths=0.4,
            alpha=0.85,
            zorder=3,
        )

    ax.axvline(0, color="k", linewidth=0.7, linestyle="--", alpha=0.5)
    ax.axhline(0, color="k", linewidth=0.7, linestyle="--", alpha=0.5)
    ax.set_xlabel("\u0394LogP (variant \u2212 parent)")
    ax.set_ylabel("\u0394TPSA (variant \u2212 parent)")
    ax.set_title("Bioisostere Impact: \u0394LogP vs \u0394TPSA per Transform")
    ax.legend(fontsize=8, bbox_to_anchor=(1.01, 1), loc="upper left", borderaxespad=0)

    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def plot_mean_impact(scored_df: pd.DataFrame, output_path: str) -> None:
    """1x2 bar chart: mean ΔLogP and mean ΔTPSA per transform."""
    if scored_df.empty:
        print("  No variants to plot (mean_impact skipped)")
        return

    summary = (
        scored_df.groupby("transform_name")
        .agg(mean_delta_logp=("delta_logp", "mean"), mean_delta_tpsa=("delta_tpsa", "mean"))
        .reset_index()
        .sort_values("mean_delta_logp")
    )

    def _bar_colors(values, positive_color="#e53935", negative_color="#43a047"):
        return [negative_color if v < 0 else positive_color for v in values]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), squeeze=False)

    # Left: mean ΔLogP
    ax = axes[0, 0]
    ax.barh(
        summary["transform_name"],
        summary["mean_delta_logp"],
        color=_bar_colors(summary["mean_delta_logp"]),
        edgecolor="k",
        linewidth=0.5,
    )
    ax.axvline(0, color="k", linewidth=0.8)
    ax.set_xlabel("Mean \u0394LogP")
    ax.set_title("Mean \u0394LogP per Transform\n(green = more hydrophilic)")

    # Right: mean ΔTPSA
    ax = axes[0, 1]
    ax.barh(
        summary["transform_name"],
        summary["mean_delta_tpsa"],
        color=_bar_colors(summary["mean_delta_tpsa"]),
        edgecolor="k",
        linewidth=0.5,
    )
    ax.axvline(0, color="k", linewidth=0.8)
    ax.set_xlabel("Mean \u0394TPSA (\u00c5\u00b2)")
    ax.set_title("Mean \u0394TPSA per Transform\n(green = lower polarity)")

    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Bioisostere Replacement Suggester — SMIRKS enumeration + ΔLogP/ΔTPSA scoring",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", required=True, help="Compounds CSV (compound_name, smiles)")
    parser.add_argument("--output-dir", default="output", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("Loading compounds...")
    compounds_df = load_compounds(args.input)

    if compounds_df.empty:
        print("No valid compounds found. Writing empty output and exiting.")
        pd.DataFrame(columns=["compound_name", "parent_smiles", "transform_name",
                               "transform_description", "variant_smiles"]).to_csv(
            os.path.join(args.output_dir, "bioisostere_variants.csv"), index=False
        )
        return

    print("Building bioisostere library...")
    library = build_bioisostere_library()
    print(f"  {len(library)} transforms loaded: {', '.join(library.keys())}")

    print("Applying transforms...")
    variants_df = apply_transforms(compounds_df, library)
    print(f"  {len(variants_df)} raw variant rows generated")

    print("Computing parent descriptors...")
    parent_desc_df = compute_parent_descriptors(compounds_df)

    print("Scoring variants...")
    scored_df = score_variants(variants_df, parent_desc_df)

    csv_path = os.path.join(args.output_dir, "bioisostere_variants.csv")
    scored_df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")

    # Summary table
    if not scored_df.empty:
        summary = (
            scored_df.groupby("transform_name")
            .agg(
                n_variants=("variant_smiles", "count"),
                mean_delta_logp=("delta_logp", "mean"),
                mean_delta_tpsa=("delta_tpsa", "mean"),
                n_ro5_pass=("ro5_pass", "sum"),
            )
            .reset_index()
            .sort_values("mean_delta_logp")
        )
        summary["mean_delta_logp"] = summary["mean_delta_logp"].round(2)
        summary["mean_delta_tpsa"] = summary["mean_delta_tpsa"].round(1)
        print("\nTransform summary:")
        print(summary.to_string(index=False))

    scatter_path = os.path.join(args.output_dir, "delta_scatter.png")
    impact_path = os.path.join(args.output_dir, "mean_impact.png")

    print("\nPlotting delta scatter...")
    plot_delta_scatter(scored_df, scatter_path)

    print("Plotting mean impact...")
    plot_mean_impact(scored_df, impact_path)

    print("\nOutputs:")
    for p in [csv_path, scatter_path, impact_path]:
        print(f"  {p}")
    print("\nDone.")


if __name__ == "__main__":
    main()
