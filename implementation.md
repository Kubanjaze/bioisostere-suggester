# Phase 20 — Bioisostere Replacement Suggester

**Version:** 1.1 (as-built)
**Author:** Kerwyn Medrano
**Date:** 2026-03-25
**Track:** Track 1 — Cheminformatics Core
**Tier:** Micro (2–3 hrs)
**API Cost:** $0.00 — pure RDKit + pandas + matplotlib

---

## 1. Project Overview

### Goal

Given a set of compounds, enumerate known bioisostere replacements using a
built-in SMIRKS library, apply each transform to each compound, and score
variants by predicted impact on LogP and TPSA. Rank transforms by mean |ΔLogP|
and flag variants that retain drug-likeness.

```bash
python main.py --input data/compounds.csv
```

Outputs:
- `output/bioisostere_variants.csv` — compound × transform variants with scores
- `output/delta_scatter.png` — scatter: ΔLogP (x) vs ΔTPSA (y) per transform
- `output/mean_impact.png` — bar chart: mean ΔLogP and mean ΔTPSA per transform

### What This Phase Teaches

| Concept | Detail |
|---|---|
| Bioisostere | Structural replacement that maintains similar biological properties |
| SMIRKS | Reaction SMARTS; RDKit `AllChem.ReactionFromSmarts(smirks)` |
| `RunReactants` | `rxn.RunReactants((mol,))` → tuple of product tuples |
| ΔLogP / ΔTPSA | `Descriptors.MolLogP`, `Descriptors.TPSA` on parent vs variant |
| Drug-likeness filter | Post-transform RO5 check (MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10) |

### Domain Context

Bioisosteric replacement is a core med-chem strategy for improving ADMET
properties without abolishing activity. Common examples:
- F → Cl: increases lipophilicity, improves membrane permeability
- Cl → F: decreases lipophilicity, reduces metabolic soft-spot
- OMe → OH: increases TPSA, can improve solubility
- Amide NH → N-Me: blocks H-bond donor, can improve BBB permeation
- Phenol OH → F: removes metabolic soft-spot, reduces TPSA

This phase introduces SMIRKS-driven structural enumeration — a pattern used in
lead optimization tools (e.g., BIOSTER, SwissBioisostere, ChEMBL transform miner).

---

## 2. Architecture

```
bioisostere-suggester/
├── main.py
├── requirements.txt
├── README.md
├── .gitignore
├── data/
│   └── compounds.csv
└── output/
    ├── bioisostere_variants.csv
    ├── delta_scatter.png
    └── mean_impact.png
```

---

## Key Concepts
- **Bioisosteric replacement**: Structural substitution that maintains biological activity while modifying physicochemical properties (e.g., F to Cl, OMe to OH)
- **SMIRKS reaction encoding**: `AllChem.ReactionFromSmarts(smirks)` compiles transformation rules; `rxn.RunReactants((mol,))` applies them to generate product variants
- **Delta-property scoring**: `DeltaLogP = variant_logp - parent_logp` and `DeltaTPSA = variant_tpsa - parent_tpsa` quantify the property impact of each transform
- **RO5 drug-likeness filter**: Post-transform check (MW <= 500, LogP <= 5, HBD <= 5, HBA <= 10) flags variants that retain oral bioavailability
- **Multi-product handling**: `RunReactants` returns multiple products when SMARTS matches multiple atoms (e.g., two aryl F atoms); all unique products are kept as separate rows

---

## 3. Input Format

### `data/compounds.csv`
```csv
compound_name,smiles
analog_A,C=CC(=O)Nc1ccccc1
analog_B,C=CC(=O)Nc1ccc(F)cc1
...
```
- 8–10 acrylamide-aniline analogues with diverse substituents (F, Cl, OMe, OH)
- Each has an amide NH → multiple transforms can fire per molecule

---

## 4. Bioisostere Library

6 transforms, encoded as SMIRKS with metadata:

| name | smirks | description | expected ΔLogP |
|---|---|---|---|
| `aryl_F_to_Cl` | `[c:1][F]>>[c:1]Cl` | Aromatic F → Cl | +0.4 to +0.6 |
| `aryl_Cl_to_F` | `[c:1][Cl]>>[c:1]F` | Aromatic Cl → F | −0.4 to −0.6 |
| `aryl_Cl_to_Br` | `[c:1][Cl]>>[c:1]Br` | Aromatic Cl → Br | +0.3 to +0.5 |
| `aryl_OMe_to_OH` | `[c:1][OX2H0][CH3]>>[c:1][OH]` | OMe → phenol | −1.2 to −1.8 |
| `aryl_OH_to_F` | `[c:1][OX2H1]>>[c:1]F` | Phenol → fluorobenzene | +0.9 to +1.3 |
| `amide_N_methyl` | `[NH:1][C:2](=[O:3])>>[N:1](C)[C:2](=[O:3])` | N-methylation of amide | +0.3 to +0.6 |

---

## 5. Module Specification

### `load_compounds(path)` → pd.DataFrame
- Standard loader: compound_name, smiles required; skip invalid SMILES
- Return df with compound_name, smiles, mol

### `build_bioisostere_library()` → dict[str, dict]
- Returns `{name: {"smirks": str, "description": str}}`
- 6 entries as per Section 4
- Compile each with `AllChem.ReactionFromSmarts(smirks)` and store as "rxn"

### `apply_transforms(df, library)` → pd.DataFrame
- For each compound × each transform:
  - `products = rxn.RunReactants((mol,))` — list of product tuples
  - For each product: `Chem.SanitizeMol(p)`, canonicalize SMILES, deduplicate
  - Skip invalid products (sanitization failure or None)
- Return rows: compound_name, parent_smiles, transform_name, variant_smiles

### `score_variants(variants_df)` → pd.DataFrame
- For each variant: compute LogP + TPSA + MW + HBD + HBA
- Merge parent descriptors; compute ΔLogP = variant_logp − parent_logp; ΔTPSA = variant_tpsa − parent_tpsa
- RO5 pass: MW ≤ 500 AND LogP ≤ 5 AND HBD ≤ 5 AND HBA ≤ 10
- Add `ro5_pass` bool column

### `plot_delta_scatter(scored_df, output_path)`
- One point per (compound, transform) pair
- x = ΔLogP, y = ΔTPSA
- Color by transform_name (categorical, tab10 colormap)
- Legend per transform
- Vertical/horizontal zero lines
- Save 150 dpi

### `plot_mean_impact(scored_df, output_path)`
- Group by transform_name; compute mean ΔLogP and mean ΔTPSA
- 1×2 subplot: left bar chart = mean ΔLogP per transform; right = mean ΔTPSA
- Bar color by sign: green if ΔLogP < 0 (decreasing lipophilicity), red if > 0
- Save 150 dpi

### `main()`
- `--input` (required): compounds CSV
- `--output-dir` (default: output)
- Print summary table: transform_name + n_variants + mean_delta_logp + mean_delta_tpsa + n_ro5_pass

---

## 6. Seed Data Design

### `data/compounds.csv` (~9 compounds)
```
analog_A   C=CC(=O)Nc1ccccc1            — amide NH only (N-methyl fires)
analog_B   C=CC(=O)Nc1ccc(F)cc1         — amide NH + aryl F
analog_C   C=CC(=O)Nc1ccc(Cl)cc1        — amide NH + aryl Cl
analog_D   C=CC(=O)Nc1ccc(OC)cc1        — amide NH + OMe
analog_E   C=CC(=O)Nc1ccc(C)cc1         — amide NH only
analog_F   C=CC(=O)Nc1cccc(F)c1         — amide NH + meta-F
analog_G   C=CC(=O)Nc1ccc(F)c(F)c1      — amide NH + 2× aryl F (2 products from F→Cl)
analog_OH  C=CC(=O)Nc1ccc(O)cc1         — amide NH + phenol
analog_FCl C=CC(=O)Nc1ccc(F)c(Cl)c1    — amide NH + F + Cl (3 transforms fire)
```

Expected total variants: ~25–35 rows (some transforms fire 0× or 2× per compound)

---

## 7. Verification Checklist

```bash
python main.py --input data/compounds.csv

# Expected:
# - bioisostere_variants.csv: compound_name + transform_name + variant_smiles + delta_logp + delta_tpsa + ro5_pass
# - aryl_F_to_Cl: delta_logp ~ +0.5, delta_tpsa ~ 0 for F-containing compounds
# - aryl_OMe_to_OH: delta_logp ~ -1.5, delta_tpsa ~ +20
# - amide_N_methyl: delta_logp ~ +0.5, delta_tpsa ~ -17 (removes NH donor)
# - delta_scatter.png: points clustered by transform (each transform has consistent Δ)
# - mean_impact.png: OMe→OH clearly moves TPSA; Cl→F clearly moves LogP negative
```

---

## 8. Risks / Assumptions / Next Step

**Risks:**
- `RunReactants` can return multiple products if the SMARTS matches multiple atoms
  (e.g., analog_G has 2× F → aryl_F_to_Cl fires twice, yielding 2 mono-Cl products)
  Handle: include all unique products as separate rows
- SMIRKS `[c:1][OX2H0][CH3]` for OMe: requires explicit OMe (no branching). If OMe
  is written as `OC` in SMILES it may not match `[CH3]`. Use `[c:1][OX2][CH3]` (no H count)
  to be safe.
- amide N-methyl: `[NH:1][C:2](=[O:3])` — the acrylamide amide `C=CC(=O)N` has sp3 N;
  use `[NH1:1][C:2](=[O:3])` to avoid matching ring NH and secondary amines incorrectly

**Assumptions:**
- All descriptors computed with RDKit (not experimental); ΔLogP reflects clogP, not measured
- Only aromatic F/Cl targeted (not aliphatic); `[c:1]` SMARTS enforces this
- No stereochemistry handling needed for these transforms
- Parent compound always passes RO5 (by design of acrylamide-aniline dataset)

**Next step:** Phase 21 — QSAR Model with Scaffold Split. Build a regression model
(RandomForest) using ECFP4 fingerprints on a curated IC50 dataset, evaluate with
scaffold-based train/test split, and compare to random split.

---

## 9. Actual Results (v1.1)

### Run command
```bash
PYTHONUTF8=1 python main.py --input data/compounds.csv
```

### Transform summary (9 compounds, 6 transforms)

| transform_name | n_variants | mean_delta_logp | mean_delta_tpsa | n_ro5_pass |
|---|---|---|---|---|
| aryl_Cl_to_F | 2 | −0.51 | 0.0 | 2 |
| aryl_OMe_to_OH | 1 | −0.30 | +11.0 | 1 |
| amide_N_methyl | 9 | +0.02 | −8.8 | 9 |
| aryl_Cl_to_Br | 2 | +0.11 | 0.0 | 2 |
| aryl_OH_to_F | 1 | +0.43 | −20.2 | 1 |
| aryl_F_to_Cl | 5 | +0.51 | 0.0 | 5 |

Total: 20 unique variants; all 20 pass RO5.

### Key insights
- analog_G (two aryl F) correctly yields 2 distinct mono-chloro regioisomers via aryl_F_to_Cl ✓
- aryl_Cl_to_F and aryl_F_to_Cl are symmetric inverses (ΔLogP ≈ ±0.51) ✓
- amide_N_methyl fires on all 9 compounds (all have secondary amide NH1); ΔTPSA = −8.8 Å²
  (NH donor at ~29.1 Å² replaced by N-methyl at ~20.3 Å², difference ≈ 8.8) ✓
- aryl_OH_to_F: ΔTPSA = −20.2 Å² (removes phenol OH H-bond donor entirely) ✓

### Deviations from plan
- **aryl_OMe_to_OH ΔTPSA**: spec estimated +20 Å²; actual = +11.0 Å². Correct: OMe TPSA ≈ 9 Å²,
  phenol OH TPSA ≈ 20 Å², delta ≈ +11. The spec estimate was too high.
- **aryl_OMe_to_OH ΔLogP**: spec estimated −1.2 to −1.8; actual = −0.30. RDKit cLogP for
  this specific acrylamide-aniline scaffold gives a smaller OMe→OH delta. Direction is correct.
- **amide_N_methyl ΔLogP ≈ +0.02** (spec expected +0.3 to +0.6): N-methylation of this
  particular acrylamide amide makes little cLogP change in RDKit — the nitrogen in the
  conjugated C=CC(=O)N system already has limited H-bond donor character.
