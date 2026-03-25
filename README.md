# bioisostere-suggester

**Phase 20 — Bioisostere Replacement Suggester**
Track 1 — Cheminformatics Core

Enumerates bioisostere variants of input compounds using a built-in SMIRKS library,
scores each variant by ΔLogP and ΔTPSA, and flags variants that pass Lipinski RO5.

## The 6 transforms

| Transform | SMIRKS | Rationale |
|---|---|---|
| `aryl_F_to_Cl` | `[c:1][F]>>[c:1]Cl` | F→Cl increases lipophilicity, improves membrane permeability |
| `aryl_Cl_to_F` | `[c:1][Cl]>>[c:1]F` | Cl→F reduces lipophilicity, removes metabolic soft-spot |
| `aryl_Cl_to_Br` | `[c:1][Cl]>>[c:1]Br` | Cl→Br heavier halogen, slight LogP increase |
| `aryl_OMe_to_OH` | `[c:1][OX2][CH3]>>[c:1][OH]` | OMe→phenol increases TPSA, improves solubility |
| `aryl_OH_to_F` | `[c:1][OX2H1]>>[c:1]F` | Phenol→F removes H-bond donor, reduces TPSA |
| `amide_N_methyl` | `[NH1:1][C:2](=[O:3])>>[N:1](C)[C:2](=[O:3])` | N-methylation blocks H-bond donor, may improve BBB |

## Quickstart

```bash
python -m venv .venv
.venv\Scripts\pip install -r requirements.txt

PYTHONUTF8=1 .venv\Scripts\python main.py --input data/compounds.csv
```

On Windows PowerShell:

```powershell
$env:PYTHONUTF8=1
.venv\Scripts\python main.py --input data/compounds.csv
```

> **RDKit install note:** `pip install rdkit` works on most systems. If it fails,
> use conda: `conda install -c conda-forge rdkit`

## CLI flags

| Flag | Default | Description |
|---|---|---|
| `--input` | required | Compounds CSV (compound_name, smiles) |
| `--output-dir` | output | Output directory |

## Outputs

| File | Description |
|---|---|
| `output/bioisostere_variants.csv` | One row per unique variant: parent + transform + variant SMILES + ΔLogP + ΔTPSA + RO5 flag |
| `output/delta_scatter.png` | Scatter ΔLogP vs ΔTPSA, colored by transform |
| `output/mean_impact.png` | 1×2 bar chart: mean ΔLogP and mean ΔTPSA per transform |

## Notes

- **Computed descriptors**: RDKit cLogP and TPSA are calculated (not experimental).
  ΔLogP/ΔTPSA reflect predicted direction; magnitudes may differ from measured values.
- **Multi-hit transforms**: if a molecule has multiple matching atoms, each unique
  product is included as a separate row. Example: `analog_G` has two aryl fluorines;
  `aryl_F_to_Cl` yields two distinct mono-chloro products.
- **SMIRKS safeguards**:
  - `aryl_OMe_to_OH` uses `[OX2]` (not `[OX2H0]`) to correctly match OMe written
    as `OC` in SMILES input.
  - `amide_N_methyl` uses `[NH1:1]` to target secondary amide NH only, avoiding
    ring NHs and primary amines.
- **RO5 flag**: `ro5_pass` is an informational flag (MW ≤ 500, LogP ≤ 5, HBD ≤ 5,
  HBA ≤ 10). It is a guideline, not a hard filter.
- **Products identical to parent** are skipped automatically.
