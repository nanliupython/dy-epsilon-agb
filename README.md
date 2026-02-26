# Dy epsilon (ε) comparison: AGB models vs meteoritic data

This repo contains a small, self-contained script that:

1. Reads an Excel workbook where **each sheet is one AGB model**.
2. Uses the **INI** row as the reference composition.
3. Uses the **last TDU_\*** row as the processed composition.
4. Mixes (dilutes) processed material with INI using a parameter **f**:
   - `mix = f*TDU_last + (1-f)*INI`
5. Computes **ε(Dy)** with internal normalization **162/164** (exponential law),
   using the **INI ratios from the same sheet** as the reference.
6. Fits **f** for each model by minimizing **weighted χ²** using meteoritic uncertainties.
7. Produces a data-model comparison plot and prints a pandas table of best-fit **f** and ε values.

## Input format

The Excel workbook should have one sheet per model with columns:

- `#Isot` (stage label): includes `INI` and `TDU_1 ... TDU_N`
- `Dy160`, `Dy161`, `Dy162`, `Dy163`, `Dy164`

Values can be in any **internally consistent** abundance unit, because all ε values are
computed relative to the INI composition from the same sheet.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate     # (Windows: .venv\Scripts\activate)
pip install -r requirements.txt
python dy_epsilon_agb.py
```

By default the script looks for `AGB_models.xlsx` in the same folder and runs the
four sheets: `2M_ref`, `2M_new`, `3M_ref`, `3M_new`.

## Files

- `dy_epsilon_agb.py` — main script
- `AGB_models.xlsx` — example input workbook
- `requirements.txt` — minimal dependencies
