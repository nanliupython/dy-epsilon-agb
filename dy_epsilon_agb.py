"""
Dy ε-notation comparison: AGB models vs meteoritic data (Excel input)
--------------------------------------------------------------------

What this does (and only this):
- Reads an Excel workbook where each sheet is one AGB model.
- Uses the row "INI" as the reference composition.
- Uses the last "TDU_*" row as the processed (last TDU) composition.
- Mixes (dilutes) last-TDU with INI by f:
      mix = f*TDU_last + (1-f)*INI
- Computes ε(Dy) with internal normalization 162/164 (exponential law),
  using the INI ratios from the same sheet as the reference.
- Fits the best f for each model by minimizing weighted chi^2 using
  meteoritic uncertainties.
- Makes a data-model comparison plot.
- Outputs a pandas table with best-fit f and ε values (so users can copy/paste).

Input workbook expectation (per sheet):
- A column "#Isot" with stage labels including "INI" and "TDU_1"... "TDU_N"
- Columns: Dy160, Dy161, Dy162, Dy163, Dy164 (values can be in any consistent
  abundance unit; mass vs number is OK because we normalize to INI from the same sheet)

Author: Nan Liu + collaborators
"""

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------------- settings ----------------
MIN_A = 160
ANCHOR_NUM = 162
DENOM = 164  # build ratios to 164, internal-normalize using 162/164

# Meteoritic epsilon data: (epsilon, 1σ). Anchors have None uncertainties.
METEORITE = {
    160: ( 6.07, 0.18),
    161: (-0.57, 0.06),
    162: ( 0.00, None),
    163: (-1.04, 0.03),
    164: ( 0.00, None),
}

# Which isotopes to include in chi^2 (exclude anchors unless you explicitly want them)
FIT_A = [160, 161, 163]

# Bounds for f (dilution / mixing fraction of processed material)
F_BOUNDS = (0.0, 1.0)

# Search resolution
N_GRID = 5001
REFINE_ITERS = 6
REFINE_N = 2001
REFINE_HALF_WINDOW_STEPS = 10

DY_COLS = [f"Dy{a}" for a in range(160, 165)]


def _last_tdu_stage(stages) -> str:
    """Return the last TDU_* stage label from a list-like of stage strings."""
    tdu = []
    for s in stages:
        m = re.match(r"^\s*TDU[_\s]?(\d+)\s*$", str(s))
        if m:
            tdu.append((int(m.group(1)), str(s).strip()))
    if not tdu:
        raise ValueError("No TDU_* rows found.")
    return sorted(tdu)[-1][1]


def load_model_from_sheet(excel_path: str, sheet_name: str):
    """
    Returns:
      ini (pd.Series indexed by A), tdu_last (pd.Series indexed by A), last_stage (str)
    """
    df = pd.read_excel(excel_path, sheet_name=sheet_name)
    if "#Isot" not in df.columns:
        raise ValueError(f"Sheet '{sheet_name}' missing '#Isot' column.")
    for c in DY_COLS:
        if c not in df.columns:
            raise ValueError(f"Sheet '{sheet_name}' missing column '{c}'.")

    df["#Isot"] = df["#Isot"].astype(str).str.strip()
    if not (df["#Isot"] == "INI").any():
        raise ValueError(f"Sheet '{sheet_name}' has no 'INI' row.")

    last_stage = _last_tdu_stage(df["#Isot"].tolist())

    ini_row = df.loc[df["#Isot"] == "INI", DY_COLS].iloc[0].astype(float)
    tdu_row = df.loc[df["#Isot"] == last_stage, DY_COLS].iloc[0].astype(float)

    ini = pd.Series({a: float(ini_row[f"Dy{a}"]) for a in range(160, 165)}).sort_index()
    tdu_last = pd.Series({a: float(tdu_row[f"Dy{a}"]) for a in range(160, 165)}).sort_index()
    return ini, tdu_last, last_stage


def mix_linear(tdu_last: pd.Series, ini: pd.Series, f: float) -> pd.Series:
    """mix = f*TDU_last + (1-f)*INI (same units in both)."""
    f = float(f)
    return f * tdu_last + (1.0 - f) * ini


def compute_eps(mix: pd.Series, ini: pd.Series) -> dict:
    """
    ε computed with internal normalization to ANCHOR_NUM/DENOM using exponential law.
    Works with any consistent abundance unit, because everything is ratioed to INI
    from the same sheet (mass-vs-number factors cancel).
    """
    A_list = sorted([a for a in set(mix.index).intersection(ini.index) if a >= MIN_A])
    if DENOM not in A_list or ANCHOR_NUM not in A_list:
        raise ValueError(f"Missing anchors: need {ANCHOR_NUM} and {DENOM}.")

    Rm = {a: mix[a] / mix[DENOM] for a in A_list}
    Rs = {a: ini[a] / ini[DENOM] for a in A_list}

    # keep finite
    A_list = [
        a for a in A_list
        if np.isfinite(Rs[a]) and Rs[a] != 0.0 and np.isfinite(Rm[a]) and Rm[a] != 0.0
    ]

    beta = np.log(Rs[ANCHOR_NUM] / Rm[ANCHOR_NUM]) / np.log(ANCHOR_NUM / DENOM)

    eps = {}
    for a in A_list:
        Rcorr = Rm[a] * (a / DENOM) ** beta
        eps[a] = (Rcorr / Rs[a] - 1.0) * 10000.0
    return eps


def chi2_for_f(ini: pd.Series, tdu_last: pd.Series, f: float, meteorite_dict=METEORITE, fit_A=FIT_A):
    mix = mix_linear(tdu_last, ini, f)
    eps = compute_eps(mix, ini)
    chi2 = 0.0
    npts = 0
    for a in fit_A:
        obs, sig = meteorite_dict[a]
        if sig is None or sig == 0:
            continue
        pred = eps.get(int(a), np.nan)
        if not np.isfinite(pred):
            continue
        chi2 += ((pred - obs) / sig) ** 2
        npts += 1
    return chi2, npts, eps


def fit_f_by_chi2(
    ini: pd.Series,
    tdu_last: pd.Series,
    meteorite_dict=METEORITE,
    f_bounds=F_BOUNDS,
    fit_A=FIT_A,
    n_grid=N_GRID,
    refine_iters=REFINE_ITERS,
    refine_n=REFINE_N,
    refine_half_window_steps=REFINE_HALF_WINDOW_STEPS,
):
    fmin, fmax = f_bounds
    fs = np.linspace(fmin, fmax, n_grid)

    chi = np.empty_like(fs, dtype=float)
    for i, f in enumerate(fs):
        chi[i], npts, _ = chi2_for_f(ini, tdu_last, f, meteorite_dict, fit_A=fit_A)

    j = int(np.nanargmin(chi))
    f_best = float(fs[j])
    chi_best = float(chi[j])

    # refine around best
    for _ in range(refine_iters):
        step = (fmax - fmin) / (n_grid - 1)
        w = refine_half_window_steps * step
        a = max(fmin, f_best - w)
        b = min(fmax, f_best + w)

        fs2 = np.linspace(a, b, refine_n)
        chi2_vals = np.empty_like(fs2, dtype=float)
        for k, f in enumerate(fs2):
            chi2_vals[k], npts, _ = chi2_for_f(ini, tdu_last, f, meteorite_dict, fit_A=fit_A)

        kbest = int(np.nanargmin(chi2_vals))
        f_best = float(fs2[kbest])
        chi_best = float(chi2_vals[kbest])

        fmin, fmax = a, b

    dof_eff = max(1, npts - 1)
    chi_red = chi_best / dof_eff
    _, _, eps_best = chi2_for_f(ini, tdu_last, f_best, meteorite_dict, fit_A=fit_A)
    return f_best, chi_best, chi_red, eps_best


def run(excel_path: str, models: list[dict]):
    """
    models: list of dicts like:
      {"sheet": "2M_ref", "label": "2 Msun Z=0.01 default", "ls": "--"}
    """
    # meteorite arrays for plotting
    mx, my, myerr = [], [], []
    for a in sorted(METEORITE):
        if a < MIN_A:
            continue
        yy, err = METEORITE[a]
        mx.append(a)
        my.append(yy)
        myerr.append(0.0 if err is None else err)

    results = []
    for spec in models:
        ini, tdu_last, last_stage = load_model_from_sheet(excel_path, spec["sheet"])
        f_best, chi2_min, chi2_red, eps_best = fit_f_by_chi2(ini, tdu_last)

        row = {
            "model": spec.get("label", spec["sheet"]),
            "sheet": spec["sheet"],
            "stage_used": last_stage,
            "f_best": f_best,
            "chi2_min": chi2_min,
            "chi2_red": chi2_red,
        }
        for a in range(160, 165):
            row[f"eps{a}"] = float(eps_best.get(a, np.nan))
        results.append(row)

    out = pd.DataFrame(results).sort_values("chi2_min").reset_index(drop=True)

    # plot
    fig, ax = plt.subplots(figsize=(10.5, 4.6))
    fig.subplots_adjust(right=0.76)

    for spec, r in zip(models, results):
        # recompute curve for plotting using best f
        ini, tdu_last, _ = load_model_from_sheet(excel_path, spec["sheet"])
        mix = mix_linear(tdu_last, ini, r["f_best"])
        eps = compute_eps(mix, ini)
        x = np.array(sorted(eps.keys()), dtype=int)
        y = np.array([eps[a] for a in x], dtype=float)

        ax.plot(
            x, y,
            linestyle=spec.get("ls", "-"),
            label=f"{r['model']} (f={r['f_best']:.3e}, χ²ᵣ={r['chi2_red']:.2f})"
        )

    ax.errorbar(mx, my, yerr=myerr, fmt="s", capsize=2, label="Meteorite (ε)")
    ax.axhline(0.0, linewidth=1)
    ax.set_xlabel("Dy isotope mass number (A)")
    ax.set_ylabel("ε(Dy) (internal norm: 162/164; ref: INI in same sheet)")
    ax.set_xlim(159.5, 164.5)
    ax.set_xticks(sorted(set(mx)))
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0)

    plt.show()
    return out


if __name__ == "__main__":
    # Example configuration (edit labels/linestyles as you like)
    EXCEL = "AGB_models.xlsx"
    MODELS = [
        {"sheet": "2M_ref", "label": "2 Msun default", "ls": "--"},
        {"sheet": "2M_new", "label": "2 Msun new", "ls": "-"},
        {"sheet": "3M_ref", "label": "3 Msun default", "ls": "--"},
        {"sheet": "3M_new", "label": "3 Msun new", "ls": "-"},
    ]
    df = run(EXCEL, MODELS)
    print(df.to_string(index=False))
