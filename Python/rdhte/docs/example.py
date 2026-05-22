"""Empirical illustration for the Python rdhte package.

Uses the bundled Granzier-Pons-Tricaud (2023, AEJ) extract distributed at
`R/rdhte_dataset.csv`.

Granzier, Pons, and Tricaud study coordination and bandwagon effects in
French two-round elections, where candidates must clear a
qualifying-vote threshold in the first round to advance to the runoff.
The rule creates a sharp regression-discontinuity design on every
candidate's first-round margin against that threshold: candidates just
above the cutoff advance, those just below do not. The authors document
substantive heterogeneity across ideology and candidate-strength
dimensions, which this script illustrates via every covariate-form
rdhte exposes.

The bundled extract has 39,534 candidate-race observations with:

    y                outcome: 1 if the candidate advances to the runoff
    x                running variable: first-round margin against the
                     qualifying threshold (cutoff at zero)
    cluster_var      district identifier (cluster-robust inference)
    w_left           1 if the candidate's party is left of center
    w_ideology       unordered party-ideology bucket (4 levels)
    w_strength       continuous proxy for ex-ante candidate strength
    w_strong         1 if above-median strength
    w_strength_qrt   ordered quartile bucket of w_strength

Sections:
  1. Single binary variable
  2. Single categorical variable -- unordered
  3. Single categorical variable -- ordered
  4. Two binary variables -- interaction
  5. Single continuous variable
  6. Interaction effect: binary x continuous
  7. Standalone bandwidth selection (rdbwhte)
  8. Efficiency-improving covariates (covs_eff)
  9. Plotting (plot.rdhte)
 10. Replicating rdrobust
 11. Tables (per-cell + spec comparison + LaTeX export)

Run from the package root:

    python Python/rdhte_illustration.py            # prints results
    python Python/rdhte_illustration.py --save     # also writes example_plot.png

The dataset path is resolved relative to this file so the script works
from any working directory.
"""
from __future__ import annotations

import pathlib
import sys
import os

import pandas as pd

from rdhte import rdhte, rdbwhte, rdhte_lincom
from rdhte.plot import plot

os.environ.setdefault("MPLBACKEND", "Agg")

try:
    from rdrobust import rdrobust
except ImportError:  # pragma: no cover
    rdrobust = None


HERE = pathlib.Path(__file__).resolve().parent
DEFAULT_DATA = HERE.parents[2] / "R" / "rdhte_dataset.csv"


def section(title: str) -> None:
    print(f"\n{'=' * 70}\n{title}\n{'=' * 70}")


def load_data(path: pathlib.Path | None = None) -> pd.DataFrame:
    p = pathlib.Path(path) if path else DEFAULT_DATA
    if not p.exists():
        sys.exit(
            f"Could not find rdhte_dataset.csv at {p}. Pass the path as the "
            f"first positional argument, or copy the file alongside this script."
        )
    return pd.read_csv(p)


def s1_single_binary(df: pd.DataFrame) -> None:
    section("1. Single binary variable")
    # Pass as Categorical (or array) to get one CATE per cell. Passing a
    # string column name routes through patsy's continuous-formula path
    # and produces a slope coefficient instead.
    m = rdhte(y="y", x="x",
              covs_hte=pd.Categorical(df["w_left"]),
              cluster="cluster_var", data=df)
    print(m.summary())
    # Advantage for left-of-center candidates.
    res = rdhte_lincom(m, linfct="`1` - `0` = 0")
    print(res["individual"])
    # bw.joint pitfall.
    print("\n  -- bwjoint=True (cell SEs may distort) --")
    print(rdhte(y="y", x="x",
                covs_hte=pd.Categorical(df["w_left"]),
                cluster="cluster_var", bwjoint=True, data=df).summary())


def s2_categorical_unordered(df: pd.DataFrame) -> None:
    section("2. Single categorical variable -- unordered")
    m = rdhte(y="y", x="x",
              covs_hte=pd.Categorical(df["w_ideology"]),
              cluster="cluster_var", data=df)
    print(m.summary())
    # All non-reference cells indistinguishable from zero?
    res = rdhte_lincom(m, linfct=["`2` = 0", "`3` = 0", "`4` = 0"])
    print("\n  individual tests:")
    print(res["individual"])
    print("\n  joint Wald:")
    print(res["joint"])


def s3_categorical_ordered(df: pd.DataFrame) -> None:
    section("3. Single categorical variable -- ordered")
    m = rdhte(y="y", x="x",
              covs_hte=pd.Categorical(df["w_strength_qrt"], ordered=True),
              cluster="cluster_var", data=df)
    print(m.summary())
    print("  -- bare numeric (no Categorical) is treated as continuous --")
    print(rdhte(y="y", x="x", covs_hte="w_strength_qrt",
                cluster="cluster_var", data=df).summary())


def s4_two_binary(df: pd.DataFrame) -> None:
    section("4. Two binary variables -- interaction")
    m = rdhte(y="y", x="x",
              covs_hte=[pd.Categorical(df["w_left"]),
                        pd.Categorical(df["w_strong"])],
              cluster="cluster_var", data=df)
    print(m.summary())


def s5_single_continuous(df: pd.DataFrame) -> None:
    section("5. Single continuous variable")
    # Passing a 1-D numeric array (or float numpy column) makes rdhte
    # switch to the linear-in-W parametrization tau(w) = beta_T + beta_TW * w.
    m = rdhte(y="y", x="x",
              covs_hte=df["w_strength"].to_numpy(),
              kernel="uni", cluster="cluster_var", data=df)
    print(m.summary())
    res = rdhte_lincom(m, linfct="`T:covs.hte` = 0")
    print(res["individual"])


def s6_binary_by_continuous(df: pd.DataFrame) -> None:
    section("6. Interaction effect: binary x continuous")
    m = rdhte(y="y", x="x", covs_hte="w_left * w_strength",
              cluster="cluster_var", data=df)
    print(m.summary())
    print("\n  joint test (each cell coef = 0):")
    # Formula interaction terms use dot-separated component labels.
    # Coef name for the cross term is 'T:w_left.w_strength'.
    res = rdhte_lincom(m, linfct=[
        "`T` = 0",
        "`T:w_left` = 0",
        "`T:w_strength` = 0",
        "`T:w_left.w_strength` = 0",
    ])
    print(res["joint"])
    # Fixed-bandwidth match against per-cell estimation.
    print("\n  -- fix h=0.1 and reproduce per-cell estimates --")
    m_fix = rdhte(y="y", x="x", covs_hte="w_left * w_strength",
                  h=0.1, cluster="cluster_var", data=df)
    print(m_fix.summary())


def s7_rdbwhte(df: pd.DataFrame) -> None:
    section("7. Standalone bandwidth selection (rdbwhte)")
    bw = rdbwhte(y="y", x="x",
                 covs_hte=pd.Categorical(df["w_ideology"]),
                 cluster="cluster_var", data=df)
    print(f"  cells: {bw.W_lev}")
    print(f"  per-cell h:\n{bw.h}")
    bw_j = rdbwhte(y="y", x="x",
                   covs_hte=pd.Categorical(df["w_ideology"]),
                   cluster="cluster_var", bwjoint=True, data=df)
    print(f"  joint h: {bw_j.h}")


def s8_covs_eff(df: pd.DataFrame) -> None:
    section("8. Efficiency-improving covariates (covs_eff)")
    m_no = rdhte(y="y", x="x",
                 covs_hte=pd.Categorical(df["w_ideology"]),
                 cluster="cluster_var", data=df)
    m_yes = rdhte(y="y", x="x",
                  covs_hte=pd.Categorical(df["w_ideology"]),
                  covs_eff="w_strength",
                  cluster="cluster_var", data=df)
    se_no = m_no.se_rb.round(4)
    se_yes = m_yes.se_rb.round(4)
    pct = (100 * (se_yes - se_no) / se_no).round(1)
    out = pd.DataFrame({
        "cell":       list(m_no.W_lev),
        "se_no_eff":  se_no,
        "se_eff":     se_yes,
        "pct_change": pct,
    })
    print(out.to_string(index=False))


def s9_plot(df: pd.DataFrame, save_to: pathlib.Path | None = None) -> None:
    section("9. Plotting (plot.rdhte)")
    m = rdhte(y="y", x="x",
              covs_hte=pd.Categorical(df["w_ideology"]),
              cluster="cluster_var", data=df)
    p = plot(m, sort=True,
             title="Heterogeneity by ideology bucket",
             ylab="Sharp RD ITT")
    if save_to is not None:
        save_to.parent.mkdir(parents=True, exist_ok=True)
        p.save(filename=str(save_to), width=6, height=4, dpi=120,
               verbose=False)
        print(f"  saved figure to {save_to}")
    else:
        print(f"  built {type(p).__name__}; pass --save to write a PNG")


def s10_replicate_rdrobust(df: pd.DataFrame) -> None:
    section("10. Replicating rdrobust")
    if rdrobust is None:
        print("  (skipping: install `rdrobust` to run this section)")
        return
    print("  -- defaults differ; outputs will not match --")
    print(rdhte(y="y", x="x", data=df).summary())
    print(rdrobust(y=df["y"].to_numpy(), x=df["x"].to_numpy()))
    print("\n  -- match exactly: rho=1, vce=hc3, h fixed --")
    print(rdhte(y="y", x="x", h=0.1, vce="hc3", data=df).summary())
    print(rdrobust(y=df["y"].to_numpy(), x=df["x"].to_numpy(),
                   h=0.1, rho=1, vce="hc3"))


def s11_tables(df: pd.DataFrame, latex_to: pathlib.Path | None = None) -> None:
    """Publication-ready tables from rdhte output.

    Three patterns:
      (A) per-cell table from a single rdhte call -- pulls from .tidy()
      (B) cell x specification comparison (varying vce across columns)
      (C) LaTeX export of (A)
    """
    section("11. Tables")

    # ---- (A) per-cell table via tidy() ----
    m = rdhte(y="y", x="x",
              covs_hte=pd.Categorical(df["w_ideology"]),
              cluster="cluster_var", data=df)
    tab_a = m.tidy()
    pretty_a = (tab_a[["term", "estimate", "std.error",
                       "conf.low", "conf.high",
                       "n.eff.left", "n.eff.right",
                       "h.left", "h.right"]]
                .rename(columns={"term": "Cell",
                                 "estimate": "Estimate",
                                 "std.error": "SE",
                                 "conf.low": "CI_lo",
                                 "conf.high": "CI_hi",
                                 "n.eff.left": "Nh-",
                                 "n.eff.right": "Nh+",
                                 "h.left": "h-",
                                 "h.right": "h+"}))
    print("Table A: per-cell estimates  (covs_hte = w_ideology, vce(cluster))")
    print(pretty_a.to_string(
        index=False,
        formatters={"Estimate": "{:.3f}".format,
                    "SE":       "{:.3f}".format,
                    "CI_lo":    "{:.3f}".format,
                    "CI_hi":    "{:.3f}".format,
                    "h-":       "{:.3f}".format,
                    "h+":       "{:.3f}".format,
                    "Nh-":      "{:.0f}".format,
                    "Nh+":      "{:.0f}".format}))

    # glance() gives a one-row summary (n, p, q, kernel, vce, bwselect, ...).
    print("\nFit summary (glance):")
    print(m.glance().to_string(index=False))

    # ---- (B) cell x specification comparison ----
    specs = {
        "HC0": dict(vce="hc0"),
        "HC2": dict(vce="hc2"),
        "HC3": dict(vce="hc3"),
        "CR1": dict(vce="cr1", cluster="cluster_var"),
    }
    col_pt, col_se = {}, {}
    cell_labels = None
    for label, kw in specs.items():
        mk = rdhte(y="y", x="x",
                   covs_hte=pd.Categorical(df["w_ideology"]),
                   data=df, **kw)
        tdy = mk.tidy()
        if cell_labels is None:
            cell_labels = list(tdy["term"])
        col_pt[label] = list(tdy["estimate"])
        col_se[label] = list(tdy["std.error"])

    mat_pt = pd.DataFrame(col_pt, index=cell_labels)
    mat_se = pd.DataFrame(col_se, index=cell_labels)
    print("\nTable B: per-cell point estimates by variance option")
    print(mat_pt.to_string(float_format="{:.3f}".format))
    print("\n(per-cell SEs)")
    print(mat_se.to_string(float_format="{:.3f}".format))

    # ---- (C) LaTeX export of (A) ----
    # Build the LaTeX by hand (pandas DataFrame.to_latex() in recent
    # pandas pulls in jinja2 via Styler, which is an optional dep we'd
    # rather not require for this illustration).
    lines = [
        r"\begin{table}[h]",
        r"\centering",
        r"\caption{Per-cell RD effects by ideology bucket}",
        r"\label{tab:rdhte-ideology}",
        r"\begin{tabular}{lcccc}",
        r"\toprule",
        r"Cell & Estimate & SE & CI lo & CI hi \\",
        r"\midrule",
    ]
    for _, row in pretty_a.iterrows():
        lines.append(
            f"{row['Cell']} & {row['Estimate']:.3f} & {row['SE']:.3f} "
            f"& {row['CI_lo']:.3f} & {row['CI_hi']:.3f} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    tex = "\n".join(lines) + "\n"

    if latex_to is not None:
        latex_to.write_text(tex, encoding="utf-8")
        print(f"\n  wrote LaTeX table to {latex_to}")
    else:
        print("\nLaTeX source (pass --save to write rdhte_table_A.tex):")
        print(tex)


def main(argv: list[str] | None = None) -> None:
    argv = argv if argv is not None else sys.argv[1:]
    save = "--save" in argv
    argv = [a for a in argv if a != "--save"]
    path = pathlib.Path(argv[0]) if argv else None
    df = load_data(path)

    s1_single_binary(df)
    s2_categorical_unordered(df)
    s3_categorical_ordered(df)
    s4_two_binary(df)
    s5_single_continuous(df)
    s6_binary_by_continuous(df)
    s7_rdbwhte(df)
    s8_covs_eff(df)
    s9_plot(df, HERE / "example_plot.png" if save else None)
    s10_replicate_rdrobust(df)
    s11_tables(df, HERE / "rdhte_table_A.tex" if save else None)


if __name__ == "__main__":
    main()
