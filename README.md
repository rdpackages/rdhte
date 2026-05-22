# Robust Local Polynomial Methods for Heterogenous Treatment Effects in RD Designs

The package `rdhte` implements estimation, inference, and bandwidth selection procedures for heterogeneous treatment effects in Regression Discontinuity (RD) designs using local polynomial methods.

- `rdhte`: point estimation and robust bias-corrected inference for conditional RD treatment effects.
- `rdbwhte`: data-driven bandwidth selection for RD heterogeneous treatment effect estimation.
- `rdhte_lincom`: post-estimation tests for linear combinations of heterogeneous treatment effect parameters.


## Python Implementation

To install/update in Python type:
```
pip install rdhte
```

- Help: [PYPI repository](https://pypi.org/project/rdhte/).

- Examples/data: [rdhte illustration](Python/rdhte_illustration.py), [rdhte data](Python/rdhte_dataset.csv).


## R Implementation

To install/update in R type:
```
install.packages('rdhte')
```

- Help: [R Manual](https://cran.r-project.org/web/packages/rdhte/rdhte.pdf), [CRAN repository](https://cran.r-project.org/package=rdhte).

- Examples/data: [rdhte illustration](R/rdhte_illustration.R), [rdhte data](R/rdhte_dataset.csv).


## Stata Implementation

To install/update in Stata type:
```
net install rdhte, from(https://raw.githubusercontent.com/rdpackages/rdhte/main/stata) replace
```

- Help: [rdhte](stata/rdhte.pdf), [rdbwhte](stata/rdbwhte.pdf), [rdhte_lincom](stata/rdhte_lincom.pdf), [rdhte_plot](stata/rdhte_plot.pdf).

- Replication: [rdhte illustration](stata/rdhte_illustration.do), [rdhte data](stata/rdhte_dataset.dta).


## References

For overviews and introductions, see the [rdpackages website](https://rdpackages.github.io).

### Software and Implementation

- Calonico, Cattaneo, Farrell, Palomba and Titiunik (2026): [rdhte: Conditional Average Treatment Effects in RD Designs](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Palomba-Titiunik_2026_Stata.pdf). Working paper.

### Technical and Methodological

- Calonico, Cattaneo, Farrell, Palomba and Titiunik (2026): [Treatment Effect Heterogeneity in Regression Discontinuity Designs](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Palomba-Titiunik_2026_HTERD.pdf). Working paper. [Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Palomba-Titiunik_2026_HTERD--Supplement.pdf).

### Dataset for Replication

- Granzier, Pons and Tricaud (2023): [Coordination and Bandwagon Effects: How Past Rankings Shape the Behavior of Voters and Candidates](https://doi.org/10.1257/app.20210840), _American Economic Journal: Applied Economics_ 15(4): 177-217.

## Funding

This work was supported in part by the National Science Foundation through grants [SES-2019432](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2019432) and [SES-2241575](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2241575).
