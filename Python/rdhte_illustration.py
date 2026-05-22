###############################################################################
# RDHTE Package
# Empirical Illustration
###############################################################################

###############################################################################
# Data and design
###############################################################################
# This script illustrates `rdhte` on a real-data extract from Granzier,
# Pons, and Tricaud (2023, AEJ: Applied), "Coordination and Bandwagon
# Effects: How Past Rankings Shape the Behavior of Voters and
# Candidates." The authors study French two-round elections, where
# candidates must clear a qualifying-vote threshold in the first round
# to advance to the runoff. The institutional rule creates a sharp RD
# design on every candidate's first-round margin against that threshold:
# candidates just above the cutoff advance, those just below do not.
#
# The bundled extract `rdhte_dataset.csv` has 39,534 candidate-race
# observations with the following variables:
#
#   y                outcome: 1 if the candidate advances to the runoff
#   x                running variable: first-round margin against the
#                    qualifying threshold (cutoff at zero)
#   cluster_var      district identifier (cluster-robust inference)
#   w_left           1 if the candidate's party is left of center
#   w_ideology       unordered party-ideology bucket (4 levels)
#   w_strength       continuous proxy for ex-ante candidate strength
#   w_strong         1 if above-median strength
#   w_strength_qrt   ordered quartile bucket of w_strength

import pandas as pd
from rdrobust import rdrobust
from rdhte import rdbwhte, rdhte, rdhte_lincom
from rdhte.plot import plot

rd_data = pd.read_csv("rdhte_dataset.csv")


###############################################################################
# Single binary variable
###############################################################################

rd_left = rdhte(
    y="y",
    x="x",
    covs_hte=pd.Categorical(rd_data["w_left"]),
    cluster="cluster_var",
    data=rd_data,
)
print(rd_left.summary())

print(rdhte_lincom(rd_left, "`1` - `0` = 0")["individual"])

print(
    rdhte(
        y="y",
        x="x",
        covs_hte=pd.Categorical(rd_data["w_left"]),
        cluster="cluster_var",
        bwjoint=True,
        data=rd_data,
    ).summary()
)


###############################################################################
# Single categorical variable -- unordered
###############################################################################

rd_ideology = rdhte(
    y="y",
    x="x",
    covs_hte=pd.Categorical(rd_data["w_ideology"]),
    cluster="cluster_var",
    data=rd_data,
)
print(rd_ideology.summary())

print(rdhte_lincom(rd_ideology, ["`2` = 0", "`3` = 0", "`4` = 0"])["joint"])


###############################################################################
# Single categorical variable -- ordered
###############################################################################

print(
    rdhte(
        y="y",
        x="x",
        covs_hte=pd.Categorical(rd_data["w_strength_qrt"], ordered=True),
        cluster="cluster_var",
        data=rd_data,
    ).summary()
)

print(
    rdhte(
        y="y",
        x="x",
        covs_hte="w_strength_qrt",
        cluster="cluster_var",
        data=rd_data,
    ).summary()
)


###############################################################################
# Two binary variables -- interaction
###############################################################################

print(
    rdhte(
        y="y",
        x="x",
        covs_hte=[
            pd.Categorical(rd_data["w_left"]),
            pd.Categorical(rd_data["w_strong"]),
        ],
        cluster="cluster_var",
        data=rd_data,
    ).summary()
)


###############################################################################
# Single continuous variable
###############################################################################

rd_strength = rdhte(
    y="y",
    x="x",
    covs_hte="w_strength",
    kernel="uni",
    cluster="cluster_var",
    data=rd_data,
)
print(rd_strength.summary())

bw = rd_strength.h[0, 0]
ols_data = rd_data.loc[rd_data["x"].abs() <= bw].copy()
ols_data["T"] = ols_data["x"] > 0
print(pd.Series({"bandwidth": bw, "observations": len(ols_data)}))


###############################################################################
# Interaction effect: binary x continuous
###############################################################################

rd_interaction = rdhte(
    y="y",
    x="x",
    covs_hte="w_left * w_strength",
    cluster="cluster_var",
    data=rd_data,
)
print(rd_interaction.summary())

print(
    rdhte_lincom(
        rd_interaction,
        [
            "`T` = 0",
            "`T:w_left` = 0",
            "`T:w_strength` = 0",
            "`T:w_left.w_strength` = 0",
        ],
    )["joint"]
)

rd_interaction = rdhte(
    y="y",
    x="x",
    covs_hte="w_left * w_strength",
    h=0.1,
    cluster="cluster_var",
    data=rd_data,
)
print(rd_interaction.summary())

print(
    rdhte(
        y="y",
        x="x",
        covs_hte="w_strength",
        h=0.1,
        cluster="cluster_var",
        subset=rd_data["w_left"] == 1,
        data=rd_data,
    ).summary()
)

print(
    rdhte(
        y="y",
        x="x",
        covs_hte="w_strength",
        h=0.1,
        cluster="cluster_var",
        subset=rd_data["w_left"] == 0,
        data=rd_data,
    ).summary()
)


###############################################################################
# Standalone bandwidth selection (rdbwhte)
###############################################################################

rd_bw = rdbwhte(
    y="y",
    x="x",
    covs_hte=pd.Categorical(rd_data["w_ideology"]),
    cluster="cluster_var",
    data=rd_data,
)
print(rd_bw.h)

print(
    rdbwhte(
        y="y",
        x="x",
        covs_hte=pd.Categorical(rd_data["w_ideology"]),
        cluster="cluster_var",
        bwjoint=True,
        data=rd_data,
    ).h
)


###############################################################################
# Efficiency-improving covariates (covs_eff)
###############################################################################

rd_eff_off = rdhte(
    y="y",
    x="x",
    covs_hte=pd.Categorical(rd_data["w_ideology"]),
    cluster="cluster_var",
    data=rd_data,
)
rd_eff_on = rdhte(
    y="y",
    x="x",
    covs_hte=pd.Categorical(rd_data["w_ideology"]),
    covs_eff="w_strength",
    cluster="cluster_var",
    data=rd_data,
)

print(pd.DataFrame({"no_covs_eff": rd_eff_off.se_rb, "covs_eff": rd_eff_on.se_rb}))


###############################################################################
# Plotting
###############################################################################

print(plot(rd_ideology))
print(plot(rd_ideology, sort=True))
print(plot(rd_ideology, sort=True,
           title="Heterogeneity by ideology bucket",
           ylab="Sharp RD ITT"))


###############################################################################
# Building publication-ready tables
###############################################################################

print(rd_ideology.tidy())
print(rd_ideology.glance())

rd_hc1 = rdhte(y="y", x="x", covs_hte=pd.Categorical(rd_data["w_ideology"]),
               vce="hc1", data=rd_data)
rd_hc2 = rdhte(y="y", x="x", covs_hte=pd.Categorical(rd_data["w_ideology"]),
               vce="hc2", data=rd_data)
rd_hc3 = rdhte(y="y", x="x", covs_hte=pd.Categorical(rd_data["w_ideology"]),
               vce="hc3", data=rd_data)
rd_cr1 = rdhte(y="y", x="x", covs_hte=pd.Categorical(rd_data["w_ideology"]),
               vce="cr1", cluster="cluster_var", data=rd_data)

print(pd.DataFrame({"HC1": rd_hc1.coef, "HC2": rd_hc2.coef,
                    "HC3": rd_hc3.coef, "CR1": rd_cr1.coef}))
print(pd.DataFrame({"HC1": rd_hc1.se_rb, "HC2": rd_hc2.se_rb,
                    "HC3": rd_hc3.se_rb, "CR1": rd_cr1.se_rb}))


###############################################################################
# Replicating rdrobust
###############################################################################

print(rdhte(y="y", x="x", data=rd_data).summary())
print(rdrobust(y=rd_data["y"].to_numpy(), x=rd_data["x"].to_numpy()))

print(rdhte(y="y", x="x", h=0.1, vce="hc3", data=rd_data).summary())
print(rdrobust(y=rd_data["y"].to_numpy(), x=rd_data["x"].to_numpy(),
               h=0.1, rho=1, vce="hc3"))

print(
    rdhte(
        y="y",
        x="x",
        h=0.116,
        vce="hc3",
        subset=rd_data["w_left"] == 1,
        data=rd_data,
    ).summary()
)

print(
    rdrobust(
        y=rd_data.loc[rd_data["w_left"] == 1, "y"].to_numpy(),
        x=rd_data.loc[rd_data["w_left"] == 1, "x"].to_numpy(),
        h=0.116,
        rho=1,
        vce="hc3",
    )
)
