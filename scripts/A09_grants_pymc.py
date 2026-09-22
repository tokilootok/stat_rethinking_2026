#!/usr/bin/env python3
"""PyMC implementation of the NWOGrants saturated binomial model."""

from __future__ import annotations

import argparse
from pathlib import Path

import arviz as az
import numpy as np
import pandas as pd
import pymc as pm


def build_model(frame: pd.DataFrame) -> pm.Model:
    gender = frame["gender_id"].to_numpy(dtype=int) - 1
    discipline = frame["discipline_id"].to_numpy(dtype=int) - 1
    applications = frame["applications"].to_numpy(dtype=int)
    awards = frame["awards"].to_numpy(dtype=int)
    n_disciplines = int(frame["discipline_id"].max())

    pooled = frame.groupby("discipline_id", sort=True)["applications"].sum().to_numpy()
    w_pool = pooled / pooled.sum()
    female_apps = (
        frame.loc[frame.gender_id == 1]
        .groupby("discipline_id", sort=True)["applications"]
        .sum()
        .to_numpy()
    )
    male_apps = (
        frame.loc[frame.gender_id == 2]
        .groupby("discipline_id", sort=True)["applications"]
        .sum()
        .to_numpy()
    )
    w_female = female_apps / female_apps.sum()
    w_male = male_apps / male_apps.sum()

    coords = {
        "cell": np.arange(len(frame)),
        "gender": ["female", "male"],
        "discipline": np.arange(n_disciplines),
    }
    with pm.Model(coords=coords) as model:
        g = pm.Data("gender_id", gender, dims="cell")
        d = pm.Data("discipline_id", discipline, dims="cell")
        n = pm.Data("applications", applications, dims="cell")

        logit_award = pm.Normal("logit_award", 0.0, 1.0, dims=("gender", "discipline"))
        p_award = pm.Deterministic("p_award", pm.math.sigmoid(logit_award), dims=("gender", "discipline"))
        pm.Binomial("awards", n=n, p=p_award[g, d], observed=awards, dims="cell")

        total_female = pm.Deterministic("total_female", pm.math.dot(p_award[0], w_female))
        total_male = pm.Deterministic("total_male", pm.math.dot(p_award[1], w_male))
        pm.Deterministic("total_risk_difference", total_female - total_male)

        standardized_female = pm.Deterministic("standardized_female", pm.math.dot(p_award[0], w_pool))
        standardized_male = pm.Deterministic("standardized_male", pm.math.dot(p_award[1], w_pool))
        pm.Deterministic(
            "standardized_direct_risk_difference",
            standardized_female - standardized_male,
        )
    return model


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", default="artifacts/A09/nwo_grants.csv")
    parser.add_argument("--output", default="artifacts/A09")
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    frame = pd.read_csv(args.data)
    model = build_model(frame)
    draws = 250 if args.quick else 1000
    tune = 250 if args.quick else 1000
    chains = 2 if args.quick else 4

    with model:
        prior = pm.sample_prior_predictive(draws=draws, random_seed=202609)
        idata = pm.sample(
            draws=draws,
            tune=tune,
            chains=chains,
            cores=chains,
            target_accept=0.95,
            random_seed=202609,
            return_inferencedata=True,
            progressbar=not args.quick,
        )
        idata = pm.sample_posterior_predictive(
            idata, random_seed=202609, extend_inferencedata=True
        )

    key = [
        "total_female",
        "total_male",
        "total_risk_difference",
        "standardized_female",
        "standardized_male",
        "standardized_direct_risk_difference",
    ]
    summary = az.summary(idata, var_names=key, round_to=4)
    print(summary)
    if not args.quick:
        if float(summary["r_hat"].max()) >= 1.01:
            raise RuntimeError("At least one estimand has R-hat >= 1.01")
        if int(idata.sample_stats["diverging"].sum()):
            raise RuntimeError("Detected divergent transitions")

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    prior.to_netcdf(output / "pymc_grants_prior.nc")
    summary.to_csv(output / "pymc_summary.csv")
    idata.to_netcdf(output / "pymc_grants.nc")


if __name__ == "__main__":
    main()

