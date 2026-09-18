#!/usr/bin/env python3
"""PyMC implementation of the same identified wine-judge IRT model as Stan."""

from __future__ import annotations

import argparse
from pathlib import Path

import arviz as az
import numpy as np
import pandas as pd
import pymc as pm


def constant_by_id(frame: pd.DataFrame, id_col: str, value_col: str) -> np.ndarray:
    grouped = frame.groupby(id_col, sort=True)[value_col]
    if not (grouped.nunique() == 1).all():
        raise ValueError(f"{value_col} must be constant within {id_col}")
    return grouped.first().to_numpy(dtype=float)


def build_model(frame: pd.DataFrame) -> pm.Model:
    score_z = frame["score_z"].to_numpy(dtype=float)
    wine_id = frame["wine_id"].to_numpy(dtype=int) - 1
    judge_id = frame["judge_id"].to_numpy(dtype=int) - 1
    wine_origin_c = constant_by_id(frame, "wine_id", "wine_origin_c")
    judge_country_c = constant_by_id(frame, "judge_id", "judge_country_c")

    n_wines = len(wine_origin_c)
    n_judges = len(judge_country_c)
    coords = {
        "observation": np.arange(len(frame)),
        "wine": np.arange(n_wines),
        "judge": np.arange(n_judges),
    }

    with pm.Model(coords=coords) as model:
        w = pm.Data("wine_id", wine_id, dims="observation")
        j = pm.Data("judge_id", judge_id, dims="observation")
        x_wine = pm.Data("wine_origin_c", wine_origin_c, dims="wine")
        x_judge = pm.Data("judge_country_c", judge_country_c, dims="judge")

        delta_origin = pm.Normal("delta_origin", 0.0, 0.5)
        quality = pm.Normal("quality", delta_origin * x_wine, 1.0, dims="wine")

        alpha_harshness = pm.Normal("alpha_harshness", 0.0, 0.5)
        delta_harshness = pm.Normal("delta_harshness", 0.0, 0.5)
        sigma_harshness = pm.Exponential("sigma_harshness", 2.0)
        z_harshness = pm.Normal("z_harshness", 0.0, 1.0, dims="judge")
        harshness = pm.Deterministic(
            "harshness",
            alpha_harshness
            + delta_harshness * x_judge
            + sigma_harshness * z_harshness,
            dims="judge",
        )

        alpha_log_discrimination = pm.Normal(
            "alpha_log_discrimination", 0.0, 0.35
        )
        delta_log_discrimination = pm.Normal(
            "delta_log_discrimination", 0.0, 0.35
        )
        sigma_log_discrimination = pm.Exponential(
            "sigma_log_discrimination", 2.0
        )
        z_log_discrimination = pm.Normal(
            "z_log_discrimination", 0.0, 1.0, dims="judge"
        )
        log_discrimination = pm.Deterministic(
            "log_discrimination",
            alpha_log_discrimination
            + delta_log_discrimination * x_judge
            + sigma_log_discrimination * z_log_discrimination,
            dims="judge",
        )
        discrimination = pm.Deterministic(
            "discrimination", pm.math.exp(log_discrimination), dims="judge"
        )
        pm.Deterministic(
            "discrimination_ratio_us_fr", pm.math.exp(delta_log_discrimination)
        )

        sigma_score = pm.Exponential("sigma_score", 1.0)
        mu = pm.Deterministic(
            "mu",
            discrimination[j] * (quality[w] - harshness[j]),
            dims="observation",
        )
        pm.Normal("score_z", mu, sigma_score, observed=score_z, dims="observation")

    return model


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", default="artifacts/A08/wines2012.csv")
    parser.add_argument("--output", default="artifacts/A08")
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--prior-only", action="store_true")
    args = parser.parse_args()

    frame = pd.read_csv(args.data)
    model = build_model(frame)
    draws = 250 if args.quick else 1000
    tune = 250 if args.quick else 1000
    chains = 2 if args.quick else 4

    with model:
        prior = pm.sample_prior_predictive(draws=draws, random_seed=202608)
        if args.prior_only:
            print(prior.prior["score_z"].shape)
            return
        idata = pm.sample(
            draws=draws,
            tune=tune,
            chains=chains,
            cores=chains,
            target_accept=0.99,
            random_seed=202608,
            return_inferencedata=True,
            progressbar=not args.quick,
        )
        idata = pm.sample_posterior_predictive(
            idata, random_seed=202608, extend_inferencedata=True
        )

    key = [
        "delta_harshness",
        "delta_log_discrimination",
        "discrimination_ratio_us_fr",
        "delta_origin",
        "sigma_score",
    ]
    summary = az.summary(idata, var_names=key, round_to=3)
    print(summary)

    if not args.quick:
        if float(summary["r_hat"].max()) >= 1.01:
            raise RuntimeError("At least one key parameter has R-hat >= 1.01")
        divergences = int(idata.sample_stats["diverging"].sum())
        if divergences:
            raise RuntimeError(f"Detected {divergences} divergent transitions")

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output / "pymc_summary.csv")
    idata.to_netcdf(output / "pymc_irt.nc")


if __name__ == "__main__":
    main()

