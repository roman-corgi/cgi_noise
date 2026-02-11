from cgi_noise.tsnr_core import tsnr_pipeline
from pathlib import Path
import sys
import yaml
from datetime import datetime
import os
import argparse
import cgi_noise.cginoiselib as fl


def run_snr_scenario(obs_params):
    DATA_DIR = Path(os.environ["CGI_NOISE_DATA_DIR"])
    SCEN_DIR = DATA_DIR / "Scenarios"

    scenario_filename = obs_params["scenario"]
    scenario_path = SCEN_DIR / scenario_filename
    print(f"Looking for config at: {scenario_path.resolve()}")

    try:
        with open(SCEN_DIR / scenario_filename, "r") as file:
            config = yaml.safe_load(file)
    except FileNotFoundError:
        print(f"Error: {scenario_path.resolve()} not found")
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"YAML error: {e}")
        sys.exit(1)

    planetSignalRate, nvRatesCore, residSpecSdevRate = tsnr_pipeline(
        config, DATA_DIR, obs_params["target_params"], obs_params["verbose"]
    )
    timeToSNR, criticalSNR = fl.compute_tsnr(
        obs_params["snr"], planetSignalRate, nvRatesCore, residSpecSdevRate
    )
    print(
        f"Integration time to SNR {obs_params['snr']:.1f}: {timeToSNR:.1f} sec, Critical SNR = {criticalSNR:.2f}"
    )


def main():
    scenarios = [
        "OPT_IMG_NFB1_HLC.yml",
        "CON_IMG_NFB1_HLC.yml",
        "OPT_SPEC_NFB3_SPC.yml",
        "CON_SPEC_NFB3_SPC.yml",
        "OPT_IMG_WFB4_SPC.yml",
        "CON_IMG_WFB4_SPC.yml",
        "OPT_IMG_WFB1_SPC.yml",
        "CON_IMG_WFB1_SPC.yml",
    ]

    parser = argparse.ArgumentParser(
        description="Run cgi_noise integration time calculation."
    )
    parser.add_argument(
        "-s",
        "--scenario",
        nargs="?",
        type=str,
        help="Scenario Name (string).",
        default=scenarios[2],
        choices=scenarios,
    )
    parser.add_argument(
        "--sma", nargs="?", type=float, default=4.1536, help="Planet sma in AU (float)."
    )
    parser.add_argument(
        "--radius",
        nargs="?",
        type=float,
        default=5.6211,
        help="Planet radius in R_jupiter (float).",
    )

    args = parser.parse_args()
    scenario = args.scenario
    assert scenario in scenarios, f"Scenario must be one of: {','.join(scenarios)}."
    sma = args.sma
    radius = args.radius

    print(f"Run started at: {datetime.now()}")

    obs_params = {
        "scenario": scenario,
        "target_params": {
            "v_mag": 5.0,
            "dist_pc": 10.0,
            "specType": "g0v",
            "phaseAng_deg": 65,
            "sma_AU": sma,
            "radius_Rjup": radius,
            "geomAlb_ag": 0.44765,
            "exoZodi": 1,
        },
        "snr": 5.0,
        "verbose": True,
    }

    run_snr_scenario(obs_params)
    print(f"Run completed at: {datetime.now()}")


if __name__ == "__main__":
    main()
