from cgi_noise.main_core import run_pipeline
from pathlib import Path
import sys
import yaml
from datetime import datetime
import os


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

    run_pipeline(config, DATA_DIR, obs_params["target_params"], obs_params["snr"])


if __name__ == "__main__":
    print(f"Run started at: {datetime.now()}")

    obs_params = {
        "scenario": 'CON_IMG_WFB4_SPC.yml', #'OPT_IMG_NFB1_HLC.yml', #   "OPT_IMG_NFB1_HLC.yml", # 'CON_IMG_NFB1_HLC.yml' # CON_SPEC_NFB3_SPC.yml
        "target_params": {
            "v_mag": 5.0,
            "dist_pc": 10.0,
            "specType": "g0v",
            "phaseAng_deg": 65,
            "sma_AU": 8.1536,
            "radius_Rjup": 5.6211,
            "geomAlb_ag": 0.44765,
            "exoZodi": 1,
        },
        "snr": 5.0,
    }

    run_snr_scenario(obs_params)
    print(f"Run completed at: {datetime.now()}")
