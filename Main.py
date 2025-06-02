from cgi_noise.main_core import run_pipeline
import argparse
from pathlib import Path
import sys
import yaml
from datetime import datetime

def run_snr_scenario(obs_params):
    ROOT_DIR = Path(__file__).resolve().parent
    DATA_DIR = ROOT_DIR / "cgi_noise" / "data"
    SCEN_DIR = DATA_DIR / "Scenarios"

    scenario_filename = obs_params["scenario"]
    print(f"Looking for config at: {(SCEN_DIR / scenario_filename).resolve()}")

    try:
        with open(SCEN_DIR / scenario_filename, "r") as file:
            config = yaml.safe_load(file)
    except FileNotFoundError:
        print(f"Error: {scenario_filename} not found in Scenarios/")
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"YAML error: {e}")
        sys.exit(1)

    run_pipeline(config, DATA_DIR, obs_params["target_params"], obs_params["snr"])


if __name__ == "__main__":
    print(f"Run started at: {datetime.now()}")

    obs_params = {
        "scenario": "OPT_IMG_NFB1_HLC.yml",
        "target_params": {
            'v_mag': 5.0,
            'dist_pc': 10.0,
            'specType': 'g0v',
            'phaseAng_deg': 65,
            'sma_AU': 4.1536,
            'radius_Rjup': 5.6211,
            'geomAlb_ag': 0.44765,
            'exoZodi': 1,
        },
        "snr": 5.0
    }

    run_snr_scenario(obs_params)
    print(f"Run completed at: {datetime.now()}")
