"""
This script simulates an astronomical observation scenario, specifically for exoplanet imaging. 
It calculates various parameters related to the target star and planet, instrument performance, 
and ultimately estimates the time required to achieve a desired signal-to-noise ratio (SNR) 
for detecting the exoplanet.

The script performs the following major steps:
1.  Loads a scenario configuration from a YAML file.
2.  Defines astronomical and instrument constants.
3.  Sets up the target host star and exoplanet parameters.
4.  Loads required CSV data files for coronagraph, detector, throughput, etc.
5.  Calculates the planet's working angle.
6.  Determines contrast stability parameters.
7.  Sets up coronagraph and focal plane parameters.
8.  Calculates star and planet flux.
9.  Estimates background zodi and speckle flux.
10. Computes core photon rates (planet, speckle, zodi).
11. Determines optimal frame time and differential quantum efficiency (dQE).
12. Accounts for Reference Differential Imaging (RDI) noise penalties.
13. Calculates detector noise rates.
14. Computes variance rates and residual speckle rates.
15. Calculates the time to achieve the desired SNR.

Uses external libraries: 'unitsConstants' for physical unit conversions and 
'cginoiselib' for specialized astronomical and instrument calculations.
"""
# from cgi_noise import cginoiselib as fl
from cgi_noise.main_core import run_pipeline
# import cgi_noise.unitsConstants as uc
import argparse
from pathlib import Path
import sys
import yaml
from datetime import datetime


def run_snr_scenario(config_filename: str):
    ROOT_DIR = Path(__file__).resolve().parent
    DATA_DIR = ROOT_DIR / "cgi_noise" / "data"
    SCEN_DIR = DATA_DIR / "Scenarios"
    print(f"Looking for config at: {(SCEN_DIR / config_filename).resolve()}")
    # Load config
    try:
        with open(SCEN_DIR / config_filename, "r") as file:
            config = yaml.safe_load(file)
    except FileNotFoundError:
        print(f"Error: {config_filename} not found in Scenarios/")
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"YAML error: {e}")
        sys.exit(1)

    # Run the full simulation pipeline
    run_pipeline(config, DATA_DIR)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run SNR time estimation scenario.")
    parser.add_argument(
        "--config", type=str, default="OPT_IMG_NFB1_HLC.yml",
        help="YAML configuration file in data/Scenarios/"
    )
    args = parser.parse_args()
    print(f"Run started at: {datetime.now()}")
    run_snr_scenario(args.config)
    print(f"Run completed at: {datetime.now()}")
