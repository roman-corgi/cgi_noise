from cgi_noise.sens_core import sens_pipeline
from pathlib import Path
import sys
import yaml
from datetime import datetime
import os
import argparse
import cgi_noise.cginoiselib as fl
import cgi_noise.unitsConstants as uc
import numpy as np
import matplotlib.pyplot as plt

def run_sensitivity_scenario(obs_params):
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
        
    filenameList = fl.getScenFileNames(config, DATA_DIR)
    CG_Data, QE_Data, DET_CBE_Data, STRAY_FRN_Data, THPT_Data, CAL_Data, CS_Data = fl.loadCSVs(filenameList)

    SNR = 5
    Thrs = 100000
    WAset, sep, Sensitivity, KC, tauPk, intMpix, tauC, coreArea, dC = sens_pipeline(config, DATA_DIR, obs_params["target_params"], False, Thrs*uc.hour, SNR   )
        

    plt.figure(figsize=(8, 6))
    plt.plot(sep/uc.mas, Sensitivity/uc.ppb, 'b-', marker='o', label=f'Sensitivity  ppb')
    plt.xlabel('Separation (mas)')
    plt.ylabel(f'{SNR:.1f}-sigma Sensitivity, ppb')
    plt.title(f'SNR={SNR:.1f} {Thrs:.1f} Hrs {scenario_filename}')
    plt.grid(True)
    plt.legend()
 
    mdC = np.mean(dC)
    mSen = np.mean(Sensitivity)
    mKC  = np.mean(KC)
    mtp = np.mean(tauPk)
    mta = np.mean(tauPk*coreArea)
    mitc = np.mean(1/tauC)
    plt.figure(figsize=(8, 6))
    plt.rcParams['mathtext.fontset'] = 'stix'  # Alternatives: 'cm', 'dejavusans'
    plt.rcParams['text.usetex'] = False  # Ensure external LaTeX is disabled
    plt.plot(WAset, Sensitivity/mSen, 'k-', marker='o', label=f'Sensitivity / {mSen/uc.ppb:.2f} ppb')
    plt.plot(WAset, dC/mdC, 'y-', marker='d', label=f'delta C (CS) / {mdC/uc.ppb:.2f} ppb')
    plt.plot(WAset, KC/mKC, 'r-.', label=f'KappaC / {mKC:.2f}')
    plt.plot(WAset, tauPk/mtp, 'b:', label=f'tau_pk / {mtp:.2e}')
    plt.plot(WAset, tauPk*coreArea/mta, 'c--', label=f'tau_pk * Omega_core / {mta:.2e}')
    plt.plot(WAset, (1/tauC)/mitc , 'm--', label=f'(1/tauCore) / {mitc:.2e}')
    
    plt.xlabel('Working Angle, lam/D')
    # plt.ylabel(f'{SNR:.1f}-sigma Sensitivity, ppb')
    plt.title(f'SNR={SNR:.1f}    {Thrs:.1e} Hrs    {scenario_filename}')
    plt.grid(True)
    plt.legend()
 

 

def main():
    scenarios = [
        "OPT_IMG_NFB1_HLC.yml",
        "CON_IMG_NFB1_HLC.yml",
        "OPT_SPEC_NFB3_SPC.yml",
        "CON_SPEC_NFB3_SPC.yml",
        "OPT_IMG_WFB4_SPC.yml",
        "CON_IMG_WFB4_SPC.yml",
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
        "verbose": True
    }

    run_sensitivity_scenario(obs_params)
    print(f"Run completed at: {datetime.now()}")


if __name__ == "__main__":
    main()
