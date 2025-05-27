# -*- coding: utf-8 -*-
"""
Created on Mon May 26 21:22:12 2025

@author: bngen
"""

# main_driver.py - Clean version of Main.py

import os
import math
import pandas as pd
from datetime import datetime
from pathlib import Path

from loadXLcol import loadXLcol
from loadCSVrow import loadCSVrow
import unitsConstants as uc
import func_library as fl


def main():
    current_dir = os.getcwd()
    print(f"Working directory: {current_dir}")
    current_datetime = datetime.now()
    print(f"Run started at: {current_datetime}")

    # === Scenario Selection ===
    scenario_filename = 'SCEN_IMG_NFOV_B1_HLC.xlsx'
    scenFolder = fl.open_folder("EBcsvData", "Scenarios")
    scenarioDF = loadXLcol(scenFolder[scenario_filename], 30).df
    scenario = scenarioDF.at['Scenario', 'Latest']

    # === Constants ===
    DPM = 2.363 * uc.meter
    lam = scenarioDF.at['CenterLambda_nm', 'Latest'] * uc.nm
    lamD = lam / DPM
    intTimeDutyFactor = scenarioDF.at['DutyFactor_CBE', 'Latest']
    usableTinteg = 100 * intTimeDutyFactor * uc.hour

    print(f"Wavelength: {lam / uc.nm} nm")
    print(f"Lambda/D: {lamD / uc.mas:.3f} mas")

    # === Define Host Star and Planet ===
    target = fl.Target(
        v_mag=5.05,
        dist_pc=13.8,
        specType='g0v',
        phaseAng_deg=65,
        sma_AU=5,
        radius_Rjup=1,
        geomAlb_ag=0.3,
        exoZodi=1,
    )

    sep_mas = target.phaseAng_to_sep(target.sma_AU, target.dist_pc, target.phaseAng_deg)
    target.albedo = target.fluxRatio_SMA_rad_to_albedo(
        fl.Target.deltaMag_to_fluxRatio(
            fl.Target.fluxRatio_to_deltaMag(5e-9)
        ), target.sma_AU, target.radius_Rjup)

    print(f"Separation: {sep_mas:.0f} mas")
    print(f"Albedo: {target.albedo:.3f}")

    # === Load Required CSV Files ===
    filenameList = fl.getScenFileNames_DRM(scenarioDF)
    CG_Data, QE_Data, DET_CBE_Data, STRAY_FRN_Data, THPT_Data, CAL_Data, CS_Data = fl.loadCSVs(filenameList)

    # === Planet Working Angle ===
    IWA, OWA = fl.workingAnglePars(CG_Data, CS_Data)
    planetWA = fl.Target.separation_to_planetWA(sep_mas, lamD, IWA, OWA)
    print(f"Planet Working Angle: {planetWA:.2f} λ/D")

    # === Throughput Parameters ===
    CSprefix = 'MCBE_'
    selDeltaC, selContrast, SystematicCont, initStatRawContrast, \
        rawContrast, IntContStab, ExtContStab = fl.contrastStabilityPars(CSprefix, planetWA, CS_Data)

    print(f"Raw Contrast: {rawContrast:.3e}")
    print(f"Selected Delta Contrast: {selDeltaC:.3e}")

    # === Placeholder Output ===
    print("Setup complete. Ready to proceed with detailed calculations.")


if __name__ == '__main__':
    main()
