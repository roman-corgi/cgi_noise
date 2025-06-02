"""
This script simulates an astronomical observation scenario, specifically for exoplanet imaging.
It calculates various parameters related to the target star and planet, instrument performance,
and estimates the time required to achieve a desired signal-to-noise ratio (SNR) for detecting
the exoplanet.

The script performs the following major steps:
1. Loads a scenario configuration from a YAML file in data/Scenarios/.
2. Defines astronomical and instrument constants.
3. Sets up the target host star and exoplanet parameters.
4. Loads required CSV data files for coronagraph, detector, throughput, etc.
5. Calculates the planet's working angle.
6. Determines contrast stability parameters.
7. Sets up coronagraph and focal plane parameters.
8. Calculates star and planet flux.
9. Estimates background zodi and speckle flux.
10. Computes core photon rates (planet, speckle, zodi).
11. Determines optimal frame time and differential quantum efficiency (dQE).
12. Accounts for Reference Differential Imaging (RDI) noise penalties.
13. Calculates detector noise rates.
14. Computes variance rates and residual speckle rates.
15. Calculates the time to achieve the desired SNR.

Run this script from the project root with:
    python -m cgi_noise.cgi_noise_main --scenario <scenario_file.yml>
Or in Spyder with %runfile, ensuring the working directory is D:\Repos\cgi_noise3.
"""
import os
import sys
import argparse
import math
from datetime import datetime
from dataclasses import dataclass
import numpy as np
import yaml
from pathlib import Path
from prettytable import PrettyTable

# Add project root to sys.path for Spyder %runfile execution
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from cgi_noise import check_data_path
from cgi_noise.unitsConstants import mas, nm, hour, second, jupiterRadius, AU, sunAbsMag, mm
from cgi_noise.cginoiselib import (
    Target, open_folder, loadCSVs, getScenFileNames, workingAnglePars,
    contrastStabilityPars, coronagraphParameters, getFocalPlaneAttributes,
    getSpectra, getStrayLightfromfile, compute_throughputs, rdi_noise_penalty,
    compute_frame_time_and_dqe, detector_noise_rates, noiseVarianceRates,
    compute_tsnr
)

def main():
    """Main function to compute tSNR for a specified scenario."""
    # Parse command-line argument for scenario file
    parser = argparse.ArgumentParser(description="Compute tSNR for an exoplanet imaging scenario.")
    parser.add_argument(
        "--scenario",
        default="OPT_IMG_NFB1_HLC.yml",
        help="Name of the YAML scenario file in data/Scenarios/ (default: OPT_IMG_NFB1_HLC.yml)"
    )
    args = parser.parse_args()

    # Ensure data/ folder is accessible
    data_path = check_data_path()

    current_dir = os.getcwd()
    print(f"Working directory: {current_dir}")
    current_datetime = datetime.now()
    print(f"Run started at: {current_datetime}")

    # Step 1: Load scenario configuration
    scenario_filename = args.scenario
    scenFolder = data_path / "Scenarios"
    try:
        with open(scenFolder / scenario_filename, "r") as file:
            config = yaml.safe_load(file)
    except FileNotFoundError:
        print(f"Error: {scenario_filename} not found in {scenFolder}! "
              "Ensure the 'data/' folder from https://github.com/roman-corgi/cgi_noise is present.")
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"Error parsing YAML: {e}")
        sys.exit(1)

    ObservationType = config['DataSpecification']['ObservationCase']
    print(f'\nObservation type: {ObservationType}')

    # Step 2: Define constants
    DPM = config['instrument']['Diam']
    lam = config['instrument']['wavelength']
    lamD = lam / DPM
    intTimeDutyFactor = config['instrument']['dutyFactor']
    opMode = config['instrument']['OpMode']
    bandWidth = config['instrument']['bandwidth']

    print(f"Central wavelength: {lam / nm} nm, with {bandWidth*100:.0f}% BW\n")
    print(f"Lambda/D: {lamD / mas:.3f} mas")

    # Step 3: Define host star and planet
    target = Target(
        v_mag=5.0,
        dist_pc=10.0,
        specType='g0v',
        phaseAng_deg=65,
        sma_AU=4.1536,
        radius_Rjup=5.6211,
        geomAlb_ag=0.44765,
        exoZodi=1,
    )

    # User-defined inputs
    allocatedTime = 100 * hour
    monthsAtL2 = 21
    frameTime = 10.0 * second
    isPhotonCounting = True

    usableTinteg = intTimeDutyFactor * allocatedTime

    sep_mas = target.phaseAng_to_sep(target.sma_AU, target.dist_pc, target.phaseAng_deg)
    target.albedo = target.albedo_from_geomAlbedo(target.phaseAng_deg, target.geomAlb_ag)

    print(f"Separation: {sep_mas:.0f} mas")
    print(f"Albedo: {target.albedo:.3f}")

    # Step 4: Load required CSV files
    filenameList = getScenFileNames(config, data_path=data_path)
    CG_Data, QE_Data, DET_CBE_Data, STRAY_FRN_Data, THPT_Data, CAL_Data, CS_Data = loadCSVs(filenameList)
    CS_Type = config['DataSpecification']['CS_Type']

    # Step 5: Planet working angle
    IWA, OWA = workingAnglePars(CG_Data, CS_Data)
    planetWA = sep_mas * mas / lamD

    table = PrettyTable()
    table.field_names = ['planet WA', 'phase, deg', 'dist, pc', 'sma, AU', 'sep, mas', 'lam/D, mas', "IWA", "OWA"]
    table.add_row([f'{planetWA:.2f}', f'{target.phaseAng_deg:.2f}', f"{target.dist_pc:.2f}",
                   f'{target.sma_AU:.2f}', f'{sep_mas:.2f}', f'{lamD/mas:.2f}', f'{IWA:.2f}', f'{OWA:.2f}'])
    print(table)

    tolerance = 0.05
    if (IWA - tolerance) <= planetWA <= IWA:
        planetWA = IWA
    elif OWA <= planetWA <= (OWA + tolerance):
        planetWA = OWA
    elif planetWA < (IWA - tolerance) or planetWA > (OWA + tolerance):
        raise ValueError(f"Planet WA={planetWA:.1f} while IWA = {IWA:.1f} and OWA = {OWA:.1f} lam/D.")

    print(f"Planet Working Angle: {planetWA:.2f} λ/D")

    # Step 6: Contrast stability parameters
    selDeltaC, rawContrast, SystematicCont, initStatRawContrast, \
        rawContrast, IntContStab, ExtContStab = contrastStabilityPars(CS_Type, planetWA, CS_Data)

    print(f"Raw Contrast: {rawContrast:.3e}")
    print(f"Selected Delta Contrast: {selDeltaC:.3e}")

    # Step 7: Coronagraph slice parameters
    cg = coronagraphParameters(CG_Data.df, config, planetWA, DPM)

    # Step 8: Focal plane setup
    f_SR, CritLam, detPixSize_m, mpix, pixPlateSc = getFocalPlaneAttributes(
        opMode, config, DET_CBE_Data, lam, bandWidth, DPM, cg.CGdesignWL, cg.omegaPSF
    )

    # Step 9: Star flux
    inBandFlux0_sum, inBandZeroMagFlux, starFlux = getSpectra(target, lam, bandWidth)
    print(f"Star Flux = {starFlux:.3e} ph/s/m^2")
    TimeonRefStar_tRef_per_tTar = 0.25

    # Step 10: Planet flux
    fluxRatio = target.albedo * (target.radius_Rjup * jupiterRadius / (target.sma_AU * AU)) ** 2
    planetFlux = fluxRatio * starFlux
    print(f"Planet Flux Ratio = {fluxRatio:.2e}\nPlanet Flux = {planetFlux:.3f} ph/s/m^2")

    # Step 11: Background zodi and speckle flux
    magLocalZodi = config['instrument']['LocZodi_magas2']
    magExoZodi_1AU = config['instrument']['ExoZodi_magas2']
    absMag = target.v_mag - 5 * math.log10(target.dist_pc / 10)

    locZodiAngFlux = inBandZeroMagFlux * 10 ** (-0.4 * magLocalZodi)
    exoZodiAngFlux = inBandZeroMagFlux * 10 ** (-0.4 * (absMag - sunAbsMag + magExoZodi_1AU)) / target.sma_AU**2 * target.exoZodi

    exoZodiDistrib = "uniform"  # Options: "lumpy", "uniform", "falloff"

    thput, throughput_rates = compute_throughputs(THPT_Data, cg, exoZodiDistrib)

    planetThroughput = throughput_rates["planet"]
    speckleThroughput = throughput_rates["speckle"]
    locZodiThroughput = throughput_rates["local_zodi"]
    exoZodiThroughput = throughput_rates["exo_zodi"]

    Acol = (np.pi / 4.0) * DPM**2

    scenario = 'Threshold IMG NF B1'
    perfLevel = 'CBE'
    stray_ph_s_mm2 = getStrayLightfromfile(scenario, perfLevel, STRAY_FRN_Data)
    stray_ph_s_pix = stray_ph_s_mm2 * (1/mm**2) * detPixSize_m**2

    # Step 12: RDI noise penalties
    RefStarSpecType = 'a0v'
    RefStarVmag = 2.26

    rdi_penalty = rdi_noise_penalty(
        target, inBandFlux0_sum, starFlux, TimeonRefStar_tRef_per_tTar, RefStarSpecType, RefStarVmag
    )
    k_sp = rdi_penalty["k_sp"]
    k_det = rdi_penalty["k_det"]
    k_lzo = rdi_penalty["k_lzo"]
    k_ezo = rdi_penalty["k_ezo"]

    @dataclass
    class corePhotonRates:
        planet: float = planetFlux * planetThroughput * Acol
        speckle: float = starFlux * rawContrast * cg.PSFpeakI * cg.CGintmpix * speckleThroughput * Acol
        locZodi: float = locZodiAngFlux * cg.omegaPSF * locZodiThroughput * Acol
        exoZodi: float = exoZodiAngFlux * cg.omegaPSF * exoZodiThroughput * Acol
        straylt: float = stray_ph_s_pix * mpix
        total: float = planet + speckle + locZodi + exoZodi + straylt
    cphrate = corePhotonRates()

    # Step 13: Frame time and dQE calculation
    desiredRate = 0.1  # e-/pix/frame
    tfmin = 3  # min frame time (s)
    tfmax = 100  # max frame time (s)

    ENF, effReadnoise, frameTime, dQE, QE_img = compute_frame_time_and_dqe(
        desiredRate, tfmin, tfmax, isPhotonCounting, QE_Data, DET_CBE_Data, lam, mpix, cphrate.total
    )
    print(f"Calculated Frame Time: {frameTime:.2f} s")
    print(f"QE in the image area: {QE_img:.3f}")
    print(f"Detected Quantum Efficiency (dQE): {dQE:.3f}")
    print(f"Excess Noise Factor (ENF): {ENF:.2f}")
    print(f"Core fraction used in the SNR region for mode {ObservationType}: f_SR: {f_SR:.3f}")

    # Step 14: Detector noise calculation
    detNoiseRate = detector_noise_rates(DET_CBE_Data, monthsAtL2, frameTime, mpix, isPhotonCounting)

    # Step 15: Variance and SNR calculation
    k_pp = config['instrument']['pp_Factor_CBE']
    nvRatesCore, residSpecRate = noiseVarianceRates(
        cphrate=cphrate, QE=QE_img, dQE=dQE, ENF=ENF, detNoiseRate=detNoiseRate,
        k_sp=k_sp, k_det=k_det, k_lzo=k_lzo, k_ezo=k_ezo, f_SR=f_SR,
        starFlux=starFlux, selDeltaC=selDeltaC, k_pp=k_pp, cg=cg,
        speckleThroughput=speckleThroughput, Acol=Acol
    )

    SNRdesired = 5.0
    planetSignalRate = f_SR * cphrate.planet * dQE
    timeToSNR, criticalSNR = compute_tsnr(SNRdesired, planetSignalRate, nvRatesCore, residSpecRate)

    print(f"\nTarget SNR = {SNRdesired:.1f} \nCritical SNR = {criticalSNR:.2f}")
    print(f"Time to SNR = {timeToSNR:.1f} seconds or {timeToSNR/hour:.3f} hours")

    if timeToSNR > usableTinteg:
        print(f"Warning: Time to SNR ({timeToSNR/hour:.2f} hrs) exceeds usable integration time ({usableTinteg/hour:.2f} hrs).")
    elif timeToSNR <= 0:
        print(f"Warning: Calculated Time to SNR is not positive ({timeToSNR/hour:.2f} hrs). Check input parameters.")

    print(f"\nRun completed at: {datetime.now()}")

if __name__ == "__main__":
    main()