# main_driver.py - Converted to a script version (no main function)

import os
import math
from datetime import datetime
from loadXLcol import loadXLcol
import unitsConstants as uc
import cginoiselib as fl
from dataclasses import dataclass, asdict
import numpy as np
import yaml
import sys

current_dir = os.getcwd()
print(f"Working directory: {current_dir}")
current_datetime = datetime.now()
print(f"Run started at: {current_datetime}")

# === Scenario Selection ===
scenario_filename = 'SCEN_IMG_NFOV_B1_HLC.yaml'
scenFolder = fl.open_folder("EBcsvData", "Scenarios")
# scenarioDF = loadXLcol(scenFolder[scenario_filename], 30).df
# scenario = scenarioDF.at['Scenario', 'Latest']

# Load YAML config
try:
    with open(scenFolder[scenario_filename], "r") as file:
        config = yaml.safe_load(file)
except FileNotFoundError:
    print("Error: config.yaml not found!")
    sys.exit(1)
except yaml.YAMLError as e:
    print(f"Error parsing YAML: {e}")
    sys.exit(1)
 
# # Extract parameters from YAML
# target = config.get("target", {})
# targetName = target.get("name", None)  # Target title/description
# stellar_type = target.get("stellarType", None)  # Get stellar type
# absmag = float(target.get("absMag", None))
# exoSolar = float(target.get("exoSolar", None))
# exoZodiSurfBright = float(target.get("exoZodiMagas2"))
# locZodiSurfBright = float(target.get("locZodiMagas2"))


# === Constants ===
DPM = 2.363 * uc.meter
lam = config['instrument']['wavelength']
lamD = lam / DPM
intTimeDutyFactor = config['instrument']['dutyFactor'] # scenarioDF.at['DutyFactor_CBE', 'Latest']

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

# === User-defined inputs ===
allocatedTime = 100 * uc.hour
monthsAtL2 = 21
frameTime = 10.0 * uc.second
isPhotonCounting = True  # or False, depending on the mode



usableTinteg = intTimeDutyFactor * allocatedTime

sep_mas = target.phaseAng_to_sep(target.sma_AU, target.dist_pc, target.phaseAng_deg)
target.albedo = target.fluxRatio_SMA_rad_to_albedo(
    fl.Target.deltaMag_to_fluxRatio(
        fl.Target.fluxRatio_to_deltaMag(5e-9)
    ), target.sma_AU, target.radius_Rjup)

print(f"Separation: {sep_mas:.0f} mas")
print(f"Albedo: {target.albedo:.3f}")

# === Load Required CSV Files ===
filenameList = fl.getScenFileNames_YML(config)
CG_Data, QE_Data, DET_CBE_Data, STRAY_FRN_Data, THPT_Data, CAL_Data, CS_Data = fl.loadCSVs(filenameList)

# === Planet Working Angle ===
IWA, OWA = fl.workingAnglePars(CG_Data, CS_Data)
planetWA = sep_mas * uc.mas / lamD
tolerance = 0.05
if (IWA - tolerance) <= planetWA <= IWA:
    planetWA = IWA
elif OWA <= planetWA <= (OWA + tolerance):
    planetWA = OWA
elif planetWA < (IWA - tolerance) or planetWA > (OWA + tolerance):
    raise ValueError(f"Planet WA={planetWA:.1f} while IWA = {IWA:.1f} and OWA = {OWA:.1f} lam/D.")

print(f"Planet Working Angle: {planetWA:.2f} λ/D")

# === Contrast Stability Parameters ===
CSprefix = 'MCBE_'
selDeltaC, rawContrast, SystematicCont, initStatRawContrast, \
    rawContrast, IntContStab, ExtContStab = fl.contrastStabilityPars(CSprefix, planetWA, CS_Data)

print(f"Raw Contrast: {rawContrast:.3e}")
print(f"Selected Delta Contrast: {selDeltaC:.3e}")

# === Coronagraph Slice Parameters ===
cg = fl.coronagraphParameters(CG_Data.df, planetWA, DPM)

# === Focal Plane Setup ===
opMode = config['instrument']['OpMode'] #scenarioDF.at['OPMODE_IMG_SPEC', 'Latest']
bandWidth = config['instrument']['bandwidth'] #scenarioDF.at['BW', 'Latest']
f_SR, CritLam, detPixSize_m, mpix, pixPlateSc = fl.getFocalPlaneAttributes(
    opMode,
    config,
    DET_CBE_Data,
    lam,
    bandWidth,
    DPM,
    cg.CGdesignWL,
    cg.omegaPSF
)

# === Star Flux ===
inBandFlux0_sum, inBandZeroMagFlux, starFlux = fl.getSpectra(target, lam, bandWidth)
print(f"Star Flux = {starFlux:.3e} ph/s/m^2")

# === Planet Flux ===
fluxRatio = target.albedo * (target.radius_Rjup * uc.jupiterRadius / (target.sma_AU * uc.AU)) ** 2
planetFlux = fluxRatio * starFlux
print(f"Planet Flux = {planetFlux:.3e} ph/s/m^2")

# === Background Zodi and Speckle Flux ===

magLocalZodi = scenarioDF.at['LocZodi_magas2','Latest']
magExoZodi_1AU = scenarioDF.at['ExoZodi_magas2','Latest']
absMag = target.v_mag - 5 * math.log10(target.dist_pc / 10)

locZodiAngFlux = inBandZeroMagFlux * 10 ** (-0.4 * magLocalZodi)
exoZodiAngFlux = inBandZeroMagFlux * 10 ** (-0.4 * (absMag - uc.sunAbsMag + magExoZodi_1AU)) / target.sma_AU**2 * target.exoZodi

exoZodiDistrib = "uniform"  # Options: "lumpy", "uniform", "falloff"

thput, throughput_rates = fl.compute_throughputs(THPT_Data, cg, exoZodiDistrib)

planetThroughput  = throughput_rates["planet"]
speckleThroughput = throughput_rates["speckle"]
locZodiThroughput = throughput_rates["local_zodi"]
exoZodiThroughput = throughput_rates["exo_zodi"]

# Unobscured Collecting Aperture Area (obscuration accounted separately)
Acol = (np.pi / 4.0) * DPM**2

@dataclass
class corePhotonRates:
    planet:  float = planetFlux * planetThroughput * Acol 
    speckle: float = starFlux * rawContrast * cg.PSFpeakI * cg.CGintmpix * speckleThroughput * Acol 
    locZodi: float = locZodiAngFlux * cg.omegaPSF * locZodiThroughput * Acol 
    exoZodi: float = exoZodiAngFlux * cg.omegaPSF * exoZodiThroughput * Acol  
    total:   float = planet + speckle + locZodi + exoZodi
cphrate = corePhotonRates()


# === Frame Time and dQE Calculation ===
desiredRate = 0.1  # e-/pix/frame
tfmin = 1          # min frame time (s)
tfmax = 100        # max frame time (s)

ENF, effReadnoise, frameTime, dQE = fl.compute_frame_time_and_dqe(
    desiredRate, tfmin, tfmax,
    isPhotonCounting, QE_Data, DET_CBE_Data,
    lam, mpix, cphrate.total
)


# Account for the additional noise due to RDI through penalty factors
rdi_penalty = fl.rdi_noise_penalty(target, inBandFlux0_sum, starFlux, scenarioDF)
k_sp  = rdi_penalty["k_sp"]
k_det = rdi_penalty["k_det"]
k_lzo = rdi_penalty["k_lzo"]
k_ezo = rdi_penalty["k_ezo"]

# === Detector noise calculation ===
detNoiseRate = fl.detector_noise_rates(DET_CBE_Data, monthsAtL2, frameTime, mpix, isPhotonCounting)


k_pp = scenarioDF.at['pp_Factor_CBE', 'Latest']
eRatesCore, residSpecRate = fl.compute_variance_rates(
    cphrate=cphrate,
    dQE=dQE,
    ENF=ENF,
    detNoiseRate=detNoiseRate,
    k_sp=k_sp,
    k_det=k_det,
    k_lzo=k_lzo,
    k_ezo=k_ezo,
    f_SR=f_SR,
    starFlux=starFlux,
    selDeltaC=selDeltaC,
    k_pp=k_pp,
    cg=cg,
    speckleThroughput=speckleThroughput,
    Acol=Acol
)

SNRdesired = 5.0
timeToSNR, criticalSNR = fl.compute_tsnr(SNRdesired, eRatesCore, residSpecRate)

print(f"\nTarget SNR = {SNRdesired:.1f} \nCritical SNR = {criticalSNR:.2f}")
print(f"Time to SNR = {timeToSNR/uc.hour:.2f} hours")





