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

import os
import math
from datetime import datetime
import unitsConstants as uc
import cginoiselib as fl
from dataclasses import dataclass
import numpy as np
import yaml
import sys
from prettytable import PrettyTable

current_dir = os.getcwd()
print(f"Working directory: {current_dir}")
current_datetime = datetime.now()
print(f"Run started at: {current_datetime}")

# === Scenario Selection ===
# Load scenario configuration from a YAML file.
scenario_filename = 'SCEN_IMG_NFOV_B1_HLC_TVAC.yaml'
scenFolder = fl.open_folder("EBcsvData", "Scenarios")
try:
    with open(scenFolder[scenario_filename], "r") as file:
        config = yaml.safe_load(file)
except FileNotFoundError:
    print("Error: config.yaml not found!")
    sys.exit(1)
except yaml.YAMLError as e:
    print(f"Error parsing YAML: {e}")
    sys.exit(1)

# === Constants ===
# Define fundamental physical and instrument constants.
DPM = 2.363 * uc.meter
lam = config['instrument']['wavelength']
lamD = lam / DPM
intTimeDutyFactor = config['instrument']['dutyFactor']

print(f"Wavelength: {lam / uc.nm} nm")
print(f"Lambda/D: {lamD / uc.mas:.3f} mas")

# === Define Host Star and Planet ===
# Uses the Target class from cginoiselib to define the exoplanetary system.
target = fl.Target(
    v_mag=5.0,
    dist_pc=10.0,
    specType='g0v',
    phaseAng_deg=65,
    sma_AU=4.1535,
    radius_Rjup=5.6211,
    geomAlb_ag=0.44765,
    exoZodi=1,
)

# === User-defined inputs ===
# Define observation parameters and constraints.
allocatedTime = 100 * uc.hour
monthsAtL2 = 21
frameTime = 10.0 * uc.second
isPhotonCounting = True  # or False, depending on the mode


usableTinteg = intTimeDutyFactor * allocatedTime

sep_mas = target.phaseAng_to_sep(target.sma_AU, target.dist_pc, target.phaseAng_deg)
target.albedo = target.albedo_from_geomAlbedo(target.phaseAng_deg, target.geomAlb_ag)

print(f"Separation: {sep_mas:.0f} mas")
print(f"Albedo: {target.albedo:.3f}")

# === Load Required CSV Files ===
# Load various instrument and atmospheric data from CSV files.
filenameList = fl.getScenFileNames(config)
CG_Data, QE_Data, DET_CBE_Data, STRAY_FRN_Data, THPT_Data, CAL_Data, CS_Data = fl.loadCSVs(filenameList)

# === Planet Working Angle ===
# Calculate the planet's working angle in units of lambda/D.
IWA, OWA = fl.workingAnglePars(CG_Data, CS_Data) # Inner and Outer Working Angles
planetWA = sep_mas * uc.mas / lamD  # Planet working angle

# Create a table with column headers
table = PrettyTable()
table.field_names = ['planet WA', 'phase, deg', 'dist, pc', 'sma, AU', 'sep, mas', 'lam/D, mas', "IWA", "OWA"]
table.add_row( [f'{planetWA:.2f}',f'{target.phaseAng_deg:.2f}', f"{target.dist_pc:.2f}", f'{target.sma_AU:.2f}', f'{sep_mas:.2f}', f'{lamD/uc.mas:.2f}', f'{IWA:.2f}',f'{OWA:.2f}' ] )
print(table)



tolerance = 0.05  # Tolerance for working angle adjustment

# Adjust planetWA to be exactly IWA or OWA if it's very close, or raise error if outside range.
if (IWA - tolerance) <= planetWA <= IWA:
    planetWA = IWA
elif OWA <= planetWA <= (OWA + tolerance):
    planetWA = OWA
elif planetWA < (IWA - tolerance) or planetWA > (OWA + tolerance):
    raise ValueError(f"Planet WA={planetWA:.1f} while IWA = {IWA:.1f} and OWA = {OWA:.1f} lam/D.")

print(f"Planet Working Angle: {planetWA:.2f} λ/D")

# === Contrast Stability Parameters ===
# Determine contrast stability parameters from loaded data.
CSprefix = 'ICBE_'
selDeltaC, rawContrast, SystematicCont, initStatRawContrast, \
    rawContrast, IntContStab, ExtContStab = fl.contrastStabilityPars(CSprefix, planetWA, CS_Data)

print(f"Raw Contrast: {rawContrast:.3e}")
print(f"Selected Delta Contrast: {selDeltaC:.3e}")

# === Coronagraph Slice Parameters ===
# Extract coronagraph parameters for the calculated working angle.
cg = fl.coronagraphParameters(CG_Data.df, planetWA, DPM)

# === Focal Plane Setup ===
# Configure focal plane attributes based on operational mode and instrument config.
opMode = config['instrument']['OpMode']
bandWidth = config['instrument']['bandwidth']
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

# === Sensor Inner Workings ===
# Caveats concerning the capabilities of the sensor being used
det_FWCserial = config['detector']['FWC_serial']

# === Star Flux ===
# Calculate the flux from the host star.
inBandFlux0_sum, inBandZeroMagFlux, starFlux = fl.getSpectra(target, lam, bandWidth)
print(f"Star Flux = {starFlux:.3e} ph/s/m^2")
TimeonRefStar_tRef_per_tTar = 0.25

# === Planet Flux ===
# Calculate the flux from the exoplanet.
fluxRatio = target.albedo * (target.radius_Rjup * uc.jupiterRadius / (target.sma_AU * uc.AU)) ** 2
planetFlux = fluxRatio * starFlux
print(f"Planet Flux = {planetFlux:.3e} ph/s/m^2")

# === Background Zodi and Speckle Flux ===
# Calculate flux from local and exo-zodiacal light, and speckles.

magLocalZodi = config['instrument']['LocZodi_magas2']    # Magnitude of local zodiacal light
magExoZodi_1AU = config['instrument']['ExoZodi_magas2']  # Magnitude of exo-zodiacal light at 1 AU
absMag = target.v_mag - 5 * math.log10(target.dist_pc / 10) # Absolute magnitude of the host star

locZodiAngFlux = inBandZeroMagFlux * 10 ** (-0.4 * magLocalZodi)
exoZodiAngFlux = inBandZeroMagFlux * 10 ** (-0.4 * (absMag - uc.sunAbsMag + magExoZodi_1AU)) / target.sma_AU**2 * target.exoZodi

exoZodiDistrib = "uniform"  # Distribution model for exo-zodiacal light. Options: "lumpy", "uniform", "falloff"

# Compute various throughputs based on the coronagraph and exo-zodi distribution.
thput, throughput_rates = fl.compute_throughputs(THPT_Data, cg, exoZodiDistrib)

planetThroughput  = throughput_rates["planet"]
speckleThroughput = throughput_rates["speckle"]
locZodiThroughput = throughput_rates["local_zodi"]
exoZodiThroughput = throughput_rates["exo_zodi"]

# Unobscured Collecting Aperture Area (obscuration accounted for separately in throughputs)
Acol = (np.pi / 4.0) * DPM**2

@dataclass
class corePhotonRates:
    """
    Dataclass to store the core photon rates from different sources at the detector.

    Attributes:
        planet (float): Photon rate from the exoplanet (photons/second).
        speckle (float): Photon rate from residual star light (speckles) (photons/second).
        locZodi (float): Photon rate from local zodiacal light (photons/second).
        exoZodi (float): Photon rate from exo-zodiacal light (photons/second).
        total (float): Sum of all core photon rates (photons/second).
    """
    planet:  float = planetFlux * planetThroughput * Acol 
    speckle: float = starFlux * rawContrast * cg.PSFpeakI * cg.CGintmpix * speckleThroughput * Acol 
    locZodi: float = locZodiAngFlux * cg.omegaPSF * locZodiThroughput * Acol 
    exoZodi: float = exoZodiAngFlux * cg.omegaPSF * exoZodiThroughput * Acol  
    total:   float = planet + speckle + locZodi + exoZodi
cphrate = corePhotonRates()


# === Frame Time and dQE Calculation ===
# Determine the optimal frame time and detector's differential Quantum Efficiency (dQE).
desiredRate = 0.1  # e-/pix/frame
tfmin = 3          # min frame time (s)
tfmax = 100        # max frame time (s)

ENF, effReadnoise, frameTime, dQE = fl.compute_frame_time_and_dqe(
    desiredRate, tfmin, tfmax,
    isPhotonCounting, QE_Data, DET_CBE_Data,
    lam, mpix, cphrate.total, 
    det_FWCserial
)
print(f"Calculated Frame Time: {frameTime:.2f} s")
print(f"Differential Quantum Efficiency (dQE): {dQE:.3f}")
print(f"Excess Noise Factor (ENF): {ENF:.2f}")


# Account for the additional noise due to RDI through penalty factors
# These factors increase the variance from different noise sources.
rdi_penalty = fl.rdi_noise_penalty(target, inBandFlux0_sum, starFlux, TimeonRefStar_tRef_per_tTar)
k_sp  = rdi_penalty["k_sp"]
k_det = rdi_penalty["k_det"]
k_lzo = rdi_penalty["k_lzo"]
k_ezo = rdi_penalty["k_ezo"]

# === Detector noise calculation ===
# Calculate various detector noise rates.
detNoiseRate = fl.detector_noise_rates(DET_CBE_Data, monthsAtL2, frameTime, mpix, isPhotonCounting)

# === Variance and SNR Calculation ===
# Compute electron rates for variance calculation and residual speckle rate.
k_pp = config['instrument']['pp_Factor_CBE']
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

if timeToSNR > usableTinteg:
    print(f"Warning: Time to SNR ({timeToSNR/uc.hour:.2f} hrs) exceeds usable integration time ({usableTinteg/uc.hour:.2f} hrs).")
elif timeToSNR <= 0:
    print(f"Warning: Calculated Time to SNR is not positive ({timeToSNR/uc.hour:.2f} hrs). Check input parameters and intermediate calculations.")
else:
    print("Observation is feasible within the allocated time.")

print(f"Run completed at: {datetime.now()}")

