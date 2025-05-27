# main_driver.py - Converted to a script version (no main function)

import os
import math
from datetime import datetime
from loadXLcol import loadXLcol
import unitsConstants as uc
import library as fl

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
filenameList = fl.getScenFileNames_DRM(scenarioDF)
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

# === Throughput Parameters ===
CSprefix = 'MCBE_'
selDeltaC, rawContrast, SystematicCont, initStatRawContrast, \
    rawContrast, IntContStab, ExtContStab = fl.contrastStabilityPars(CSprefix, planetWA, CS_Data)

print(f"Raw Contrast: {rawContrast:.3e}")
print(f"Selected Delta Contrast: {selDeltaC:.3e}")

# === Coronagraph Slice Parameters ===
cg = fl.make_cg_parameters(CG_Data.df, planetWA, DPM)

# === Focal Plane Setup ===
opMode = scenarioDF.at['OPMODE_IMG_SPEC', 'Latest']
bandWidth = scenarioDF.at['BW', 'Latest']
f_SR, CritLam, detPixSize_m, mpix, pixPlateSc = fl.getFocalPlaneAttributes(
    opMode,
    scenarioDF,
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

loZodiAngFlux = inBandZeroMagFlux * 10 ** (-0.4 * magLocalZodi)
exoZodiAngFlux = target.exoZodi * inBandZeroMagFlux * 10 ** (-0.4 * (absMag - uc.sunAbsMag + magExoZodi_1AU)) / target.sma_AU**2

ZodiFlux = (loZodiAngFlux + exoZodiAngFlux) * cg.omegaPSF
loZodiFlux = loZodiAngFlux * cg.omegaPSF
exoZodiFlux = exoZodiAngFlux * cg.omegaPSF

rate_speckleBkg = f_SR * starFlux * rawContrast * cg.PSFpeakI * cg.CGintmpix * cg.CGtauPol * cg.omegaPSF

print(f"Zodi Flux = {ZodiFlux:.3e} ph/s/m^2")
print(f"  Local Zodi Flux  = {loZodiFlux:.3e}")
print(f"  Exo-Zodi Flux    = {exoZodiFlux:.3e}")
print(f"Speckle Rate = {rate_speckleBkg:.3e} e-/s")

# === Placeholder Output ===
print("Speckle rate calculation complete. Ready for detector noise and SNR estimation.")




# === Frame Time and dQE Calculation ===
desiredRate = 0.1  # e/pix/frame — user-defined efficiency target
tfmin = 1          # minimum allowable frame time (s)
tfmax = 100        # maximum allowable frame time (s)

det_QE = QE_Data.df.loc[QE_Data.df['lambda_nm'] <= (lam / uc.nm), 'QE_at_neg100degC'].iloc[-1]
det_EMgain = DET_CBE_Data.df.at[0, 'EMGain']
det_readnoise = DET_CBE_Data.df.at[0, 'ReadNoise_e']
det_PCthresh = DET_CBE_Data.df.at[0, 'PCThresh_nsigma']
det_FWCserial = 90000

if isPhotonCounting:
    ENF = 1.0
    effReadnoise = 0.0
    rate_total = rate_speckleBkg  # or use sum of all incident background rates
    frameTime = round(min(tfmax, max(tfmin, desiredRate / (rate_total / mpix * det_QE))), 1)
    approxPerPixelPerFrame = frameTime * rate_total * det_QE / mpix
    eff_coincidence = (1 - math.exp(-approxPerPixelPerFrame)) / approxPerPixelPerFrame if approxPerPixelPerFrame > 0 else 1.0
    eff_thresholding = math.exp(-det_PCthresh * det_readnoise / det_EMgain)
    dQE = det_QE * eff_coincidence * eff_thresholding
else:
    ENF = math.sqrt(2)
    effReadnoise = det_readnoise / det_EMgain
    Nsigma = 3
    NEE = Nsigma * ENF * det_EMgain
    y_crit = ((NEE**2 + 2 * det_FWCserial) - math.sqrt(NEE**4 + 4 * NEE**2 * det_FWCserial)) / 2
    rate_total = rate_speckleBkg
    tfr_crit = y_crit / (rate_total / mpix * det_QE)
    frameTime = min(tfmax, max(tfmin, math.floor(tfr_crit)))
    dQE = det_QE





# === Detector noise calculation ===
det_noise = fl.compute_detector_noise(DET_CBE_Data, monthsAtL2, frameTime, mpix, isPhotonCounting)

print(f"Total Detector Noise Rate = {det_noise.total_noise_rate:.3e} e-/s")
print(f"  Dark Current Rate  = {det_noise.dark_current_per_s:.3e} e-/s")
print(f"  CIC Noise Rate     = {det_noise.CIC_noise_per_s:.3e} e-/s")
print(f"  Read Noise Rate    = {det_noise.read_noise_per_s:.3e} e-/s")





















