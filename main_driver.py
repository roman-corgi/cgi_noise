# main_driver.py - Converted to a script version (no main function)

import os
import math
from datetime import datetime
from loadXLcol import loadXLcol
import unitsConstants as uc
import library as fl
from dataclasses import dataclass, asdict
import numpy as np

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

# === Contrast Stability Parameters ===
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

locZodiAngFlux = inBandZeroMagFlux * 10 ** (-0.4 * magLocalZodi)
exoZodiAngFlux = inBandZeroMagFlux * 10 ** (-0.4 * (absMag - uc.sunAbsMag + magExoZodi_1AU)) / target.sma_AU**2 * target.exoZodi

exoZodiDistrib = "falloff"  # Options: "lumpy", "uniform", "falloff"

thput, throughput_rates = fl.compute_throughputs(THPT_Data, cg, exoZodiDistrib)

planetThroughput  = throughput_rates["planet"]
speckleThroughput = throughput_rates["speckle"]
locZodiThroughput = throughput_rates["local_zodi"]
exoZodiThroughput = throughput_rates["exo_zodi"]

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
    frameTime = round(min(tfmax, max(tfmin, desiredRate / (cphrate.total * det_QE / mpix))), 1)
    approxPerPixelPerFrame = frameTime * cphrate.total * det_QE / mpix
    eff_coincidence = (1 - math.exp(-approxPerPixelPerFrame)) / approxPerPixelPerFrame if approxPerPixelPerFrame > 0 else 1.0
    eff_thresholding = math.exp(-det_PCthresh * det_readnoise / det_EMgain)
    dQE = det_QE * eff_coincidence * eff_thresholding
else:
    ENF = math.sqrt(2)
    effReadnoise = det_readnoise / det_EMgain
    Nsigma = 3
    NEE = Nsigma * ENF * det_EMgain
    y_crit = ((NEE**2 + 2 * det_FWCserial) - math.sqrt(NEE**4 + 4 * NEE**2 * det_FWCserial)) / 2
    tfr_crit = y_crit / (cphrate.total * det_QE / mpix)
    frameTime = min(tfmax, max(tfmin, math.floor(tfr_crit)))
    dQE = det_QE


# === Reference Star Noise Penalty in RDI ===
RefStarSpecType ='a0v'
RefStarDist = 10 # pc
RefStarVmag = 3.0

RefStarinBandZeroMagFlux = inBandFlux0_sum.at[RefStarSpecType]
RefStarAbsMag = RefStarVmag - 5*math.log10(RefStarDist/10)
RefStarDeltaMag = target.v_mag - RefStarVmag
RefStarFlux = RefStarinBandZeroMagFlux*(10**((-0.4)*RefStarVmag))
BrightnessRatio = RefStarFlux/starFlux

timeRatio = scenarioDF.at['TimeonRefStar_tRef_per_tTar', 'Latest']

# RDI normalization factor 
betaRDI = 1 / (BrightnessRatio*timeRatio)

k_sp  = 1 + betaRDI
k_det = 1 + betaRDI**2 * timeRatio
k_lzo = k_det
k_ezo = k_sp

# === Detector noise calculation ===
detNoiseRate = fl.detector_noise_rates(DET_CBE_Data, monthsAtL2, frameTime, mpix, isPhotonCounting)

# noise variance rates class
@dataclass
class varianceRates:
    planet:  float  = ENF**2 * cphrate.planet * dQE 
    speckle: float  = ENF**2 * cphrate.speckle * dQE * k_sp
    locZodi: float  = ENF**2 * cphrate.locZodi * dQE * k_lzo
    exoZodi: float  = ENF**2 * cphrate.exoZodi * dQE * k_ezo
    detDark: float  = ENF**2 * detNoiseRate.dark * k_det
    detCIC:  float  = ENF**2 * detNoiseRate.CIC  * k_det
    detRead: float =           detNoiseRate.read * k_det

    @property
    def total(self):
        return sum(asdict(self).values())

    def __repr__(self):
        fields = asdict(self)
        fields_str = ", ".join([f"{k}={v}" for k, v in fields.items()])
        return f"varianceRates({fields_str}, total={self.total})"

# object instance: electron rates in core, including a "total" method that is a computed attribute
eRatesCore = varianceRates();



import sys
sys.exit()






# # contrast instability causes a residual post-differential imaging speckle, with associated contrast
# residSpecRate = starFlux * cstab_REQ.rss * cg.PSF_pk * mpixModel * speckleThroughput * Acol * dQE


# timeToSNR = SNRdesired**2 * eRatesCore.total / (eRatesCore.planet**2 - SNRdesired**2 * residSpecRate**2)


# criticalSNR = eRatesCore.planet / residSpecRate

# print(f"\nTarget SNR = {SNRdesired:.1f} \nCritical SNR = {criticalSNR:.2f}")
# print(f"Time to SNR = {timeToSNR/uc.hour:.2f} hours")



 















