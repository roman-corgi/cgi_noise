"""
Library of core functions for the EB performance modeling pipeline.

This module provides structured access to scenario loading, throughput computation,
optical and detector models, noise variance calculations, and astrophysical fluxes.
All units are assumed to follow SI unless otherwise noted, and helper constants
are provided in the 'unitsConstants' module. 
"""
from dataclasses import dataclass
from pathlib import Path
import os
import unitsConstants as uc
import math
from loadCSVrow import loadCSVrow
from dataclasses import dataclass, asdict

def open_folder(*folders):
    """Opens a directory and returns a dictionary of file paths keyed by filenames."""
    filenamedir = Path(os.getcwd())
    folder = Path(filenamedir, *folders)
    return {file.name: file for file in folder.iterdir() if file.is_file()}


def getScenFileNames_DRM(scenarioData):
    filenamedir = Path(os.getcwd())
    filenameList = []
    ffList = [
        ("Photometry", "CoronagraphFile"),
        ('Photometry', 'QE_Curve_file'),
        ('Photometry', 'DetModelFile_CBE'),
        ('Photometry', 'StrayLightFRNfile'),
        ('Photometry', 'ThroughputFile'),
        ('Calibration', 'CalibrationFile'),
        ('Cstability', 'ContrastStabilityFile')
    ]
    for folder, key in ffList:
        name = scenarioData.loc[key, 'Latest'] + ".csv"
        path = filenamedir / "EBcsvData" / folder / name
        filenameList.append(str(path))
    return filenameList


def loadCSVs(filenameList):
    return [loadCSVrow(f) for f in filenameList]


def workingAnglePars(CG_Data, CS_Data):
    IWAc = CG_Data.df.at[0, 'rlamD']
    IWAs = CS_Data.df.at[0, 'r_lam_D']
    OWAc = CG_Data.df['rlamD'].iloc[-1]
    OWAs = CS_Data.df['r_lam_D'].iloc[-1]
    return max(IWAs, IWAc), min(OWAs, OWAc)


def contrastStabilityPars(CSprefix, planetWA, CS_Data):
    tol = 0.05
    indCS = CS_Data.df['r_lam_D'].searchsorted(planetWA + tol) - 1

    headers = CS_Data.df.columns.tolist()
    nCols = len(headers)
    fnARC = CSprefix + "AvgRawContrast"
    fnECS = CSprefix + "ExtContStab"
    fnICS = CSprefix + "IntContStab"
    fnSC  = CSprefix + "SystematicC"
    fnISRC = CSprefix + "InitStatContrast"

    ExtContStab = CS_Data.df.at[indCS, fnECS] * uc.ppb
    IntContStab = CS_Data.df.at[indCS, fnICS] * uc.ppb
    rawContrast = CS_Data.df.at[indCS, fnARC] * uc.ppb
    initStatRawContrast = CS_Data.df.at[indCS, fnISRC] * uc.ppb

    if nCols == 16 and 'SystematicC' in headers[13]:
        SystematicCont = CS_Data.df.at[indCS, fnSC] * uc.ppb
        selDeltaC = math.sqrt((ExtContStab**2) + (IntContStab**2) + (SystematicCont**2))
    elif nCols == 13:
        SystematicCont = 0
        selDeltaC = math.sqrt((ExtContStab**2) + (IntContStab**2))
    else:
        raise IndexError('The contrast stability file referenced is not formatted as expected.')

    return selDeltaC, rawContrast, SystematicCont, initStatRawContrast, rawContrast, IntContStab, ExtContStab


def getFocalPlaneAttributes(opMode, scenarioData, DET_CBE_Data, lam, bandWidth, DPM, CGdesignWL, omegaPSF):
    FocalPlaneAtt = loadCSVrow(Path(os.getcwd(), 'EBcsvData', 'CONST_SNR_FPattributes.csv'))
    AmiciPar = loadCSVrow(Path(os.getcwd(), 'EBcsvData', 'CONST_Amici_parameters.csv'))

    detPixSize_m = DET_CBE_Data.df.at[0, 'PixelSize_m']

    if opMode == "SPEC":
        try:
            resolution = scenarioData.at['R_required', 'Latest']
            f_SR = 1 / (resolution * bandWidth)
        except:
            resolution = 0.0001
            f_SR = -1

        CritLam = FocalPlaneAtt.df.at[1, 'Critical_Lambda_m']
        compbeamD_m = AmiciPar.df.loc[0, 'compressd_beam_diamtr_m']
        fnlFocLen_m = AmiciPar.df.loc[0, 'final_focal_length_m']
        Fno = fnlFocLen_m / compbeamD_m
        pixPerlamD = lam * Fno / detPixSize_m
        PSF_x_lamD = AmiciPar.df.loc[0, 'PSF_core_x_extent_lamD']
        PSF_y_lamD = AmiciPar.df.loc[0, 'PSF_core_y_extent_lamD']
        xpixPerCor = PSF_x_lamD * 2 * pixPerlamD
        ypixPerCor = PSF_y_lamD * 2 * pixPerlamD
        Rlamsq = AmiciPar.df.loc[0, 'lam_squared']
        Rlam = AmiciPar.df.loc[0, 'lam']
        Rconst = AmiciPar.df.loc[0, 'constant']
        ResPowatPSF = Rconst + Rlam * lam + Rlamsq * lam**2
        dpix_dlam = ResPowatPSF * xpixPerCor / lam
        xpixPerSpec = dpix_dlam * lam / resolution
        mpix = xpixPerSpec * ypixPerCor
        pixPlateSc = CritLam / DPM / 2 / uc.mas

    elif opMode == "IMG":
        f_SR = 1
        CritLam = FocalPlaneAtt.df.at[0, 'Critical_Lambda_m']
        mpix = omegaPSF * uc.arcsec**2 * (lam / CGdesignWL)**2 * (2 * DPM / CritLam)**2
        pixPlateSc = CritLam / DPM / 2 / uc.mas
    else:
        raise Exception("getFocalPlaneAttributes: Valid Operational Modes are IMG and SPEC")

    return f_SR, CritLam, detPixSize_m, mpix, pixPlateSc


@dataclass
class CGParameters:
    CGcoreThruput: float
    PSFpeakI: float
    omegaPSF: float
    CGintSamp: float
    CGradius_arcsec: float
    CGdesignWL: float
    CGintmpix: float
    CG_PSFarea_sqlamD: float
    CGintensity: float
    CG_occulter_transmission: float
    CGcontrast: float
    CGtauPol: float = 1.0


@dataclass
class Target:
    v_mag: float
    dist_pc: float
    specType: str
    phaseAng_deg: float
    sma_AU: float
    radius_Rjup: float
    geomAlb_ag: float
    exoZodi: float
    albedo: float = None

    @staticmethod
    def phaseAng_to_sep(sma_AU, dist_pc, phaseAng_deg):
        sep_mas = ((sma_AU * uc.AU * math.sin(math.radians(phaseAng_deg))) / (dist_pc * uc.pc)) / uc.mas
        return sep_mas

    @staticmethod
    def fluxRatio_to_deltaMag(fluxRatio):
        return (-2.5) * math.log10(fluxRatio)

    @staticmethod
    def deltaMag_to_fluxRatio(deltaMag):
        return 10 ** (-0.4 * deltaMag)

    @staticmethod
    def fluxRatio_SMA_rad_to_albedo(fluxRatio, sma_AU, radius_Rjup):
        return fluxRatio * (sma_AU * uc.AU / (radius_Rjup * uc.jupiterRadius)) ** 2


def coronagraphParameters(cg_df, planetWA, DPM):
    CGtauPol = 1
    indWA = cg_df[(cg_df.rlamD <= planetWA)]['rlamD'].idxmax()

    CGcoreThruput = cg_df.loc[indWA, 'coreThruput'] * CGtauPol
    PSFpeakI = cg_df.loc[indWA, 'PSFpeak'] * CGtauPol
    omegaPSF = cg_df.loc[indWA, 'area_sq_arcsec']
    CGintSamp = cg_df.loc[2, 'rlamD'] - cg_df.loc[1, 'rlamD']
    CGradius_arcsec = cg_df.at[indWA, 'r_as']

    CGdesignWL = DPM * cg_df.iloc[0, 1] * uc.arcsec / cg_df.iloc[0, 0]
    CGintmpix = omegaPSF * (uc.arcsec**2) / ((CGintSamp * CGdesignWL / DPM)**2)
    CG_PSFarea_sqlamD = omegaPSF / (CGdesignWL / uc.arcsec)**2

    CGintensity = cg_df.loc[indWA, 'I']
    CG_occulter_transmission = cg_df.at[indWA, 'occTrans'] * CGtauPol
    CGcontrast = cg_df.loc[indWA, 'contrast']

    return CGParameters(
        CGcoreThruput=CGcoreThruput,
        PSFpeakI=PSFpeakI,
        omegaPSF=omegaPSF,
        CGintSamp=CGintSamp,
        CGradius_arcsec=CGradius_arcsec,
        CGdesignWL=CGdesignWL,
        CGintmpix=CGintmpix,
        CG_PSFarea_sqlamD=CG_PSFarea_sqlamD,
        CGintensity=CGintensity,
        CG_occulter_transmission=CG_occulter_transmission,
        CGcontrast=CGcontrast,
        CGtauPol=CGtauPol
    )

def getSpectra(target, lam, bandWidth):
    spectra_path = Path(os.getcwd(), 'EBcsvData', 'Spectra', 'SPECTRA_ALL_BPGS.csv')
    SPECTRA_Data = loadCSVrow(spectra_path)

    bandRange = SPECTRA_Data.df[abs(SPECTRA_Data.df['Wavelength_m'] - lam) <= (0.5 * bandWidth * lam)]
    onlySpec = bandRange.drop(['Wavelength_m', 'E_ph_J'], axis=1)

    Ephot = uc.h_planck * uc.c_light / lam
    onlySpecEphot = onlySpec.apply(lambda x: x / Ephot, axis=1, result_type='broadcast')

    deltaLambda = SPECTRA_Data.df.at[2, 'Wavelength_m'] - SPECTRA_Data.df.at[1, 'Wavelength_m']
    inBandFlux0_sum = onlySpecEphot.sum(axis=0) * deltaLambda

    inBandZeroMagFlux = inBandFlux0_sum.at[target.specType]
    starFlux = inBandZeroMagFlux * 10 ** (-0.4 * target.v_mag)

    return inBandFlux0_sum, inBandZeroMagFlux, starFlux

@dataclass
class DetNoiseRates:
    dark: float
    CIC: float
    read: float

def detector_noise_rates(DET_CBE_Data, monthsAtL2, frameTime, mpix, isPhotonCounting):
    missionFraction = monthsAtL2 / DET_CBE_Data.df.at[0, 'DetEOL_mos']
    detDarkBOM = DET_CBE_Data.df.at[0, 'DarkBOM_e_per_pix_per_hr']
    detDarkEOM = DET_CBE_Data.df.at[0, 'DarkEOM_e_per_pix_per_hr']
    dark_per_hr = detDarkBOM + missionFraction * (detDarkEOM - detDarkBOM)
    dark_per_s = dark_per_hr / 3600

    detCIC1 = DET_CBE_Data.df.at[0, 'CICatGain1BOM_e_per_pix_per_fr']
    detCIC2 = DET_CBE_Data.df.at[0, 'CICatGain2BOM_e_per_pix_per_fr']
    gain1 = DET_CBE_Data.df.at[0, 'Gain1BOM']
    gain2 = DET_CBE_Data.df.at[0, 'Gain2BOM']
    CIC_degradation = DET_CBE_Data.df.at[0, 'CICdegradationEOM']
    EMgain = DET_CBE_Data.df.at[0, 'EMGain']

    CIC_rate = ((detCIC2 - detCIC1) / (gain2 - gain1)) * EMgain + (
        detCIC1 - ((detCIC2 - detCIC1) / (gain2 - gain1)) * gain1
    ) * (1 + missionFraction * (CIC_degradation - 1))
    CIC_per_s = CIC_rate / frameTime

    if isPhotonCounting:
        readNoise = 0
    else:
        detCamRead = DET_CBE_Data.df.at[0, 'ReadNoise_e']
        EMgain = DET_CBE_Data.df.at[0, 'EMGain']
        readNoise = detCamRead / EMgain

    read_noise_per_s = (mpix / frameTime) * (readNoise ** 2)

    return DetNoiseRates(
        dark = dark_per_s,
        CIC  = CIC_per_s,
        read = read_noise_per_s
    )



@dataclass
class Throughput:
    refl: float
    filt: float
    polr: float
    core: float
    occt: float

def compute_throughputs(THPT_Data, cg, ezdistrib="falloff"):
    """
    Compute optical throughputs and exozodi factors.

    Parameters:
    - THPT_Data: loaded CSV row with throughput data.
    - cg: CGParameters object.
    - ezdistrib: Exo-Zodi distribution: one of {"lumpy", "uniform", "falloff"}

    Returns:
    - Throughput instance.
    - Dictionary with total throughputs: planet, speckle, local_zodi, exo_zodi
    """

    # Select the appropriate distribution factor for exozodi
    dist_map = {
        "lumpy": 0.49,
        "uniform": 1.00,
        "falloff": 0.74
    }

    if ezdistrib not in dist_map:
        raise ValueError(f"Invalid ezodistribution: {ezdistrib}. Must be 'lumpy', 'uniform', or 'falloff'.")

    distFactor = dist_map[ezdistrib]

    thput = Throughput(
        refl=THPT_Data.df.at[0, 'Pupil_Transmission']
             * THPT_Data.df.at[0, 'CBE_OTAplusTCA']
             * THPT_Data.df.at[0, 'CBE_CGI'],
        filt=1.0,
        polr=1.0,
        core=THPT_Data.df.at[0, 'CBE_Core'],
        occt=cg.CG_occulter_transmission
    )

    planetThroughput  = thput.refl * thput.filt * thput.polr * thput.core
    speckleThroughput = thput.refl * thput.filt * thput.polr * thput.polr
    locZodiThroughput = thput.refl * thput.filt * thput.polr * thput.occt
    exoZodiThroughput = locZodiThroughput * distFactor

    return thput, {
        "planet": planetThroughput,
        "speckle": speckleThroughput,
        "local_zodi": locZodiThroughput,
        "exo_zodi": exoZodiThroughput
    }


def rdi_noise_penalty(target, inBandFlux0_sum, starFlux, scenarioDF,
                               RefStarSpecType='a0v', RefStarDist=10, RefStarVmag=3.0):
    """
    Compute noise penalty factors for Reference Differential Imaging (RDI).
    
    Parameters:
    - target: Target dataclass instance.
    - inBandFlux0_sum: zero-magnitude flux per spectral type (Series).
    - starFlux: target star flux.
    - scenarioDF: scenario DataFrame with timing info.
    - RefStarSpecType: spectral type of the reference star (default 'a0v').
    - RefStarDist: distance to the reference star in parsecs (default 10).
    - RefStarVmag: V magnitude of the reference star (default 3.0).
    
    Returns:
    - Dictionary of penalty factors: k_sp, k_det, k_lzo, k_ezo
    """

    RefStarinBandZeroMagFlux = inBandFlux0_sum.at[RefStarSpecType]
    RefStarAbsMag = RefStarVmag - 5 * math.log10(RefStarDist / 10)
    RefStarFlux = RefStarinBandZeroMagFlux * (10 ** (-0.4 * RefStarVmag))
    BrightnessRatio = RefStarFlux / starFlux

    timeRatio = scenarioDF.at['TimeonRefStar_tRef_per_tTar', 'Latest']
    betaRDI = 1 / (BrightnessRatio * timeRatio)

    k_sp = 1 + betaRDI
    k_det = 1 + betaRDI**2 * timeRatio
    k_lzo = k_det
    k_ezo = k_sp

    return {
        "k_sp": k_sp,
        "k_det": k_det,
        "k_lzo": k_lzo,
        "k_ezo": k_ezo
    }


def compute_frame_time_and_dqe(
    desiredRate, tfmin, tfmax,
    isPhotonCounting, QE_Data, DET_CBE_Data,
    lam, mpix, cphrate_total
):
    """
    Compute frame time and effective quantum efficiency (dQE) based on photon counting mode.

    Parameters:
    - desiredRate: target e-/pix/frame
    - tfmin, tfmax: min/max allowed frame time (seconds)
    - isPhotonCounting: True if PC mode, else false
    - QE_Data: QE curve CSV data
    - DET_CBE_Data: detector model CSV data
    - lam: central wavelength (meters)
    - mpix: number of pixels integrated
    - cphrate_total: total core photon rate (e-/s)

    Returns:
    - ENF: excess noise factor
    - effReadnoise: effective read noise (e-/s)
    - frameTime: calculated frame time (seconds)
    - dQE: effective quantum efficiency
    """

    det_QE = QE_Data.df.loc[QE_Data.df['lambda_nm'] <= (lam / uc.nm), 'QE_at_neg100degC'].iloc[-1]
    det_EMgain = DET_CBE_Data.df.at[0, 'EMGain']
    det_readnoise = DET_CBE_Data.df.at[0, 'ReadNoise_e']
    det_PCthresh = DET_CBE_Data.df.at[0, 'PCThresh_nsigma']
    det_FWCserial = 90000

    if isPhotonCounting:
        ENF = 1.0
        effReadnoise = 0.0
        frameTime = round(min(tfmax, max(tfmin, desiredRate / (cphrate_total * det_QE / mpix))), 1)
        approxPerPixelPerFrame = frameTime * cphrate_total * det_QE / mpix
        eff_coincidence = (1 - math.exp(-approxPerPixelPerFrame)) / approxPerPixelPerFrame if approxPerPixelPerFrame > 0 else 1.0
        eff_thresholding = math.exp(-det_PCthresh * det_readnoise / det_EMgain)
        dQE = det_QE * eff_coincidence * eff_thresholding
    else:
        ENF = math.sqrt(2)
        effReadnoise = det_readnoise / det_EMgain
        Nsigma = 3
        NEE = Nsigma * ENF * det_EMgain
        y_crit = ((NEE**2 + 2 * det_FWCserial) - math.sqrt(NEE**4 + 4 * NEE**2 * det_FWCserial)) / 2
        tfr_crit = y_crit / (cphrate_total * det_QE / mpix)
        frameTime = min(tfmax, max(tfmin, math.floor(tfr_crit)))
        dQE = det_QE

    return ENF, effReadnoise, frameTime, dQE


@dataclass
class VarianceRates:
    planet: float
    speckle: float
    locZodi: float
    exoZodi: float
    detDark: float
    detCIC: float
    detRead: float

    @property
    def total(self):
        return sum(asdict(self).values())

    def __repr__(self):
        fields = asdict(self)
        fields_str = ", ".join([f"{k}={v:.3e}" for k, v in fields.items()])
        return f"VarianceRates({fields_str}, total={self.total:.3e})"

def compute_variance_rates(cphrate, dQE, ENF, detNoiseRate, k_sp, k_det, k_lzo, k_ezo,
                           f_SR, starFlux, selDeltaC, k_pp, cg, speckleThroughput, Acol):
    """
    Compute variance rates and residual speckle rate after post-processing.

    Parameters:
    - All as before, plus:
      - f_SR: spectral resolution factor
      - starFlux: flux of the target star
      - selDeltaC: selected delta contrast (unitless)
      - k_pp: post-processing factor (e.g., 30 for 30x speckle suppression)
      - cg: CGParameters object
      - speckleThroughput: total system throughput for speckle
      - Acol: collecting area (m^2)

    Returns:
    - VarianceRates object
    - residSpecRate: residual speckle rate (e-/s)
    """

    residSpecRate = (
        f_SR * starFlux * (selDeltaC / k_pp) *
        cg.PSFpeakI * cg.CGintmpix *
        speckleThroughput * Acol * dQE
    )

    rates = VarianceRates(
        planet  = ENF**2 * cphrate.planet * dQE,
        speckle = ENF**2 * cphrate.speckle * dQE * k_sp,
        locZodi = ENF**2 * cphrate.locZodi * dQE * k_lzo,
        exoZodi = ENF**2 * cphrate.exoZodi * dQE * k_ezo,
        detDark = ENF**2 * detNoiseRate.dark * k_det,
        detCIC  = ENF**2 * detNoiseRate.CIC  * k_det,
        detRead = detNoiseRate.read * k_det,
    )

    return rates, residSpecRate


def compute_tsnr(SNRdesired, eRatesCore, residSpecRate):
    """
    Compute the required integration time and critical SNR.

    Parameters:
    - SNRdesired: Target signal-to-noise ratio (float)
    - eRatesCore: VarianceRates object (must include .planet and .total)
    - residSpecRate: residual speckle rate in electrons/sec

    Returns:
    - timeToSNR: integration time in seconds to reach SNRdesired
    - criticalSNR: maximum achievable SNR given residual speckle
    """

    denom = eRatesCore.planet**2 - SNRdesired**2 * residSpecRate**2
    if denom <= 0:
        raise ValueError("SNR condition is not achievable with given residual speckle level.")

    timeToSNR = SNRdesired**2 * eRatesCore.total / denom
    criticalSNR = eRatesCore.planet / residSpecRate

    return timeToSNR, criticalSNR
