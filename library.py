from dataclasses import dataclass
from pathlib import Path
import os
import unitsConstants as uc
import math
import pandas as pd
from loadCSVrow import loadCSVrow


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


def make_cg_parameters(cg_df, planetWA, DPM):
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

    # ENF = 1.0 if isPhotonCounting else math.sqrt(2)
    # total_noise_rate = ENF ** 2 * mpix * (dark_per_s + CIC_per_s) + read_noise_per_s

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


