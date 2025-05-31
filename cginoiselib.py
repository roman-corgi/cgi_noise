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
try:
    import unitsConstants as uc
except ModuleNotFoundError:
    from . import unitsConstants as uc
import math
try:
    from loadCSVrow import loadCSVrow
except ModuleNotFoundError:
    from .loadCSVrow import loadCSVrow
from dataclasses import dataclass, asdict
import numpy as np

def open_folder(*folders):
    """
    Constructs a path to a subfolder and returns a dictionary of file paths within it.

    Args:
        *folders: A sequence of folder names to be joined to the current working directory.

    Returns:
        A dictionary where keys are filenames (str) and values are Path objects
        for each file in the specified directory.
    """
    filenamedir = Path(os.getcwd())
    folder = Path(filenamedir, *folders)
    return {file.name: file for file in folder.iterdir() if file.is_file()}

def getScenFileNames(config):
    """
    Retrieves a list of full file paths for scenario-specific CSV data files.

    The function constructs paths based on filenames specified in the 'config'
    dictionary, organized under subfolders within an "Data" directory.

    Args:
        config: A dictionary containing a 'DataSpecification' key, which itself is a
                dictionary mapping descriptive keys (e.g., 'CoronagraphFile')
                to base filenames (without .csv extension).

    Returns:
        A list of strings, where each string is the absolute path to a CSV file.
        The order of files in the list corresponds to:
        1. CoronagraphFile
        2. QE_Curve_file
        3. DetModelFile_CBE
        4. StrayLightFRNfile
        5. ThroughputFile
        6. CalibrationFile
        7. ContrastStabilityFile
    """
    filenamedir = Path(os.getcwd())
    filenameList = []
    ffList = [
        ("Photometry", "Coronagraph_Data"),
        ('Photometry', 'QE_Curve_Data'),
        ('Photometry', 'Detector_Data'),
        ('Photometry', 'StrayLight_Data'),
        ('Photometry', 'Throughput_Data'),
        ('Calibration', 'Calibration_Data'),
        ('Cstability', 'ContrastStability_Data')
    ]
    for folder, key in ffList:
        name = config['DataSpecification'][key] + ".csv"
        path = filenamedir / "Data" / folder / name
        filenameList.append(str(path))
    return filenameList


def loadCSVs(filenameList):
    return [loadCSVrow(f) for f in filenameList]


def workingAnglePars(CG_Data, CS_Data):
    """
    Determines the effective Inner Working Angle (IWA) and Outer Working Angle (OWA).

    These are derived from coronagraph (CG_Data) and contrast stability (CS_Data)
    data. The IWA is the maximum of the individual IWAs, and the OWA is the
    minimum of the individual OWAs. Angles are typically in units of lambda/D.

    Args:
        CG_Data: Loaded CSV data for coronagraph performance. Expected to have a
                 DataFrame `df` with an 'r_lam_D' column.
        CS_Data: Loaded CSV data for contrast stability. Expected to have a
                 DataFrame `df` with an 'r_lam_D' column.

    Returns:
        A tuple (IWA, OWA), where:
            IWA (float): The maximum effective inner working angle (lambda/D).
            OWA (float): The minimum effective outer working angle (lambda/D).
    """
    IWAc = CG_Data.df.at[0, 'r_lam_D']
    IWAs = CS_Data.df.at[0, 'r_lam_D']
    
    if len(CG_Data.df) > 1:
        OWAc = CG_Data.df['r_lam_D'].iloc[-1]
    else:
        OWAc = CG_Data.df['r_lam_D'].iloc[0] + 3
    
    if len(CS_Data.df) > 1:
        OWAs = CS_Data.df['r_lam_D'].iloc[-1]
    else:
        OWAs = CS_Data.df['r_lam_D'].iloc[0] + 3
    
    # OWAs = CS_Data.df['r_lam_D'].iloc[-1]
    return max(IWAs, IWAc), min(OWAs, OWAc)


def contrastStabilityPars(CS_Type, planetWA, CS_Data):
    """
    Extracts contrast stability parameters from CSV data at a given planet working angle.

    Parameters are scaled by `uc.ppb` (parts per billion).

    Args:
        CS_Type: A string prefix (e.g., "MCBE_") used to identify relevant
                  column names in the CS_Data DataFrame.
        planetWA: The planet's working angle in units of lambda/D.
        CS_Data: Loaded CSV data for contrast stability. Expected to have a
                 DataFrame `df` with an 'r_lam_D' column and other columns
                 prefixed by `CS_Type`.

    Returns:
        A tuple containing:
            selDeltaC (float): Selected delta contrast (quadrature sum of stabilities).
            rawContrast (float): Average raw contrast at the planetWA.
            SystematicCont (float): Systematic contrast contribution.
            initStatRawContrast (float): Initial static raw contrast.
            rawContrast (float): (Repeated) Average raw contrast.
            IntContStab (float): Internal contrast stability.
            ExtContStab (float): External contrast stability.
        All returned contrast values are in parts per billion (ppb).

    Raises:
        IndexError: If the contrast stability file format is not as expected
                    (based on column count and names).
    """

    tol = 0.05
    indCS = CS_Data.df['r_lam_D'].searchsorted(planetWA + tol) - 1

    headers = CS_Data.df.columns.tolist()
    nCols = len(headers)
    fnARC = CS_Type + "_AvgRawContrast"
    fnECS = CS_Type + "_ExtContStab"
    fnICS = CS_Type + "_IntContStab"
    fnSC  = CS_Type + "_SystematicC"
    fnISRC = CS_Type + "_InitStatContrast"

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


def getFocalPlaneAttributes(opMode, config, DET_CBE_Data, lam, bandWidth, DPM, CGdesignWL, omegaPSF):
    FocalPlaneAtt = loadCSVrow(Path(os.getcwd(), 'Data', 'CONST_SNR_FPattributes.csv'))
    AmiciPar = loadCSVrow(Path(os.getcwd(), 'Data', 'CONST_Amici_parameters.csv'))
    """
    Calculates focal plane attributes based on the operational mode (Imaging or Spectroscopy).

    Loads 'CONST_SNR_FPattributes.csv' and 'CONST_Amici_parameters.csv' internally.

    Args:
        opMode: Operational mode, either "SPEC" (Spectroscopy) or "IMG" (Imaging).
        config: Scenario configuration dictionary, used for 'R_required' in SPEC mode.
        DET_CBE_Data: Loaded CSV data for detector model (Current Best Estimate - CBE).
                      Expected to have a DataFrame `df` with 'PixelSize_m'.
        lam: Observation wavelength in meters.
        bandWidth: Fractional spectral bandwidth (delta_lambda / lambda).
        DPM: Diameter of the primary mirror in meters.
        CGdesignWL: Coronagraph design wavelength in meters.
        omegaPSF: Solid angle of the PSF core in arcsec^2 (used in IMG mode).

    Returns:
        A tuple containing:
            f_SR (float): Spectral resolution factor (1/ (R * bandWidth) for SPEC, 1 for IMG).
            CritLam (float): Critical wavelength (Nyquist sampling) in meters.
            detPixSize_m (float): Detector pixel size in meters.
            mpix (float): Number of pixels in the photometric aperture (PSF core region).
            pixPlateSc (float): Pixel plate scale in milliarcseconds/pixel.

    Raises:
        Exception: If `opMode` is not "SPEC" or "IMG".
        KeyError: If required keys are missing in `config` for SPEC mode.
    """

    detPixSize_m = DET_CBE_Data.df.at[0, 'PixelSize_m']

    if opMode == "SPEC":
        try:
            resolution = config['instrument']['R_required'] #scenarioData.at['R_required', 'Latest']
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
    """
    Holds parameters related to coronagraph performance.

    Attributes:
        CGcoreThruput (float): Core throughput of the coronagraph system.
        PSFpeakI (float): Peak intensity of the coronagraphic PSF (normalized).
        omegaPSF (float): Solid angle of the PSF core in arcsec^2.
        CGintSamp (float): Sampling interval in lambda/D units from coronagraph data.
        CGradius_arcsec (float): Radius associated with the working angle in arcseconds.
        CGdesignWL (float): Design wavelength of the coronagraph in meters.
        CGintmpix (float): Number of pixels within the integration area,
                           derived from omegaPSF and sampling.
        CG_PSFarea_sqlamD (float): Area of the PSF core in (lambda/D)^2.
        CGintensity (float): Intensity value from coronagraph data at the working angle.
        CG_occulter_transmission (float): Transmission of the coronagraph occulter.
        CGcontrast (float): Raw contrast achieved by the coronagraph at the working angle.
        CGtauPol (float): Polarization throughput factor (default is 1.0).
    """
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
    """
    Represents an astrophysical target (star-planet system).

    Attributes:
        v_mag (float): V-band magnitude of the host star.
        dist_pc (float): Distance to the system in parsecs.
        specType (str): Spectral type of the host star (e.g., 'G2V').
        phaseAng_deg (float): Orbital phase angle of the planet in degrees.
        sma_AU (float): Semi-major axis of the planet's orbit in Astronomical Units.
        radius_Rjup (float): Radius of the planet in Jupiter radii.
        geomAlb_ag (float): Geometric albedo of the planet (unused if `albedo` is set).
        exoZodi (float): Level of exo-zodiacal dust, in units of "zodis"
                         (multiples of Solar System zodi brightness).
        albedo (float, optional): Effective Lambertian albedo of the planet.
                                  If None, can be calculated from flux ratio.
    """
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
        """Converts orbital parameters to on-sky separation.

        Args:
            sma_AU: Semi-major axis in AU.
            dist_pc: Distance to target in parsecs.
            phaseAng_deg: Planet's orbital phase angle in degrees.

        Returns:
            Projected separation in milliarcseconds (mas).
        """
        sep_mas = ((sma_AU * uc.AU * math.sin(math.radians(phaseAng_deg))) / (dist_pc * uc.pc)) / uc.mas
        return sep_mas
    
    @staticmethod
    def albedo_from_geomAlbedo(phaseAng_deg,geomAlb_ag):
        alpha = phaseAng_deg * uc.deg
        return 1/np.pi * (np.sin(alpha) + (np.pi - alpha)*np.cos(alpha))*geomAlb_ag
    
    @staticmethod
    def fluxRatio_to_deltaMag(fluxRatio):
        """Converts a flux ratio to a difference in magnitudes.

        Args:
            fluxRatio: Ratio of planet flux to star flux.

        Returns:
            Difference in magnitudes (delta_mag).
        """
        return (-2.5) * math.log10(fluxRatio)

    @staticmethod
    def deltaMag_to_fluxRatio(deltaMag):
        """Converts a difference in magnitudes to a flux ratio.

        Args:
            deltaMag: Difference in magnitudes.

        Returns:
            Ratio of fluxes.
        """
        return 10 ** (-0.4 * deltaMag)

    @staticmethod
    def fluxRatio_SMA_rad_to_albedo(fluxRatio, sma_AU, radius_Rjup):
        """Calculates planet albedo from flux ratio, SMA, and radius.

        Assumes Lambertian sphere phase function at 90 deg phase angle (phi(alpha)=1/pi)
        is implicitly handled if fluxRatio is defined appropriately.
        More generally, flux_ratio = albedo * (R_planet / SMA)^2 * phase_function.
        This function solves for albedo assuming phase_function is incorporated or is 1.

        Args:
            fluxRatio: Planet-to-star flux ratio.
            sma_AU: Semi-major axis in AU.
            radius_Rjup: Planet radius in Jupiter radii.

        Returns:
            The geometric albedo (Ag) if flux ratio definition is appropriate,
            or an effective Lambertian albedo.
        """
        return fluxRatio * (sma_AU * uc.AU / (radius_Rjup * uc.jupiterRadius)) ** 2


def coronagraphParameters(cg_df, config, planetWA, DPM):
    """
    Extracts and calculates coronagraph parameters for a given working angle.

    Args:
        cg_df: Pandas DataFrame containing coronagraph performance data, indexed
               by working angle ('r_lam_D'). Expected columns include 'coreThruput',
               'PSFpeak', 'area_sq_arcsec', 'r_as', 'I', 'occTrans', 'contrast'.
        planetWA: The planet's working angle in lambda/D.
        DPM: Diameter of the primary mirror in meters.

    Returns:
        A CGParameters dataclass instance populated with calculated values.
    """
 
    indWA = cg_df[(cg_df.r_lam_D <= planetWA)]['r_lam_D'].idxmax(axis=0)

        
    omegaPSF = cg_df.loc[indWA, 'area_sq_arcsec']
    CGintSamp = cg_df.loc[2, 'r_lam_D'] - cg_df.loc[1, 'r_lam_D']
    CGradius_arcsec = cg_df.at[indWA, 'r_as']

    CGdesignWL = DPM * cg_df.iloc[0, 1] * uc.arcsec / cg_df.iloc[0, 0]
    CGintmpix = omegaPSF * (uc.arcsec**2) / ((CGintSamp * CGdesignWL / DPM)**2)
    CG_PSFarea_sqlamD = omegaPSF / (CGdesignWL / uc.arcsec)**2

    CGintensity = cg_df.loc[indWA, 'I']
    CG_occulter_transmission = cg_df.at[indWA, 'occTrans']  
    CGcontrast = cg_df.loc[indWA, 'contrast']

    ObservationType = config['DataSpecification']['ObservationCase']
    if ObservationType.find("IMG_NFB1_HLC") != -1:
        # for Kappa_c, Core Throughput use TVAC measurement based on HLC Band 1
        Kappa_c_meas = config['TVACmeasured']['Kappa_c_HLCB1']
        CoreThroughput = config['TVACmeasured']['CoreThput_HLCB1']
        PSFpeakI = CoreThroughput * Kappa_c_meas / CGintmpix        
    else:        
        CoreThroughput = cg_df.loc[indWA, 'coreThruput'] 
        PSFpeakI = cg_df.loc[indWA, 'PSFpeak']  

    return CGParameters(
        CGcoreThruput = CoreThroughput,
        PSFpeakI = PSFpeakI,
        omegaPSF=omegaPSF,
        CGintSamp=CGintSamp,
        CGradius_arcsec=CGradius_arcsec,
        CGdesignWL=CGdesignWL,
        CGintmpix=CGintmpix,
        CG_PSFarea_sqlamD=CG_PSFarea_sqlamD,
        CGintensity=CGintensity,
        CG_occulter_transmission=CG_occulter_transmission,
        CGcontrast=CGcontrast
    )

def getSpectra(target, lam, bandWidth):
    """
    Calculates stellar flux based on spectral type, V-magnitude, wavelength, and bandwidth.

    Loads 'SPECTRA_ALL_BPGS.csv' which contains spectral data for various star types.

    Args:
        target: A Target dataclass instance with `specType` and `v_mag`.
        lam: Central observation wavelength in meters.
        bandWidth: Fractional bandwidth (delta_lambda/lambda).

    Returns:
        A tuple containing:
            inBandFlux0_sum (pd.Series): Integrated zero-magnitude flux (ph/s/m^2)
                                         over the band for all spectral types in the CSV.
            inBandZeroMagFlux (float): Integrated zero-magnitude flux (ph/s/m^2)
                                       for the target's spectral type.
            starFlux (float): Actual flux from the target star (ph/s/m^2) at Earth.
    """
    spectra_path = Path(os.getcwd(), 'Data', 'Spectra', 'SPECTRA_ALL_BPGS.csv')
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

def getStrayLightfromfile(scenario,perfLevel,STRAY_FRN_Data):
    # rowID = STRAY_FRN_Data.df.loc[STRAY_FRN_Data.df['PerfLevel']==perfLevel].index[0]
    rowID = STRAY_FRN_Data.df.loc[STRAY_FRN_Data.df['PerfLevel']==perfLevel].index[0]
    #scenario = scenarioData.at['Scenario','Latest']
    try:
        strayLight = STRAY_FRN_Data.df.at[rowID,scenario]
    except:
        #print(scenario)
        scenario = scenario.replace('DRM','EB')
        #print(scenario)
        try:
            strayLight = STRAY_FRN_Data.df.at[rowID,scenario]
        except:
            raise Exception('Stray needs help')
            strayLight = None
    #print(f'strayLight = {strayLight}')
    return strayLight


@dataclass
class DetNoiseRates:
    """
    Holds various detector noise rates.

    Attributes:
        dark (float): Dark current rate in electrons/pixel/second.
        CIC (float): Clock-Induced Charge rate in electrons/pixel/second.
        read (float): Read noise squared rate in (electrons/second)^2 / (pixel integrated area / s),
                      effectively (electrons_rms/frame)^2 * (pixels_in_aperture / frame_time).
                      This term is often directly added to variance.
    """
    dark: float
    CIC: float
    read: float

def detector_noise_rates(DET_CBE_Data, monthsAtL2, frameTime, mpix, isPhotonCounting):
    """
    Calculates detector noise rates considering mission lifetime degradation.

    Args:
        DET_CBE_Data: Loaded CSV data for detector model. Expected DataFrame `df`
                      with columns like 'DetEOL_mos', 'DarkBOM_e_per_pix_per_hr', etc.
        monthsAtL2: Mission duration at L2 in months, for degradation calculation.
        frameTime: Exposure time per frame in seconds.
        mpix: Number of pixels in the photometric aperture.
        isPhotonCounting: Boolean, True if operating in photon counting mode.

    Returns:
        A DetNoiseRates dataclass instance with calculated dark current, CIC,
        and read noise rates.
    """
    missionFraction = monthsAtL2 / DET_CBE_Data.df.at[0, 'DetEOL_mos']
    detDarkBOM = DET_CBE_Data.df.at[0, 'DarkBOM_e_per_pix_per_hr']
    detDarkEOM = DET_CBE_Data.df.at[0, 'DarkEOM_e_per_pix_per_hr']
    dark_per_hr = detDarkBOM + missionFraction * (detDarkEOM - detDarkBOM)
    dark_per_s = mpix * dark_per_hr / 3600  # from per pixel to per core

    detCIC1 = DET_CBE_Data.df.at[0, 'CICatGain1BOM_e_per_pix_per_fr']
    detCIC2 = DET_CBE_Data.df.at[0, 'CICatGain2BOM_e_per_pix_per_fr']
    gain1 = DET_CBE_Data.df.at[0, 'Gain1BOM']
    gain2 = DET_CBE_Data.df.at[0, 'Gain2BOM']
    CIC_degradation = DET_CBE_Data.df.at[0, 'CICdegradationEOM']
    EMgain = DET_CBE_Data.df.at[0, 'EMGain']

    CIC_rate = ((detCIC2 - detCIC1) / (gain2 - gain1)) * EMgain + (
        detCIC1 - ((detCIC2 - detCIC1) / (gain2 - gain1)) * gain1
    ) * (1 + missionFraction * (CIC_degradation - 1))
    CIC_per_s = mpix * CIC_rate / frameTime   # from per pixel to per core

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
    """
    Holds various optical throughput components.

    Attributes:
        refl (float): Combined reflectivity/transmissivity of telescope optics (OTA, CGI).
        filt (float): Filter throughput.
        polr (float): Polarizer throughput.
        core (float): Coronagraph core throughput (e.g., for planet light).
        occt (float): Coronagraph occulter transmission (e.g., for starlight/zodi).
    """
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

    Raises:
        ValueError: If `ezdistrib` is not a recognized value.
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
    speckleThroughput = thput.refl * thput.filt * thput.polr 
    locZodiThroughput = thput.refl * thput.filt * thput.polr * thput.occt
    exoZodiThroughput = locZodiThroughput * distFactor

    return thput, {
        "planet": planetThroughput,
        "speckle": speckleThroughput,
        "local_zodi": locZodiThroughput,
        "exo_zodi": exoZodiThroughput
    }


def rdi_noise_penalty(target, inBandFlux0_sum, starFlux, TimeonRefStar_tRef_per_tTar, RefStarSpecType, RefStarVmag):
    """
    Computes noise penalty factors due to Reference Differential Imaging (RDI).

    These factors (k_sp, k_det, k_lzo, k_ezo) quantify the increase in variance
    for different noise components when using RDI.

    Args:
        target: Target dataclass instance for the science target.
        inBandFlux0_sum: Pandas Series of integrated zero-magnitude flux (ph/s/m^2)
                         indexed by spectral type.
        starFlux: Flux of the science target star in ph/s/m^2.
        TimeonRefStar_tRef_per_tTar: Ratio of time spent on reference star to
                                     time spent on science target (t_ref / t_target).
        RefStarSpecType: Spectral type of the reference star (default: 'a0v').
        RefStarDist: Distance to the reference star in parsecs (default: 10 pc).
        RefStarVmag: V-band magnitude of the reference star (default: 3.0).

    Returns:
        A dictionary with RDI penalty factors:
            k_sp (float): Penalty for speckle noise.
            k_det (float): Penalty for detector noise.
            k_lzo (float): Penalty for local zodiacal light noise.
            k_ezo (float): Penalty for exo-zodiacal light noise.
    """

    RefStarinBandZeroMagFlux = inBandFlux0_sum.at[RefStarSpecType]
    RefStarFlux = RefStarinBandZeroMagFlux * (10 ** (-0.4 * RefStarVmag))
    BrightnessRatio = RefStarFlux / starFlux

    timeRatio = TimeonRefStar_tRef_per_tTar
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
    lam, mpix, cphrate_total):
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

    QE_img = QE_Data.df.loc[QE_Data.df['lambda_nm'] <= (lam / uc.nm), 'QE_at_neg100degC'].iloc[-1]
    det_EMgain = DET_CBE_Data.df.at[0, 'EMGain']
    det_readnoise = DET_CBE_Data.df.at[0, 'ReadNoise_e']
    det_PCthresh = DET_CBE_Data.df.at[0, 'PCThresh_nsigma']
    det_FWCser = DET_CBE_Data.df.at[0, 'FWCserial']
    det_FWCimg = DET_CBE_Data.df.at[0, 'FWCimage']
    CTI_clocking = DET_CBE_Data.df.at[0, 'CTI_clk']
    nTransfers = DET_CBE_Data.df.at[0, 'CTI_xfers']

    CTE = (1-CTI_clocking)**nTransfers

    if isPhotonCounting:
        ENF = 1.0
        effReadnoise = 0.0
        frameTime = round(min(tfmax, max(tfmin, desiredRate / (cphrate_total * QE_img / mpix))), 1)
        approxPerPixelPerFrame = frameTime * cphrate_total * QE_img / mpix
        eff_coincidence = (1 - math.exp(-approxPerPixelPerFrame)) / approxPerPixelPerFrame if approxPerPixelPerFrame > 0 else 1.0
        eff_thresholding = math.exp(-det_PCthresh * det_readnoise / det_EMgain)
        dQE = QE_img * eff_coincidence * eff_thresholding * CTE
    else:
        ENF = math.sqrt(2)
        effReadnoise = det_readnoise / det_EMgain
        Nsigma = 3
        NEE = Nsigma * ENF * det_EMgain
        y_crit = ((NEE**2 + 2 * det_FWCser) - math.sqrt(NEE**4 + 4 * NEE**2 * det_FWCser)) / 2
        tfr_crit = y_crit / (cphrate_total * QE_img / mpix)
        frameTime = min(tfmax, max(tfmin, math.floor(tfr_crit)))
        dQE = QE_img * CTE

    return ENF, effReadnoise, frameTime, dQE, QE_img


@dataclass
class VarianceRates:
    planet: float
    speckle: float
    locZodi: float
    exoZodi: float
    straylt: float
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

def noiseVarianceRates(cphrate, QE, dQE, ENF, detNoiseRate, k_sp, k_det, k_lzo, k_ezo,
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
        cg.PSFpeakI * cg.CGintmpix * speckleThroughput * Acol * dQE
    )

    # Note: the planet rate, following the EB model, differs from the others in using QE instead of dQE
    # this is the image area QE instead of delivered QE: this practice will need to be revisited
    rates = VarianceRates(
        planet  = f_SR * ENF**2 * cphrate.planet  * QE,  
        speckle = f_SR * ENF**2 * cphrate.speckle * dQE * k_sp,
        locZodi = f_SR * ENF**2 * cphrate.locZodi * dQE * k_lzo,
        exoZodi = f_SR * ENF**2 * cphrate.exoZodi * dQE * k_ezo,
        straylt = f_SR * ENF**2 * cphrate.straylt * dQE * k_det,
        detDark = f_SR * ENF**2 * detNoiseRate.dark     * k_det,
        detCIC  = f_SR * ENF**2 * detNoiseRate.CIC      * k_det,
        detRead = f_SR *          detNoiseRate.read     * k_det,
    )

    return rates, residSpecRate


def compute_tsnr(SNRdesired, planetSignalRate, nvRatesCore, residSpecRate):
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

    denom = planetSignalRate**2 - SNRdesired**2 * residSpecRate**2
    if denom <= 0:
        raise ValueError("compute_tsnr: SNR condition is not achievable with given residual speckle level.")

    timeToSNR = SNRdesired**2 * nvRatesCore.total / denom
    criticalSNR = planetSignalRate  / residSpecRate

    return timeToSNR, criticalSNR
