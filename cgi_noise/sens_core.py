import math
import numpy as np
from dataclasses import dataclass
import cgi_noise.cginoiselib as fl
import cgi_noise.unitsConstants as uc
from prettytable import PrettyTable
import os

@dataclass
class corePhotonRates:
    planet: float
    speckle: float
    locZodi: float
    exoZodi: float
    straylt: float
    total: float = 0.0

def sens_pipeline(config, DATA_DIR, target_params, verbose=True, totalTinteg=1000, SNR=5.0):
    ObservationCase = config['DataSpecification']['ObservationCase']

    DPM = config['instrument']['Diam']
    lam = config['instrument']['wavelength']
    lamD = lam / DPM
    intTimeDutyFactor = config['instrument']['dutyFactor']
    usableTinteg = intTimeDutyFactor * totalTinteg
  
    opMode = config['instrument']['OpMode']
    bandWidth = config['instrument']['bandwidth']

    target = fl.Target(**target_params)

    filenameList = fl.getScenFileNames(config, DATA_DIR)
    CG_Data, QE_Data, DET_CBE_Data, STRAY_FRN_Data, THPT_Data, CAL_Data, CS_Data = fl.loadCSVs(filenameList)
    CS_Type = config['DataSpecification']['CS_Type']

    inBandFlux0_sum, inBandZeroMagFlux, starFlux = fl.getSpectra(target, lam, bandWidth, DATA_DIR)
    magLocalZodi = config['instrument']['LocZodi_magas2']
    magExoZodi_1AU = config['instrument']['ExoZodi_magas2']
    absMag = target.v_mag - 5 * math.log10(target.dist_pc / 10)
    locZodiAngFlux = inBandZeroMagFlux * 10**(-0.4 * magLocalZodi)
    TimeonRefStar_tRef_per_tTar = 0.25
    
    IWA, OWA = fl.workingAnglePars(CG_Data, CS_Data)

    Acol = (np.pi / 4) * DPM**2    
    npts = 50
    WAset = np.linspace(IWA, OWA, npts)
    Sensitivity = np.zeros(npts)
    # dC = np.zeros(npts)
    sep = np.zeros(npts)
    # rsr = np.zeros(npts)
    KC = np.zeros(npts)
    tauPk = np.zeros(npts)
    intMpix = np.zeros(npts)
    tauC = np.zeros(npts)
    dC= np.zeros(npts)
    coreArea = np.zeros(npts)
    for iWA in range(npts):
        WA = WAset[iWA]
        selDeltaC, AvgRawC, SystematicC, initStatRaw, IntContStab, ExtContStab = fl.contrastStabilityPars(CS_Type, WA, CS_Data)

        cg = fl.coronagraphParameters(CG_Data.df, config, WA, DPM)
        f_SR, _, detPixSize_m, mpix = fl.getFocalPlaneAttributes(opMode, config, DET_CBE_Data, lam, bandWidth, DPM, cg.CGdesignWL, cg.omegaPSF, DATA_DIR)

        sep[iWA] = lamD * WA
        
        sma_AU = target.dist_pc * uc.pc * sep[iWA] / np.sin(target.phaseAng_deg *uc.deg) / uc.AU

        exoZodiAngFlux = inBandZeroMagFlux * 10**(-0.4 * (absMag - uc.sunAbsMag + magExoZodi_1AU)) / sma_AU**2 * target.exoZodi

        thruputComponents, throughput_rates = fl.compute_throughputs(THPT_Data, cg, "uniform")
    
        stray_ph_s_mm2 = fl.getStrayLightfromfile(ObservationCase, 'CBE', STRAY_FRN_Data)
        stray_ph_s_pix = stray_ph_s_mm2 * (1 / uc.mm**2) * detPixSize_m**2

        cphrate = corePhotonRates(
            planet=0,
            speckle=starFlux * AvgRawC * cg.PSFpeakI * cg.CGintmpix * throughput_rates["speckle"] * Acol,
            locZodi=locZodiAngFlux * cg.omegaPSF * throughput_rates["local_zodi"] * Acol,
            exoZodi=exoZodiAngFlux * cg.omegaPSF * throughput_rates["exo_zodi"] * Acol,
            straylt=stray_ph_s_pix * mpix
        )
        cphrate.total = sum([cphrate.planet, cphrate.speckle, cphrate.locZodi, cphrate.exoZodi, cphrate.straylt])

        # print(f'Cavg = {AvgRawC:.2e}, Thp =  {throughput_rates["speckle"]:.2e} ')
        ENF, effReadnoise, frameTime, dQE, QE_img = fl.compute_frame_time_and_dqe(0.1, 3, 100, True, QE_Data, DET_CBE_Data, lam, mpix, cphrate.total)

        detNoiseRate = fl.detector_noise_rates(DET_CBE_Data, 21, frameTime, mpix, True)
    
        rdi_penalty = fl.rdi_noise_penalty(inBandFlux0_sum, starFlux, TimeonRefStar_tRef_per_tTar, 'a0v', 2.26)
        k_sp  = rdi_penalty['k_sp']
        k_det = rdi_penalty['k_det']
        k_lzo = rdi_penalty['k_lzo']
        k_ezo = rdi_penalty['k_ezo']

        randomNoiseVarianceRates, residSpecStdDevRate = fl.noiseRates(
            cphrate, QE_img, dQE, ENF, detNoiseRate,
            k_sp, k_det, k_lzo, k_ezo,
            f_SR, starFlux, selDeltaC,
            config['instrument']['pp_Factor_CBE'], cg,
            throughput_rates['speckle'], Acol
        )
        randomNV = randomNoiseVarianceRates.total * usableTinteg
        resSpecNV = (residSpecStdDevRate * usableTinteg)**2
        nonPlanetNoiseVariance = randomNV + resSpecNV
        # nvr = randomNoiseVarianceRates
        # print(f'{iWA} WA = {WA:.3f}:  Rnd tot: {randomNV:.2f}  spec:{(nvr.speckle*usableTinteg):.2f}  exo:{(nvr.exoZodi*usableTinteg):.2f}   resid Spec NV: {resSpecNV:.2f}')
        Kappa = 1.0 / (f_SR * starFlux * Acol * throughput_rates["planet"] * dQE * usableTinteg)
        KC[iWA] = (cg.PSFpeakI * cg.CGintmpix) / thruputComponents.core
        Sensitivity[iWA] = (Kappa/2.0) * SNR**2 * ( 1 + np.sqrt( 1 + 4.0* nonPlanetNoiseVariance/SNR**2))
        
        dC[iWA]=selDeltaC 
        tauPk[iWA] = cg.PSFpeakI
        intMpix[iWA] = cg.CGintmpix
        tauC[iWA] = thruputComponents.core
        coreArea[iWA] = cg.CG_PSFarea_sqlamD
    return WAset, sep,  Sensitivity, KC, tauPk, intMpix, tauC, coreArea, dC