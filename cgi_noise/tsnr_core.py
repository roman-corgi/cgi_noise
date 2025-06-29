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

def tsnr_pipeline(config, DATA_DIR, target_params, verbose=True):
    ObservationCase = config['DataSpecification']['ObservationCase']

    DPM = config['instrument']['Diam']
    lam = config['instrument']['wavelength']
    lamD = lam / DPM
    intTimeDutyFactor = config['instrument']['dutyFactor']
    opMode = config['instrument']['OpMode']
    bandWidth = config['instrument']['bandwidth']

    target = fl.Target(**target_params)
    sep_mas = target.phaseAng_to_sep(target.sma_AU, target.dist_pc, target.phaseAng_deg)
    target.albedo = target.albedo_from_geomAlbedo(target.phaseAng_deg, target.geomAlb_ag)

    filenameList = fl.getScenFileNames(config, DATA_DIR)
    CG_Data, QE_Data, DET_CBE_Data, STRAY_FRN_Data, THPT_Data, CAL_Data, CS_Data = fl.loadCSVs(filenameList)
    CS_Type = config['DataSpecification']['CS_Type']

    IWA, OWA = fl.workingAnglePars(CG_Data, CS_Data)
    planetWA = sep_mas * uc.mas / lamD

    tol = 0.05
    if (IWA - tol) <= planetWA <= IWA:
        planetWA = IWA
    elif OWA <= planetWA <= (OWA + tol):
        planetWA = OWA
    elif planetWA < (IWA - tol) or planetWA > (OWA + tol):
        raise ValueError(f"Planet WA={planetWA:.1f} outside of IWA={IWA:.1f}, OWA={OWA:.1f}.")

    selDeltaC, AvgRawC, SystematicC, initStatRaw, IntContStab, ExtContStab = fl.contrastStabilityPars(CS_Type, planetWA, CS_Data)

    cg = fl.coronagraphParameters(CG_Data.df, config, planetWA, DPM)
    f_SR, _, detPixSize_m, mpix = fl.getFocalPlaneAttributes(opMode, config, DET_CBE_Data, lam, bandWidth, DPM, cg.CGdesignWL, cg.omegaPSF, DATA_DIR)

    inBandFlux0_sum, inBandZeroMagFlux, starFlux = fl.getSpectra(target, lam, bandWidth, DATA_DIR)

    TimeonRefStar_tRef_per_tTar = 0.25

    fluxRatio = target.albedo * (target.radius_Rjup * uc.jupiterRadius / (target.sma_AU * uc.AU))**2
    planetFlux = fluxRatio * starFlux

    magLocalZodi = config['instrument']['LocZodi_magas2']
    magExoZodi_1AU = config['instrument']['ExoZodi_magas2']
    absMag = target.v_mag - 5 * math.log10(target.dist_pc / 10)
    locZodiAngFlux = inBandZeroMagFlux * 10**(-0.4 * magLocalZodi)
    exoZodiAngFlux = inBandZeroMagFlux * 10**(-0.4 * (absMag - uc.sunAbsMag + magExoZodi_1AU)) / target.sma_AU**2 * target.exoZodi

    thput, throughput_rates = fl.compute_throughputs(THPT_Data, cg, "uniform")
    Acol = (np.pi / 4) * DPM**2
    stray_ph_s_mm2 = fl.getStrayLightfromfile(ObservationCase, 'CBE', STRAY_FRN_Data)
    stray_ph_s_pix = stray_ph_s_mm2 * (1 / uc.mm**2) * detPixSize_m**2

    cphrate = corePhotonRates(
        planet=planetFlux * throughput_rates["planet"] * Acol,
        speckle=starFlux * AvgRawC * cg.PSFpeakI * cg.CGintmpix * throughput_rates["speckle"] * Acol,
        locZodi=locZodiAngFlux * cg.omegaPSF * throughput_rates["local_zodi"] * Acol,
        exoZodi=exoZodiAngFlux * cg.omegaPSF * throughput_rates["exo_zodi"] * Acol,
        straylt=stray_ph_s_pix * mpix
    )
    cphrate.total = sum([cphrate.planet, cphrate.speckle, cphrate.locZodi, cphrate.exoZodi, cphrate.straylt])

    ENF, effReadnoise, frameTime, dQE, QE_img = fl.compute_frame_time_and_dqe(0.1, 3, 100, True, QE_Data, DET_CBE_Data, lam, mpix, cphrate.total)

    detNoiseRate = fl.detector_noise_rates(DET_CBE_Data, 21, frameTime, mpix, True)

    rdi_penalty = fl.rdi_noise_penalty(inBandFlux0_sum, starFlux, TimeonRefStar_tRef_per_tTar, 'a0v', 2.26)
    k_sp = rdi_penalty['k_sp']
    k_det = rdi_penalty['k_det']
    k_lzo = rdi_penalty['k_lzo']
    k_ezo = rdi_penalty['k_ezo']

    nvRatesCore, residSpecRate = fl.noiseRates(
        cphrate, QE_img, dQE, ENF, detNoiseRate,
        k_sp, k_det, k_lzo, k_ezo,
        f_SR, starFlux, selDeltaC,
        config['instrument']['pp_Factor_CBE'], cg,
        throughput_rates['speckle'], Acol
    )

    planetSignalRate = f_SR * cphrate.planet * dQE

    if verbose:
        print(f"Central wavelength: {lam / uc.nm:.1f} nm, with {bandWidth * 100:.0f}% BW\n")
        print(f"Lambda/D: {lamD / uc.mas:.3f} mas")
        print(f"Separation: {sep_mas:.0f} mas")
        print(f"Albedo: {target.albedo:.3f}")

        table = PrettyTable()
        table.field_names = ['planet WA', 'phase', 'dist', 'sma', 'sep', 'lam/D', "IWA", "OWA"]
        table.add_row([f'{planetWA:.2f}', f'{target.phaseAng_deg:.2f}', f"{target.dist_pc:.2f}", f'{target.sma_AU:.2f}', f'{sep_mas:.2f}', f'{lamD/uc.mas:.2f}', f'{IWA:.2f}', f'{OWA:.2f}'])
        print(table)

        print(f"Star Flux = {starFlux:.3e} ph/s/m^2")
        print(f"Planet Flux Ratio = {fluxRatio:.2e}\nPlanet Flux = {planetFlux:.3f} ph/s/m^2")
        print(f"Calculated Frame Time: {frameTime:.2f} s")
        print(f'QE in the image area: {QE_img:.3f}')
        print(f"Detected Quantum Efficiency (dQE): {dQE:.3f}")
        print(f"Excess Noise Factor (ENF): {ENF:.2f}")
        print(f"Core fraction used in the SNR region for mode {ObservationCase}: f_SR: {f_SR:.3f}")

        print("\nTotal noise variance rate breakdown:")
        table = PrettyTable()
        table.field_names = ['planet', 'speckle', 'local Zodi', 'exo Zodi', 'Stray']
        table.add_row([f"{nvRatesCore.planet:.4f}",f"{nvRatesCore.speckle:.4f}",f"{nvRatesCore.locZodi:.4f}",f"{nvRatesCore.exoZodi:.4f}",f"{nvRatesCore.straylt:.3e}",])
        print(table)

        print("\nContrast Stability Numbers:") 
        csfilename = None
        for filepath in filenameList:
            base = os.path.basename(filepath)
            if base.startswith("CS_"):
                csfilename = os.path.splitext(base)[0]
                break

        table = PrettyTable()
        table.field_names = ['CS case', 'DeltaC',  'External CS', 'Internal CS', 'Systematic', 'Avg Raw', 'initStatRaw']
        table.add_row([f'{csfilename}', f'{selDeltaC:.2e}', f'{ExtContStab:.2e}', f'{IntContStab:.2e}', f'{SystematicC:.2e}', f'{AvgRawC:.2e}', f"{initStatRaw:.2e}"])
        print(table)

    return planetSignalRate, nvRatesCore, residSpecRate
