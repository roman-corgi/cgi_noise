# Library of Functions Used in Various Jupyter Notebooks for CGI Perf
import os
import sys
import pandas as pd
import math
current_dir = os.getcwd()
from loadCSVrow import loadCSVrow
from loadXLcol import loadXLcol
from pathlib import Path
from datetime import datetime as dt
import xlwings as xw
import unitsConstants as uc
import numpy as np

def open_folder(*folders):
    filenamedir = current_dir
    folder = Path(filenamedir, *folders)
    filesDict = {}
    
    for file in folder.iterdir():
        filesDict[file.name]=(file)    
    return filesDict

#------------------------------------------------------------------------------
class Target:
    """Takes in these parameters v_mag, dist_pc, exoZodi_xsolar, specType,
    Calculates these parameters
    why to give:
    optional parameters:
    name
    phase angle 
    used in Ref Star, """
    
    def __init__(self, v_mag = None, phaseAng_deg = None, geomAlb_ag = None,\
                           radius_Rjup = None, dist_pc = None, exoZodi = None,\
                              specType = None, sma_AU = None, albedo = None,\
                                  eeid_mas = None):
        
        self.v_mag = v_mag
        self.phaseAng_deg = phaseAng_deg
        self.geomAlb_ag = geomAlb_ag
        self.radius_Rjup = radius_Rjup
        self.dist_pc = dist_pc
        self.exoZodi = exoZodi
        self.specType = specType 
        self.sma_AU = sma_AU
        self.albedo = albedo
        self.eeid_mas = eeid_mas
        
    @staticmethod   
    def albedo_for_planetSensitivity(geomAlb_ag, phaseAng_deg):
        """ Calculate albedo from geometric albedo and phase function"""
        phaseAng_rad = math.radians(phaseAng_deg)  
        phaseFunction = (1/math.pi)*(math.sin(phaseAng_rad) +\
                                                   (math.pi-phaseAng_rad)*\
                                                       math.cos(phaseAng_rad))
        albedo = geomAlb_ag*phaseFunction
        
        return albedo
   
    @staticmethod
    def fluxRatio_alb_sma_to_Radius(fluxRatio,albedo,sma_AU):
        """Calculate the radius of the planet from the SMA, 
        flux ratio, and albedo"""
        
        radius_Rjup = sma_AU*uc.AU*(math.sqrt(fluxRatio/albedo))/uc.jupiterRadius
        
        return radius_Rjup

    @staticmethod    
    def alb_rad_sma_to_fluxRatio(albedo,radius_Rjup,sma_AU):
        """Calculate fluxRatio from planet albedo, planet radius (uc.jupiterRadius) and 
            class attribute, semimajor axis (sma_AU).
            planetFluxRatio on SNR page"""
        
        fluxRatio = albedo * (radius_Rjup * uc.jupiterRadius / (sma_AU * uc.AU) )**2
        
        return fluxRatio


    @staticmethod    
    def planetWAtoSMA(planetWA, lamD, dist_pc, phaseAng_deg):
        """Calculate planet SMA from planet WA"""

        sma_AU = planetWA * lamD * dist_pc*(uc.pc/uc.AU)/\
            math.sin(math.radians(phaseAng_deg))
            
        return sma_AU
    
    @staticmethod            
    def planetWAtoSeparation(planetWA, lamD):
        """Calculate angular separation in milliarcseconds from 
        planet working angle in lambda/D"""
        sep_mas = planetWA * lamD / uc.mas
        
        return sep_mas
    
    @staticmethod
    def sep_and_dist_to_SMA(sep_mas,dist_pc,phaseAng_deg):
        """Calculate semi-major axis of elliptical orbit (SMA) in AU 
        from the angular separation in milliarcseconds, distance to target in 
        seconds of parallax, and phase angle in degrees"""
        
        sma_AU = float((dist_pc * uc.pc * sep_mas * uc.mas))/\
                     (math.sin(math.radians(phaseAng_deg)))/uc.AU
        return sma_AU
                        

    @staticmethod
    def fluxRatio_to_deltaMag(fluxRatio):        
        """Calculate delta Mag from flux ratio"""
        deltaMag = (-2.5)*math.log10(fluxRatio)
        return deltaMag
    
    @staticmethod
    def deltaMag_to_fluxRatio(deltaMag):
        """Calculate flux ratio from given delta mag"""
        fluxRatio = 10**((-0.4)*deltaMag)
        return fluxRatio
    
    @staticmethod    
    def fluxRatio_SMA_rad_to_albedo(fluxRatio,sma_AU,radius_Rjup):
        """Calculate albedo from radius, flux ratio, and SMA"""
        albedo = fluxRatio *(sma_AU*uc.AU/(radius_Rjup*uc.jupiterRadius))**2
        return albedo
        
    @staticmethod    
    def separation_to_planetWA(sep_mas, lamD, IWA, OWA):
        """Calculate planet working angle from separation in milliarcseconds"""
        planetWA = sep_mas*uc.mas / lamD
        tolerance = .05
        
        if (IWA - tolerance) <= planetWA <= IWA:
            planetWA = IWA
        elif OWA <= planetWA <= (OWA + tolerance):
            planetWA = OWA
        elif planetWA < (IWA - tolerance) or planetWA > (OWA + tolerance):
            raise ValueError(f" Planet WA={planetWA:.1f} while IWA = {IWA:.1f} and OWA = {OWA:.1f} lam/D.")
                
        return planetWA
        

    @staticmethod    
    def phaseAng_to_sep(sma_AU, dist_pc, phaseAng_deg):
        """Calculate separation(mas) from planet phase angle(deg),
            and class attribute dist_pc (distance in parsecs) 
            and class attribute sma_AU (semimajor axis in AU)"""
              
        sep_mas = ((sma_AU * uc.AU * math.sin\
                         (math.radians(phaseAng_deg)))\
                        /(dist_pc * uc.pc))/uc.mas
        return sep_mas
   
    
    @staticmethod
    def deltaMag_from_known_planet_location(dist_pc, phaseAng_deg, sep_mas,\
                                            geomAlb_ag, radius_Rjup):
        """Calculate delta mag from planet phase angle (deg),
        geometric albedo (ag), and planet radius (uc.jupiterRadius)   
        assuming the insta-planet's location
        Use this to know what delta-mag to plug into the
        time-to-SNR calculator. It provides a way to calculate
        delta mag given the other parameters."""
        
        phAng_rad = math.radians(phaseAng_deg)
        phasefunction = (1/math.pi)*(math.sin(phAng_rad)+(math.pi-phAng_rad)*math.cos(phAng_rad))
        sma_AU = dist_pc*uc.pc*sep_mas*uc.mas/math.sin(phAng_rad)/uc.AU
        
        #calculate flux ratio assuming lambertian scatter
        fluxRatio_lambertian = geomAlb_ag * phasefunction * (radius_Rjup*uc.jupiterRadius/(sma_AU*uc.AU))**2
        deltaMag = (-2.5)*math.log10(fluxRatio_lambertian)
        
        return deltaMag
      
#-----------------  END OF TARGET CLASS  --------------------------------------
#------------------------------------------------------------------------------



def getScenFileNames_DRM(scenarioData):
    filenamedir = current_dir #Path(os.getcwd()).parent.parent
    filenameList= []
    ffList = [("Photometry", "CoronagraphFile"),('Photometry','QE_Curve_file'),\
                      ('Photometry','DetModelFile_CBE'),('Photometry','StrayLightFRNfile'),\
                      ('Photometry','ThroughputFile'),('Calibration','CalibrationFile'),\
                       ('Cstability','ContrastStabilityFile')]
    for item in ffList:
        filename = os.path.abspath(os.path.join(filenamedir,"EBcsvData",item[0],\
                                scenarioData.loc[item[1],'Latest']+".csv"))        
        filenameList.append(filename)        
    
    return filenameList


def loadCSVs(filenameList):
    loadedFiles = []
    for filename in filenameList:
        loadedFiles += [loadCSVrow(filename)]
    return loadedFiles

def workingAnglePars(CG_Data, CS_Data):
    
    IWAc = CG_Data.df.at[0,'rlamD']
    IWAs = CS_Data.df.at[0,'r_lam_D']
    OWAc = CG_Data.df['rlamD'].iloc[-1] # idx of last data row is just before first comment row
    OWAs = CS_Data.df['r_lam_D'].iloc[-1] # idx of last data row is just before first comment row

    return max(IWAs, IWAc) , min(OWAs, OWAc) 

def contrastStabilityPars( CSprefix, planetWA, CS_Data):
    """Constrast Stability"""

    # pick an index into the C-Stability Dataframe, and pin to nearest working angle available if outside range  
    tol = .05
    indCS = np.searchsorted(CS_Data.df['r_lam_D'], planetWA + tol) - 1

    headers = CS_Data.df.columns.tolist()
    nCols = len(headers)
    fnARC = CSprefix + "AvgRawContrast"
    fnECS = CSprefix + "ExtContStab"
    fnICS = CSprefix + "IntContStab"
    fnSC  = CSprefix + "SystematicC"
    fnISRC = CSprefix + "InitStatContrast"
    
    ExtContStab = CS_Data.df.at[indCS, fnECS] 
    IntContStab = CS_Data.df.at[indCS, fnICS]
    
    rawContrast = CS_Data.df.at[indCS, fnARC] #same as average Raw contrast SelContrast in EB spreadsheet
    initStatRawContrast = CS_Data.df.at[indCS, fnISRC] 
    
    selContrast = rawContrast
    
    if nCols == 16 and 'SystematicC' in headers[13]:
        SystematicCont = CS_Data.df.at[indCS, fnSC]
        selDeltaC = math.sqrt((ExtContStab**2) + (IntContStab**2) + (SystematicCont**2))*uc.ppb
    elif nCols == 13:
        SystematicCont = 0
        selDeltaC = math.sqrt((ExtContStab**2) + (IntContStab**2))*uc.ppb
    else:
        raise IndexError('The contrast stability file referenced is not formatted as expected.')

    return selDeltaC, selContrast, SystematicCont, initStatRawContrast,\
        rawContrast, IntContStab, ExtContStab


def getDarkCurrent(DET_CBE_Data, monthsatL2):
    
        # mission End of Life
    detEOL_mos = DET_CBE_Data.df.at[0,'DetEOL_mos']
    missionFraction = monthsatL2 /detEOL_mos
        
    # Detector dark current

    detDarkBOM = DET_CBE_Data.df.at[0,'DarkBOM_e_per_pix_per_hr']
    detDarkEOM = DET_CBE_Data.df.at[0,'DarkEOM_e_per_pix_per_hr']
    DarkCur_epoch_per_hr = detDarkBOM+missionFraction*(detDarkEOM-detDarkBOM)
    DarkCur_epoch_per_s = DarkCur_epoch_per_hr/3600

    return detDarkBOM, detDarkEOM, DarkCur_epoch_per_s, DarkCur_epoch_per_hr,\
        missionFraction, detEOL_mos

def coronagraph_pars(CG_Data, planetWA, IWA, DPM, lamD):
    
    # #Coronagraph Parameters ---need planet WA and IWA and OWA
    
    CorforPol = 0;
    CGtauPol = 1 + CorforPol;
       
    indWA = CG_Data.df[(CG_Data.df.rlamD <=planetWA)]['rlamD'].idxmax(axis=0)# index of last radial slice < planet's radius
     
    # Frac light incident on primary in FWHM of PSF centered at r
    modelCoreThput = CG_Data.df.loc[indWA,'coreThruput']*CGtauPol
    
    # Frac light in peak pixel of PSF centered at r & normalized to flux
    PSFpeakI = CG_Data.df.loc[indWA,'PSFpeak']*CGtauPol
    
    # Area of FWHM of PSF centered at r
    omegaPSF = CG_Data.df.loc[indWA,'area_sq_arcsec']
      
    # radius used in CGperf as best match to target working angle
    CGradius_arcsec = CG_Data.df.at[indWA,'r_as']
    
    CGdesignWL = DPM * CG_Data.df.iloc[0,1] * uc.arcsec / CG_Data.df.iloc[0,0]
    
    # sampling used by John Krist in generating CG results (intrinsic sampling)
    CGintSamp  = CG_Data.df.loc[2,'rlamD']-CG_Data.df.loc[1,'rlamD'] 
    
    # pixels within PSF core assuming instrinsic sampling
    CGintmpix = omegaPSF * (uc.arcsec**2)/((CGintSamp*CGdesignWL/DPM)**2) #mpixIntrinsic
    
    CG_PSFarea_sqlamD = omegaPSF/(lamD/uc.arcsec)**2
    
    CGintmpix = omegaPSF * (uc.arcsec**2)/((CGintSamp*CGdesignWL/DPM)**2) #mpixIntrinsic
     
    # Speckle intensity at r (normalized to flux)
    CGintensity = CG_Data.df.loc[indWA,'I']
    
    # Occulter*lyot stop transmission at r
    CG_occulter_transmission = CG_Data.df.at[indWA, 'occTrans']*CGtauPol
   
    #  Intensity/PSF_Peak
    CGcontrast = CG_Data.df.loc[indWA,'contrast']
    
    return CGtauPol, indWA, modelCoreThput, PSFpeakI, omegaPSF, CGintSamp,\
        CGradius_arcsec, CGdesignWL, CGintmpix, CG_PSFarea_sqlamD, CGintensity,\
            CG_occulter_transmission, CGcontrast

def getFocalPlaneAttributes(opMode, scenarioData, DET_CBE_Data, lam, bandWidth, DPM, CGdesignWL, omegaPSF):
    filenamedir = current_dir
        
    FocalPlaneAtt = loadCSVrow(Path(filenamedir,'EBcsvData','CONST_SNR_FPattributes.csv'))
        
    AmiciPar = loadCSVrow(Path(filenamedir,'EBcsvData','CONST_Amici_parameters.csv'))
        
    detPixSize_m = DET_CBE_Data.df.at[0,'PixelSize_m']
    
    ## Get Focal Plane Attributes for Selected Operating Mode
    
    if opMode == "HLC_NF_IMG":
        f_SR = 1 
        CritLam = FocalPlaneAtt.df.at[0,'Critical_Lambda_m']
        mpix = omegaPSF * uc.arcsec**2 * (lam*uc.nm/CGdesignWL)**2*(2*DPM/CritLam)**2 
        pixPlateSc = CritLam/DPM/2/uc.mas 
        
    elif opMode == "SPC_Amici_SPEC":   
        try:
            resolution = scenarioData.at['R_required','Latest']
            print(f"resolution = {resolution}")
            f_SR = 1/(resolution*bandWidth)        
        except:
            print("R_required is not specified in scenario-- set to default 0")
            resolution = 0.0001
            f_SR = -1 
            
        CritLam = FocalPlaneAtt.df.at[1,'Critical_Lambda_m'] 
        mpix = mpix_Amici(AmiciPar, lam, DPM,detPixSize_m,resolution)
        pixPlateSc = CritLam/DPM/2/uc.mas 
            
    elif opMode == "SPC_WF_IMG":
        f_SR = 1 
        CritLam = FocalPlaneAtt.df.at[2,'Critical_Lambda_m']
        mpix = omegaPSF*uc.arcsec**2*(lam*uc.nm/CGdesignWL)**2*(2*DPM/CritLam)**2 
        pixPlateSc = CritLam/DPM/2/uc.mas
        
    else:
        raise Exception("Check Scenario Operating Mode in Scenario xlsx file and Operational Parameters readme")
    
    return f_SR, CritLam, detPixSize_m, mpix, pixPlateSc
          
def getSpectra( target, lam, bandWidth):
    filenamedir = current_dir
    SpectraFolder = os.path.join(filenamedir,'EBcsvData','Spectra')
    
    #### Using hardcoded csv filename
    SPECTRA_file = os.path.abspath(os.path.join(SpectraFolder,"SPECTRA_ALL_BPGS.csv"))
    SPECTRA_Data = loadCSVrow(SPECTRA_file)# Spectra data csv to dataframe
    
    # rows of the spectra within the request band
    bandRange = SPECTRA_Data.df[abs(SPECTRA_Data.df['Wavelength_m'] - lam)  <= (0.5*bandWidth*lam)]
        
    # temporarily remove the Wavelength column, 
    # since we want to divide the different spectral types by EPhot
    onlySpec = bandRange.drop(['Wavelength_m', 'E_ph_J'], axis=1)
        
    # Calculate scalar Ephot
    Ephot = uc.h_planck * uc.c_light / (lam)
    
    # Divide each element in each column by Ephot
    onlySpecEphot = onlySpec.apply(lambda x:x/Ephot, axis =1, result_type ='broadcast')
    
    # Calculate the increments of wavelength in the table
    deltaLambda = SPECTRA_Data.df.at[2,'Wavelength_m'] -SPECTRA_Data.df.at[1,'Wavelength_m']
    
    # sum of each of spectral types over multiple wavelengths within band
    inBandFlux0_sum = (onlySpecEphot.sum(axis=0)) * deltaLambda
    
    # Get spectrum for star to find inBWflux0
    inBandZeroMagFlux = inBandFlux0_sum.at[target.specType]
    
    starFlux = inBandZeroMagFlux * 10**((-0.4)*target.v_mag) # ph/s/m^2
      
    return inBandFlux0_sum, inBandZeroMagFlux, starFlux

def getFluxRatio(target, starFlux):
      
    fluxRatio = Target.alb_rad_sma_to_fluxRatio(target.albedo,\
                                                target.radius_Rjup,\
                                                    target.sma_AU)
    planetFlux = fluxRatio * starFlux
    
    return fluxRatio, planetFlux

def getrefStarRDI(target, inBandFlux0_sum, starFlux,\
                  RefStarSpecType, RefStarVmag,\
                      RefStarDist, RefStarExoZodi, timeOnRef):
    
    # Observing scenario is RDI with one observation of a brighter reference
    # and one observation of the target. The noise is computed effective after the
    # speckle subtraction. To the observation 
    
    RefStarAbsMag = RefStarVmag - 5*math.log10(RefStarDist/10)
    RefStarinBandZeroMagFlux = inBandFlux0_sum.at[RefStarSpecType]
    RefStarDeltaMag = target.v_mag - RefStarVmag
    RefStarFlux = RefStarinBandZeroMagFlux*(10**((-0.4)*RefStarVmag))
    BrightnessRatio = RefStarFlux/starFlux
    betaRDI = 1 / (BrightnessRatio*timeOnRef)
    
    k_sp = 1 + betaRDI
    k_det = 1 + betaRDI**2 * timeOnRef
    k_lzo = k_det
    k_ezo = k_sp
            
    v_sp = math.sqrt(betaRDI)
    v_det = betaRDI * math.sqrt(timeOnRef)
    v_lzo = betaRDI * math.sqrt(timeOnRef)
    v_ezo = math.sqrt(betaRDI)
       
    return RefStarVmag, RefStarAbsMag, RefStarinBandZeroMagFlux, RefStarDeltaMag,\
         RefStarFlux, BrightnessRatio, betaRDI, k_sp, k_det, k_lzo, k_ezo,\
             v_sp, v_det, v_lzo, v_ezo
                 
def getZodi(magLocalZodi, magExoZodi_1AU, target, inBandZeroMagFlux, omegaPSF):
    # Zodi brightness - assume local zodi = 22.1 mag /arcsec**2 in V band
    # assume  skynoise = (local-zodi + exo-zodi)
    absMag = target.v_mag - 5 * math.log10(target.dist_pc / 10)
               
    # zodi solid angular flux in units of ph/s/m**2/arcsec**2
    loZodiAngFlux =  inBandZeroMagFlux *10**((-0.4)* magLocalZodi) 
    exoZodiAngFlux = target.exoZodi * inBandZeroMagFlux * \
                    10**(-0.4*(absMag-uc.sunAbsMag+magExoZodi_1AU)) / (target.sma_AU**2) 
            
    ZodiFlux    = (loZodiAngFlux + exoZodiAngFlux) * omegaPSF 
    loZodiFlux  = loZodiAngFlux * omegaPSF 
    exoZodiFlux = exoZodiAngFlux * omegaPSF     
   
    return ZodiFlux, exoZodiAngFlux, loZodiAngFlux, loZodiFlux, exoZodiFlux, absMag

def getNoiseRates(f_SR, starFlux, fluxRatio, colArea, thpt_t_pnt, thpt_tau_pk, thpt_t_speckle,\
             det_QE, exoZodiFlux, loZodiFlux, thpt_t_unif, k_comp, rawContrast,\
                 CGintmpix, mpix, DarkCur_epoch_per_s):
    
    rate_planet_imgArea = f_SR * starFlux * fluxRatio * colArea * thpt_t_pnt * det_QE
    
    rate_exoZodi_incPht = f_SR * exoZodiFlux * colArea * thpt_t_unif
    
    rate_loZodi_incPht =  f_SR * loZodiFlux * colArea * thpt_t_unif 
      
    rate_Zodi_imgArea = (rate_exoZodi_incPht + rate_loZodi_incPht)*det_QE

    # speckleRate_imgArea
    rate_speckleBkg = f_SR * starFlux * rawContrast * uc.ppb * thpt_tau_pk * \
        CGintmpix * colArea * thpt_t_speckle  * det_QE 
        
    # total rate pixel with planet
    rate_photoConverted = DarkCur_epoch_per_s +\
        (rate_planet_imgArea + rate_Zodi_imgArea + rate_speckleBkg )/mpix
    
    rate_totalwithoutplanet =  DarkCur_epoch_per_s + (rate_Zodi_imgArea +\
                                                      rate_speckleBkg )/mpix

    return rate_planet_imgArea, rate_Zodi_imgArea, rate_exoZodi_incPht,\
        rate_loZodi_incPht, rate_speckleBkg, rate_photoConverted,\
            rate_totalwithoutplanet
        

def getFrameExposureTime(DET_CBE_Data, FWC_gr, rate_totalwithoutplanet,\
                         rate_photoConverted, isPhotonCounting):
    
    # Calculate Frame Time for CIC in dark current units and CR hits per frame
    # This is for Photon Counting
    
    detEMgain = DET_CBE_Data.df.at[0,'EMGain'] # Electron multiplication gain

    # maximum frame exposure time for analog
    maxANLGt_fr= 0.9*FWC_gr/(3*detEMgain*rate_totalwithoutplanet)
    
    frameTime_ANLG = min(maxANLGt_fr,100)
    
    # calculated frame exposure time for photon counting opt frame
    maxPCt_fr = round(min(80, max(3, 0.1/rate_photoConverted)))
    
    if isPhotonCounting:    
        frameTime = maxPCt_fr
    else: # analog
        frameTime = frameTime_ANLG
    
    return frameTime, frameTime_ANLG,maxANLGt_fr,maxPCt_fr, detEMgain

def getDetectorCIC(DET_CBE_Data, detEMgain, missionFraction, frameTime):
    
    # ## Instrument Parameters
    # # Detector

    detCICdegradationEOM = DET_CBE_Data.df.at[0,'CICdegradationEOM']
    detCIC1 = DET_CBE_Data.df.at[0,'CICatGain1BOM_e_per_pix_per_fr']
    detCIC2 = DET_CBE_Data.df.at[0,'CICatGain2BOM_e_per_pix_per_fr']
    detCICgain1 = DET_CBE_Data.df.at[0,'Gain1BOM']
    detCICgain2 = DET_CBE_Data.df.at[0,'Gain2BOM']
    
    det_CIC_gain = ((detCIC2 - detCIC1)/(detCICgain2 - detCICgain1))*detEMgain\
        +(detCIC1-((detCIC2 - detCIC1)/(detCICgain2 - detCICgain1))*detCICgain1)\
            *(1 + missionFraction*(detCICdegradationEOM-1))

    det_CIC_in_DC_units = det_CIC_gain/frameTime

    return det_CIC_in_DC_units,  det_CIC_gain

def getDetectorCosmicRays(perfLevel,DET_CBE_Data, detEMgain, frameTime ):
    
    # Cosmic Ray Tail Length
    detCRtailGain1 = DET_CBE_Data.df.at[0,'CRtailGain1']
    detCRtailGain2 = DET_CBE_Data.df.at[0,'CRtailGain2']
    detCRtailLen1 = DET_CBE_Data.df.at[0,'CRtailLen1_pix']
    detCRtailLen2 = DET_CBE_Data.df.at[0,'CRtailLen2_pix']
    CRtailLen_CBE = ((detCRtailLen2-detCRtailLen1)/(detCRtailGain2-detCRtailGain1))*detEMgain+(detCRtailLen2- ((detCRtailLen2-detCRtailLen1)/(detCRtailGain2-detCRtailGain1))*detCRtailGain2)
    
    CRtailLen_REQ = 150 #pixels -- Made up on 6/14/2021 -- ATTENTION!
    
    if perfLevel == "CBE":
        CRtailLen = CRtailLen_CBE
    else:
        CRtailLen = CRtailLen_REQ
        
    CRrate = 5*(100**2) # rays/s/m^2  
    detPixSize = DET_CBE_Data.df.at[0,'PixelSize_m']
    detPixAcross = DET_CBE_Data.df.at[0,'PixelsAcross_pix']
       
    CRhitsPerFrame = CRrate * frameTime * (detPixSize * detPixAcross)**2 # Cosmic Rays per frame
    
    return CRhitsPerFrame,detPixAcross,detPixSize,CRrate,CRtailLen

def getCTE(DET_CBE_Data, rate_photoConverted,frameTime,missionFraction):
    detFluxSlope = DET_CBE_Data.df.at[0,'CTE_FluxSlope_per_log_of_flux_per_e_per_pix_per_s']
    detFluxKnee = DET_CBE_Data.df.at[0,'CTE_dqeKneeFlux_e_per_pix_per_s']
    
    signalPerPixPerFrame = rate_photoConverted*frameTime

    CTE_traps = min(1,max(0,(detFluxSlope*missionFraction)\
                              *(math.log10(signalPerPixPerFrame)-math.log10(detFluxKnee)) +1))
       
    detGen = float(DET_CBE_Data.df.at[0,'Generation'])
    if detGen > 3:
        detCTEclk = DET_CBE_Data.df.at[0,'CTI_clk']
        detCTExfers = DET_CBE_Data.df.at[0,'CTI_xfers']
        CTE_clocking_efficiency = (1-detCTEclk)**detCTExfers
    else:
        CTE_clocking_efficiency = 1
    
    return CTE_clocking_efficiency, CTE_traps, signalPerPixPerFrame

def gethotPixels(DET_CBE_Data, missionFraction):
    
    hotPixFrac = DET_CBE_Data.df.at[0,'HotPixFrac']
    hotPix = hotPixFrac*missionFraction
    
    return hotPixFrac, hotPix

def getENF(isPhotonCounting):
    if isPhotonCounting:
        ENF = 1
    else:
        ENF = math.sqrt(2) # 1.414
    return ENF

def getReadNoiseandPCeffloss(detCamRead, detPCthreshold, isPhotonCounting, frameTime, detEMgain):
    """Read noise"""
    readNoise_w_gain = detCamRead/detEMgain # read noise with gain if analog
   
    if isPhotonCounting:
        readNoise = 0 # Minimal read noise with photon counting
    else:
        readNoise = readNoise_w_gain # read noise with gain if analog
    
    readNoise_leakage = 0.5*math.erfc((detPCthreshold/math.sqrt(2)))
  
    readNoise_leakage_in_current_units = readNoise_leakage/frameTime
    
    if isPhotonCounting:
        PCeffloss = 1 - math.exp( -detPCthreshold*detCamRead/detEMgain)
    else:
        PCeffloss = 0
   
    return readNoise, readNoise_leakage, readNoise_leakage_in_current_units,\
         PCeffloss, readNoise_w_gain

def getdetdQE(det_CTE, PCeffloss, hotPix, signalPerPixPerFrame, detPixAcross, CRtailLen, CRhitsPerFrame, det_QE):
    det_PC_threshold_efficiency = 1 - PCeffloss
   
    signal_region_electron_rate = signalPerPixPerFrame * det_CTE

    # Photon-counting coincidence efficiency 
    det_PC_coincid_effic = (1 - math.exp(-signal_region_electron_rate))\
        / (signal_region_electron_rate)
    
    # Efficiency after subtracting fraction lost to hot pixels
    det_hotPix = 1 - hotPix

    det_cosmicRays = 1 - CRhitsPerFrame * CRtailLen/(detPixAcross**2)
    
    estimated_dQE_CBE = det_QE * det_CTE * det_PC_threshold_efficiency * \
            det_PC_coincid_effic * det_hotPix * det_cosmicRays
    
    return signal_region_electron_rate, det_PC_threshold_efficiency,\
        det_PC_coincid_effic, det_hotPix, det_cosmicRays, estimated_dQE_CBE 

def getNoiseVarianceRatesI( f_SR, starFlux, planetFlux, colArea, rawContrast, \
              thpt_t_pnt, thpt_t_speckle, thpt_tau_pk, CGintmpix, k_pp,\
                     dQE, rate_exoZodi_incPht, rate_loZodi_incPht,\
                         selDeltaC):
    """rates for time to SNR"""
    planetRate_proc = f_SR*planetFlux*colArea*thpt_t_pnt*dQE

    speckleRate_proc = f_SR * starFlux * rawContrast * thpt_tau_pk\
         * CGintmpix * thpt_t_speckle * colArea * dQE * uc.ppb
    
    ezo_bkgRate = rate_exoZodi_incPht * dQE
    lzo_bkgRate = rate_loZodi_incPht * dQE
    zodiRate_proc = ezo_bkgRate + lzo_bkgRate
    
    residSpecRate =  f_SR * starFlux * (selDeltaC/k_pp) *\
         thpt_tau_pk * CGintmpix * thpt_t_speckle * colArea * dQE
            
    return planetRate_proc, speckleRate_proc, zodiRate_proc,\
            ezo_bkgRate, lzo_bkgRate, residSpecRate

# def intTime(dutyFactor, allocTinteg):
#     # return the available actual integration time, given a duty factor and the total allocated integration time for the observation
#     # the duty factor accounts for the fraction of the time the reaction wheel jitter meets requirements
 
#     usableTinteg =  dutyFactor * allocTinteg * uc.hour
    
#     return usableTinteg

def getStrayLightFRN(scenario, perfLevel, STRAY_FRN_Data, CG_Data, IWA, OWA,\
                     opMode, DET_CBE_Data, ENF, detPixSize, mpix, dQE, Kappa,\
                          usableTinteg, compStarSpecType, inBandFlux0_sum,\
                              f_SR, CGintmpix, colArea, thpt_t_speckle):
    try:
        strayLight = getStrayLightfromfile(scenario, perfLevel, STRAY_FRN_Data)
    except:
        photonRateFlux,_ = getStrayLight_luminescenceBackground(perfLevel,opMode,\
                                                              DET_CBE_Data,\
                                                                  mpix, dQE, False)
        
        leakageRateSL, leakageRate, strayLightBackground\
            = getStrayLight_companionStarBackground(compStarSpecType,\
                                                              IWA, OWA, CG_Data,\
                                          detPixSize, inBandFlux0_sum, f_SR,\
                                              CGintmpix, usableTinteg, mpix,\
                                              colArea, thpt_t_speckle, dQE, False)
        
        strayLight,_,_,_,_,_,_ = getStrayLightfromcalc(CG_Data, photonRateFlux,\
                                                         leakageRateSL, IWA, OWA)
    
    # This is the flux ratio noise from stray light:
    strayLight_FRN = Kappa * math.sqrt(ENF**2 * strayLight * 1000000 \
                                   * detPixSize**2 * mpix * dQE * usableTinteg) 

    strayLight_ph_pix_s = strayLight * 1000000 * detPixSize**2
    
    strayLight_ph_pix_h = strayLight_ph_pix_s * uc.hour
    
    strayLight_e_pix_h = strayLight_ph_pix_h * dQE
    
    return strayLight, strayLight_FRN, strayLight_ph_pix_s,\
        strayLight_ph_pix_h, strayLight_e_pix_h
       
def DRMinstaplanet(planetRate_proc,usableTinteg, totNoiseVarRate,\
                             residSpecRate, SNRtarget):
    """User-defined planet specifications and integration time
    for estimating time to SNR and
    possible SNR in allowed time"""
    signalCounts = planetRate_proc*usableTinteg
    varianceRandomNoise = totNoiseVarRate*usableTinteg 
    varianceResidSpeckle = (residSpecRate*usableTinteg)**2
    SNRinAllowedTime = signalCounts/math.sqrt(varianceRandomNoise+varianceResidSpeckle)
    
    timetoSNR = (SNRtarget**2)*totNoiseVarRate/\
        (planetRate_proc**2-SNRtarget**2*residSpecRate**2)/3600
        
    instaPlanetDF =\
        pd.DataFrame({'Name':[f'Time to reach SNR = {SNRtarget}',\
                              'usable integration', 'signal counts',\
                              'variance of random noise',\
                                  'variance of residual speckle',\
                                      'SNR in allowed time'],\
                      'Value':[timetoSNR,usableTinteg,signalCounts,\
                               varianceRandomNoise, varianceResidSpeckle,\
                                   SNRinAllowedTime],\
                          'Units':['hrs','sec','e-','e-','e-','SNR']})
    return instaPlanetDF

def getStrayLightfromfile(scenario,perfLevel,STRAY_FRN_Data):
    
    # rowID = STRAY_FRN_Data.df.loc[STRAY_FRN_Data.df['PerfLevel']==perfLevel].index[0]
    rowID = STRAY_FRN_Data.df.loc[STRAY_FRN_Data.df['PerfLevel']==perfLevel].index[0]
    
    try:
        strayLight = STRAY_FRN_Data.df.at[rowID,scenario]
    except:
        scenario = scenario.replace('DRM','EB')
        try:
            strayLight = STRAY_FRN_Data.df.at[rowID,scenario] 
        except:
            raise Exception('Stray needs help')
            strayLight = None
    return strayLight

def getStrayLight_luminescenceBackground(perfLevel,opMode,DET_CBE_Data,\
                                         mpix, dQE, seeOutput):
    filenamedir = current_dir
    sp = Path(filenamedir, 'EBcsvData', 'Photometry','STRAY_Assumptions.csv')
    MiscSTRAY = loadCSVrow(sp)

    GCRfluxL2 = pd.to_numeric(MiscSTRAY.df.at[0,'GCRfluxatL2_evts_per_sq_cm_per_sec']) #5
    photonsPerRelEvent = float(MiscSTRAY.df.at[0,'ph_per_rel_event_per_mm']) #250
    diamCompBeam = float(MiscSTRAY.df.at[0,'diamCompBeam_mm'])
    
    if opMode == "HLC_NF_IMG":
        luminOpticThickness = MiscSTRAY.df.at[0,'B1_LuminescOpticThickness_mm']
    elif opMode == "SPC_Amici_SPEC":   
        luminOpticThickness = MiscSTRAY.df.at[0,'B3_LuminescOpticThickness_mm']
    elif opMode == "SPC_WF_IMG":
        luminOpticThickness = MiscSTRAY.df.at[0,'B4_LuminescOpticThickness_mm']
    else:
        raise Exception("Check Scenario Operating Mode in Scenario csv file and Operational Parameters readme")
    
    luminOpticDistance = MiscSTRAY.df.at[0,'LuminescOpticDist_m'] #0.1
    sBaffling = MiscSTRAY.df.at[0,'s_baffling_requirement'] #0.001
    
    lumRatePerAngle = photonsPerRelEvent/(2*math.pi)
    luminOpticArea = (math.pi/4) * ((diamCompBeam/10)**2)
    
    detPixSize = DET_CBE_Data.df.at[0,'PixelSize_m']
    detPixAcross = DET_CBE_Data.df.at[0,'PixelsAcross_pix']
    
    omegaSignal = mpix * detPixSize**2 / luminOpticDistance**2 
    omegaIndirect = 2*math.pi * sBaffling * mpix/(detPixAcross**2)
    
    phRateDirect = GCRfluxL2 * lumRatePerAngle * luminOpticArea\
        * luminOpticThickness * omegaSignal 
        
    phRateIndirect = GCRfluxL2 * lumRatePerAngle * luminOpticArea\
        * luminOpticThickness * omegaIndirect 
    
    photonRate = phRateDirect + phRateIndirect 
    photonRatePerPix = photonRate / mpix 
    photonRateFlux = photonRatePerPix / detPixSize**2
    
    eRatePerSNR = photonRatePerPix * mpix * dQE
    
    if seeOutput:
        print(f'GCR flux at L2 =                      {GCRfluxL2:.2e}  events/cm^2/s')
        print(f'photons/relativistic event =          {photonsPerRelEvent:.2e}  ph/ev/mm/2piSr')
        print(f'diameter of compressed beam =         {diamCompBeam:.2e}  mm')
        print(f'luminesc. optic thickness =           {luminOpticThickness:.2e}  mm')
        print(f'luminesc. optic distance =            {luminOpticDistance:.2e}  m')
        print(f's_baffling (requirement) =            {sBaffling:.2e}')
        print(f'lum rate per solid angle =            {lumRatePerAngle:.2e}  ph/Sr/evt/mm')
        print(f'luminesc. optic area =                {luminOpticArea:.2e}  cm^2')
        print(f'omega_signal =                        {omegaSignal:.2e}  Sr')
        print(f'omega_indirect =                      {omegaIndirect:.2e}  Sr')
        print(f'stray light ph rate - direct =        {phRateDirect:.2e}  ph/s')
        print(f'stray light ph rate - indirect =      {phRateIndirect:.2e}  ph/s')
        print(f'stray light photon rate =             {photonRate:.2e}  ph/s')
        print(f'Luminscence Bkg phot rate per pix =   {photonRatePerPix:.2e}  ph/pix/s')
        print(f' expressed as phot flux on detector = {photonRateFlux:.2e}  ph/mm^2/s')
        
    
    return photonRateFlux,eRatePerSNR

def getStrayLight_companionStarBackground(compStarSpecType, IWA, OWA, CG_Data,\
                                          detPixSize, inBandFlux0_sum, f_SR,\
                                              CGintmpix,usableTinteg, mpix, \
                                              colArea, thpt_t_speckle, dQE,\
                                                  seeOutput):
    filenamedir = current_dir
    sp = Path(filenamedir, 'EBcsvData', 'Photometry','STRAY_Assumptions.csv')
    MiscSTRAY = loadCSVrow(sp)

    darkHoleMidRadius = (IWA + OWA) / 2 
    PSFpeakMidRadius = CG_Data.df.at[CG_Data.df[CG_Data.df['rlamD'] <= 6]['rlamD'].idxmax(),'PSFpeak'] 
    
    # companion star
    compStarLeak = MiscSTRAY.df.at[0,'CompnStarLeakage']
    compStarMag = MiscSTRAY.df.at[0,'CompnStarMag_mag']
    
    compStarIntegFlux = inBandFlux0_sum.at[compStarSpecType]
    compStarFluxInBand = 10**(-0.4 * compStarMag) * compStarIntegFlux 

    strayLightBackground = f_SR * compStarFluxInBand * compStarLeak * PSFpeakMidRadius\
    * CGintmpix * colArea * thpt_t_speckle * dQE * usableTinteg
    
    leakageRate = (strayLightBackground / usableTinteg) / mpix
    
    leakageRateSL = leakageRate / detPixSize**2
    
    if seeOutput:
        print('Companion Star Background')
        print(f'Mid-Radius in the dark hole =      {darkHoleMidRadius} lam/D')
        print(f"PSF peak at Mid-Radius =           {PSFpeakMidRadius} lam/D" )
        print(f'Companion Star Leakage =           {compStarLeak}')
        print(f'Companion Star Magnitude =         {compStarMag} mag')
        print(f'Companion Star Spectral Type =     {compStarSpecType}')
        print(f'Companion Star Integ Flux =        {compStarIntegFlux:.1e}')
        print(f'Companion Star Flux In Band =      {compStarFluxInBand:.1e}')
        
    
    return leakageRateSL, leakageRate, strayLightBackground

def getStrayLightfromcalc(CG_Data, photonRateFlux, leakageRateSL, IWA, OWA):
    # MUFs in the stray light calculation reflect uncertainties estimated in the scientific papers
    filenamedir = current_dir #Path(os.getcwd()).parent.parent
    thispath = Path(filenamedir, 'EBcsvData', 'Photometry','STRAY_MUF_table.csv')
    MUFtSTRAY = loadCSVrow(thispath) 
    
    #Light Pollution Requirements in native units for various contributors with MUF
    # Backup formulation of stray light FRN
    #1)Cherenkov + Fluorescence + Phosph.: 
    #(Cherenkov estimate only since it dominates. MUF also applied)
    luminBackground = photonRateFlux * MUFtSTRAY.df.at[0,'MUF']/1000000
    
    #2)Astronomical Sources:
    companionStar = leakageRateSL * MUFtSTRAY.df.at[1,'MUF']/1000000
    
    #3)Stray Light: (Indexed to the sum of the calculated numbers)
    strayLightwMUF = MUFtSTRAY.df.at[2,'Indx']*(luminBackground + companionStar)
    
    #4)Payload Internal Sources (also indexed)
    payloadInternal = MUFtSTRAY.df.at[3,'Indx']*(luminBackground + companionStar)
    
    #5)Rates Combined-- stray light reequirement
    strayLight = luminBackground + companionStar \
        + strayLightwMUF + payloadInternal
    
    #6)Internal Stray Light Sources: requirement
    intslreq = MUFtSTRAY.df.at[0,'CGI']*luminBackground\
        + MUFtSTRAY.df.at[1,'CGI']*companionStar \
        + MUFtSTRAY.df.at[2,'CGI']*strayLightwMUF \
        + MUFtSTRAY.df.at[3,'CGI']*payloadInternal
    
    #7)External Stray Light Sources: requirement
    extslreq = MUFtSTRAY.df.at[0,'Ext']*luminBackground\
        + MUFtSTRAY.df.at[1,'Ext']*companionStar\
        + MUFtSTRAY.df.at[2,'Ext']*strayLightwMUF\
        + MUFtSTRAY.df.at[3,'Ext']*payloadInternal
    
    return strayLight, luminBackground,companionStar,strayLightwMUF,\
        payloadInternal, intslreq, extslreq

def DRM_planetSens_Nsigma(residSpecRate, usableTinteg, ENF, k_sp, specRate_proc,\
                         k_lzo, lzo_bkgRate, k_ezo, ezo_bkgRate, k_det,\
                             darkNoiseRate, CIC_RNLK_noiseRate, readNoiseRate,\
                                SNRtarget, Kappa):
    """Calculations for planet FRN sensitivity"""
    res_Speckle = residSpecRate*usableTinteg
    
    nonPlanetVarRate = (ENF**2)*(k_sp*specRate_proc+k_lzo*lzo_bkgRate+\
                             k_ezo*ezo_bkgRate+k_det*\
                             (darkNoiseRate+CIC_RNLK_noiseRate))\
                        + k_det*readNoiseRate 
    nonpl_random = math.sqrt(nonPlanetVarRate*usableTinteg)
    tot_nonpl_noise = math.sqrt(nonpl_random**2 + res_Speckle**2)
    N_sigmaSens = (0.5*SNRtarget**2)*\
        (1+math.sqrt(1+4*tot_nonpl_noise**2/SNRtarget**2))*(Kappa*uc.ppb)
        
    return nonPlanetVarRate, nonpl_random, tot_nonpl_noise, N_sigmaSens

def DRM_planetSens_contribtoNsigmaSens(N_sigmaSens, tot_nonpl_noise, ENF,\
                                       k_sp, specRate_proc, zodiRate_proc,\
                                k_det, darkNoiseRate, CIC_RNLK_noiseRate,\
                                   readNoiseRate, usableTinteg, f_SR,\
                                     starFlux, colArea, thpt_t_pnt, dQE):
    # Contributions to N-sig sensitivity
    planet_implicit = math.sqrt((ENF**2)*f_SR*starFlux*\
                                N_sigmaSens*colArea*thpt_t_pnt*\
                                dQE*usableTinteg)
    
    speckle = math.sqrt(ENF**2*k_sp*specRate_proc*usableTinteg)
    zodi = math.sqrt(ENF**2*k_sp*zodiRate_proc*usableTinteg)
    dark = math.sqrt(ENF**2*k_det*darkNoiseRate*usableTinteg)
    cicrnleakage = math.sqrt(ENF**2*k_det*CIC_RNLK_noiseRate*usableTinteg)
    read = math.sqrt(k_det*readNoiseRate*usableTinteg)
    
    # non-planet random noise 
    nonplanetRandom = math.sqrt(speckle**2 + zodi**2 + dark**2 + cicrnleakage**2 + read**2)
    
    SNRcheck = planet_implicit**2/math.sqrt(tot_nonpl_noise**2 + planet_implicit**2)
    
    return nonplanetRandom, SNRcheck, planet_implicit, speckle, zodi, dark,\
        cicrnleakage, read

def DRM_planetSens_NsigmaSensInfiniteTime(k_det, readNoiseRate, usableTinteg,\
                                          nonpl_random, SNRtarget,\
                                              residSpecRate, Kappa):
    read = math.sqrt(k_det*readNoiseRate*usableTinteg)
    tot_nonpl_noise_inftime = math.sqrt(read**2+nonpl_random**2)
    infiniteTimeSens = SNRtarget*residSpecRate*(Kappa*uc.ppb)*usableTinteg
    
    return tot_nonpl_noise_inftime, infiniteTimeSens

def DRM_planetSens_calcFRN(ENF, k_sp, k_lzo, k_ezo, k_det, lzo_bkgRate, \
                           ezo_bkgRate, darkNoiseRate, CIC_RNLK_noiseRate,\
                               readNoiseRate, residSpecRate, specRate_proc,\
                                   usableTinteg, Kappa):
    res_Speckle = residSpecRate*usableTinteg
    
    randomNoiseRate = ENF**2*(k_sp*specRate_proc+k_lzo*lzo_bkgRate\
                              +k_ezo*ezo_bkgRate+\
                                  k_det*(darkNoiseRate+CIC_RNLK_noiseRate))\
                                    +k_det*readNoiseRate
    
    nonplrandom = math.sqrt(randomNoiseRate * usableTinteg)
    
    totnonplnoise = math.sqrt(nonplrandom**2 + res_Speckle**2)
    
    cstab = Kappa*res_Speckle
    
    bde = Kappa*nonplrandom
    
    total = Kappa*totnonplnoise
    
    return res_Speckle, randomNoiseRate, nonplrandom, totnonplnoise,\
        cstab, bde, total

def DRM_planetSens_vals(scenarioData, perfLevel, CSprefix, isPhotonCounting,\
                    target, planetWA, lamD, detPCthreshold):
    """Calculations for Planet Sensitivity noise sources 
    and contributors to sensitivity"""
    
    target.sma_AU = target.planetWAtoSMA(planetWA, lamD,\
                                         target.dist_pc, target.phaseAng_deg)
    
    planetRate_proc, usableTinteg, totNoiseVarRate, residSpecRate,\
        SNRtarget, ENF, k_sp, k_det, k_lzo, k_ezo, ezo_bkgRate, lzo_bkgRate,\
            specRate_proc, zodiRate_proc, Kappa, f_SR, starFlux, colArea, \
            darkNoiseRate,CIC_RNLK_noiseRate, readNoiseRate, thpt_t_pnt, dQE,\
                planetWA, lamD, IWA, OWA, CG_Data,thpt_t_refl,\
                    inBandZeroMagFlux,omegaPSF,allocTinteg, k_pp_CBE = \
    DRMgetSNRvals(scenarioData, perfLevel, CSprefix, isPhotonCounting, target, detPCthreshold)

    nonPlanetVarRate, nonpl_random, tot_nonpl_noise, N_sigmaSens=\
    DRM_planetSens_Nsigma(residSpecRate, usableTinteg, ENF, k_sp, specRate_proc,\
                             k_lzo, lzo_bkgRate, k_ezo, ezo_bkgRate, k_det,\
                                 darkNoiseRate, CIC_RNLK_noiseRate, readNoiseRate,\
                                    SNRtarget, Kappa)
    
    tot_nonpl_noise_inftime, infiniteTimeSens =\
    DRM_planetSens_NsigmaSensInfiniteTime(k_det, readNoiseRate, usableTinteg,\
                                             nonpl_random, SNRtarget,\
                                             residSpecRate, Kappa)
    
    nonplanetRandom, SNRcheck, planet_implicit, speckle, zodi,_,_,_ =\
    DRM_planetSens_contribtoNsigmaSens(N_sigmaSens, tot_nonpl_noise, ENF,\
                                       k_sp, specRate_proc, zodiRate_proc,\
                                k_det, darkNoiseRate, CIC_RNLK_noiseRate,\
                                   readNoiseRate, usableTinteg, f_SR,\
                                     starFlux, colArea, thpt_t_pnt, dQE)
    _,_,_,_,cstab, bde, total=\
    DRM_planetSens_calcFRN(ENF, k_sp, k_lzo, k_ezo, k_det, lzo_bkgRate, \
                           ezo_bkgRate, darkNoiseRate, CIC_RNLK_noiseRate,\
                               readNoiseRate, residSpecRate, specRate_proc,\
                                   usableTinteg, Kappa)

    separation, nsigmasens_finite, nsigmasens_inf =\
        DRM_tableofNsigma(N_sigmaSens, infiniteTimeSens, planetWA, lamD)
    
    return  nsigmasens_finite, nsigmasens_inf, separation,\
        planet_implicit, speckle, zodi, cstab, bde, total

def mpix_Amici(AmiciPar, lambda_nm, DPM, detPixSize_m, resolution):

    compbeamD_m = AmiciPar.df.loc[0,'compressd_beam_diamtr_m']# compressed beam diameter
    fnlFocLen_m = AmiciPar.df.loc[0,'final_focal_length_m']# final focal length
    Fno = fnlFocLen_m / compbeamD_m ; # Fno of final focus
    
    pixPerlamD = lambda_nm * uc.nm * Fno / detPixSize_m;# pixels per lam/D
    PSF_x_lamD = AmiciPar.df.loc[0,'PSF_core_x_extent_lamD']# PSF core x extent--invariant of wavelength
    PSF_y_lamD = AmiciPar.df.loc[0,'PSF_core_y_extent_lamD'] # PSF core y extent--invariant of wavelength
    xpixPerCor = PSF_x_lamD * 2 * pixPerlamD; # x direction pixels per core--w/o dispersion
    ypixPerCor = PSF_y_lamD * 2 * pixPerlamD; # y direction pixels per core--w/o dispersion
    
    Rlamsq = AmiciPar.df.loc[0,'lam_squared']
    Rlam = AmiciPar.df.loc[0,'lam']
    Rconst = AmiciPar.df.loc[0,'constant']
    
    ResPowatPSF = Rconst + Rlam * lambda_nm + Rlamsq * lambda_nm**2# Resolving power at PSF core width--Curve from Qian
    dpix_dlam = ResPowatPSF * xpixPerCor / lambda_nm# pix/nm Pixelwise dispersion dpix_dlam --Reverse engineered 
    
    xpixPerSpec = dpix_dlam * lambda_nm/resolution # x pixels per spectral element in dispersion (x) direction
    mpix = xpixPerSpec * ypixPerCor #pixels-- mpix (pixels per spec. elem) For the specified R value
    # mpixlim = xpixPerCor * ypixPerCor #pixels-- mpix core limited

    return mpix
    
def DRMgetSNRvals(scenarioData, perfLevel, CSprefix, isPhotonCounting, target, detPCthreshold, compStarSpecType='m5v'):
       
    lam, SNRtarget, scenario, allocTinteg, intTimeDutyFactor_CBE,\
        intTimeDutyFactor_REQ, k_pp_CBE, k_pp_REQ, monthsatL2, resolution,\
            bandWidth, thpt_t_core_REQ, tau_core_REQ,\
                opMode, fidPlnt, t_core_CBEmeas, K_c_CBEmeas \
                     = unpackScenPars(scenarioData)
    
    magLocalZodi , magExoZodi_1AU, DarkCurrent_adjust, CIC_adjust, QE_adjust,\
        SNR_for_NEFR, Channels_per_iter, Probes_per_Channel, Bandwidth_per_Channel,\
            DarkHole_Contrast, ComparisonTime_sec, FWC_gr, ENF_for_Analog,\
                DPM, thpt_t_PSFnominal, k_comp, timeOnRef,  RefStarSpecType,\
                     RefStarVmag_CBE, RefStarDist, RefStarExoZodi, RefStarVmag_req\
                         = unpackMiscPars()
    
    lamD = (lam*uc.nm)/DPM 
    
    sep_mas = Target.phaseAng_to_sep(target.sma_AU, target.dist_pc,\
                                     target.phaseAng_deg)
   
    #=-----------------------------------------
    # get list of specific csv files for selected mode
    filenameList = getScenFileNames_DRM(scenarioData)
    
    # get the data tables of the mode-specific csv files
    [CG_Data, QE_Data, DET_CBE_Data, STRAY_FRN_Data, THPT_Data, CAL_Data, CS_Data] \
    = loadCSVs(filenameList)
    
    # calculate planet working angle from given target specifications
    IWA, OWA = workingAnglePars(CG_Data, CS_Data)
    planetWA = Target.separation_to_planetWA(sep_mas, lamD, IWA, OWA)
    
    # Constrast Stability for Mission CBE
    selDeltaC, selContrast, SystematicCont, initStatRawContrast,\
        rawContrast, IntContStab, ExtContStab\
            = contrastStabilityPars( CSprefix, planetWA, CS_Data)
    
    det_QE = getQE(scenarioData, QE_Data)
    detDarkBOM, detDarkEOM, DarkCur_epoch_per_s, DarkCur_epoch_per_hr,\
        missionFraction, detEOL_mos\
        = getDarkCurrent(DET_CBE_Data, DarkCurrent_adjust, monthsatL2)
    
       
    CGtauPol, indWA, CGcoreThruput, PSFpeakI, omegaPSF, CGintSamp,\
        CGradius_arcsec, CGdesignWL, CGintmpix, CG_PSFarea_sqlamD, CGintensity,\
            CG_occulter_transmission, CGcontrast = \
                coronagraph_pars(CG_Data, planetWA, IWA, DPM, lamD)
                
    thpt_t_obsc, colArea, thpt_OTA_TCA, thpt_CGI, thpt_bkgLimitedCore,\
            thpt_coreArea_LamD, thpt_t_PSF, thpt_t_core,\
                 tau_core_REQ, thpt_t_occ, thpt_t_pnt, thpt_t_pol,\
                    thpt_t_refl, thpt_t_speckle, thpt_t_unif, thpt_tau_pk\
                        = throughput_pars(perfLevel, THPT_Data, scenarioData,\
                                          CG_occulter_transmission,\
                    CGcoreThruput, PSFpeakI, CGtauPol, omegaPSF, \
                        DPM, lamD)
    
    f_SR, CritLam, detPixSize_m, mpix, pixPlateSc =\
        getFocalPlaneAttributes(opMode, scenarioData, DET_CBE_Data, lam, bandWidth,\
                            DPM, CGdesignWL, omegaPSF)


    inBandFlux0_sum, inBandZeroMagFlux, starFlux =\
        getSpectra(target, lam, bandWidth)
    
    fluxRatio, planetFlux = getFluxRatio(target, starFlux)
    
    RefStarVmag, RefStarAbsMag, RefStarinBandZeroMagFlux, RefStarDeltaMag,\
         RefStarFlux, BrightnessRatio, betaRDI, k_sp, k_det, k_lzo, k_ezo,\
             v_sp, v_det, v_lzo, v_ezo =\
                 getrefStarRDI(target, inBandFlux0_sum, starFlux,\
                  RefStarSpecType, RefStarVmag_CBE,\
                      RefStarDist, RefStarExoZodi, timeOnRef)
    
    ZodiFlux, exoZodiAngFlux, loZodiAngFlux, loZodiFlux, exoZodiFlux, absMag = \
        getZodi(magLocalZodi, magExoZodi_1AU, target, inBandZeroMagFlux, omegaPSF)
    
    rate_planet_imgArea, rate_Zodi_imgArea, rate_exoZodi_incPht,\
        rate_loZodi_incPht, rate_speckleBkg, rate_photoConverted,\
            rate_totalwithoutplanet = \
    getNoiseRates(f_SR, starFlux, fluxRatio, colArea, thpt_t_pnt, thpt_tau_pk,\
                  thpt_t_speckle,\
             det_QE, exoZodiFlux, loZodiFlux, thpt_t_unif, k_comp, rawContrast,\
                 CGintmpix, mpix, DarkCur_epoch_per_s)
    
    frameTime, frameTime_ANLG,maxANLGt_fr,maxPCt_fr, detEMgain =\
    getFrameExposureTime(DET_CBE_Data, FWC_gr, rate_totalwithoutplanet,\
                         rate_photoConverted, isPhotonCounting)
        
    det_CIC_in_DC_units, det_CIC_at_epoch, det_CIC_gain =\
    getDetectorCIC(DET_CBE_Data, CIC_adjust, detEMgain, missionFraction, frameTime)
    
    CRhitsPerFrame,detPixAcross,detPixSize,CRrate,CRtailLen\
    = getDetectorCosmicRays(perfLevel,DET_CBE_Data, detEMgain, frameTime )
    
    CTE_clocking_efficiency, CTE_traps, signalPerPixPerFrame\
    = getCTE(DET_CBE_Data, rate_photoConverted,frameTime,missionFraction)
    
    hotPixFrac, hotPix = gethotPixels(DET_CBE_Data, missionFraction)
    
    ENF = getENF(isPhotonCounting)
    detCamRead = DET_CBE_Data.df.at[0,'ReadNoise_e']
    
    readNoise, readNoise_leakage, readNoise_leakage_in_current_units,\
         PCeffloss, readNoise_w_gain\
       = getReadNoiseandPCeffloss(detCamRead, detPCthreshold, isPhotonCounting, frameTime, detEMgain)
    det_CTE = CTE_clocking_efficiency * CTE_traps
    
    signal_region_electron_rate, det_PC_threshold_efficiency,\
        det_PC_coincid_effic, det_hotPix, det_cosmicRays, dQE\
    = getdetdQE( det_CTE, PCeffloss, hotPix,\
                signalPerPixPerFrame, detPixAcross, CRtailLen, CRhitsPerFrame,\
                    QE_adjust, det_QE)
       
    
    planetRate_proc, speckleRate_proc, zodiRate_proc,\
                ezo_bkgRate, lzo_bkgRate, residSpecRate =\
    getNoiseVarianceRatesI( f_SR, starFlux, planetFlux, colArea, rawContrast, \
              thpt_t_pnt, thpt_t_speckle, thpt_tau_pk, CGintmpix, k_pp_CBE,\
                     dQE, rate_exoZodi_incPht, rate_loZodi_incPht,\
                         selDeltaC)
    
    intTimeDutyFactor, allocTinteg, usableTinteg = intTime(perfLevel,scenarioData)
    Kappa = SNR_for_NEFR/(f_SR*starFlux*colArea*thpt_t_pnt*dQE\
                          *usableTinteg)/uc.ppb
    
    strayLight, strayLight_FRN, strayLight_ph_pix_s,\
        strayLight_ph_pix_h, strayLight_e_pix_h =\
            getStrayLightFRN(scenario, perfLevel, STRAY_FRN_Data, CG_Data, IWA, OWA,\
                     opMode, DET_CBE_Data, ENF, detPixSize, mpix, dQE, Kappa,\
                          usableTinteg, compStarSpecType, inBandFlux0_sum,\
                              f_SR, CGintmpix, colArea, thpt_t_speckle)

    totNoiseVarRate, readNoiseRate, luminesRate, noiseVarRate_perSNRregion,\
        CIC_RNLK_noiseRate, darkNoiseRate, zodi_shot, speckle_shot, planet_shot=\
            getNoiseVarianceRatesII(planetRate_proc, speckleRate_proc, zodiRate_proc,\
            ezo_bkgRate, lzo_bkgRate, residSpecRate, rate_planet_imgArea, ENF,\
                DarkCur_epoch_per_s, readNoise_leakage_in_current_units,\
                mpix, det_CIC_in_DC_units, readNoise, frameTime, isPhotonCounting,\
                    k_sp, k_lzo, k_ezo, k_det, dQE, strayLight_ph_pix_s)
                

    return planetRate_proc,usableTinteg, totNoiseVarRate, residSpecRate,\
        SNRtarget, ENF, k_sp, k_det, k_lzo, k_ezo, ezo_bkgRate, lzo_bkgRate,\
            speckleRate_proc, zodiRate_proc, Kappa, f_SR, starFlux, colArea, \
            darkNoiseRate,CIC_RNLK_noiseRate, readNoiseRate, thpt_t_pnt, dQE,\
                planetWA, lamD, IWA, OWA, CG_Data, thpt_t_refl,\
                    inBandZeroMagFlux,omegaPSF,allocTinteg, k_pp_CBE

def DRM_tableofNsigma(N_sigmaSens, infiniteTimeSens, planetWA, lamD):
    """Converts planet FRN sensitity values into ppb for the table of 
    sensitivities for """
    # Plot of Flux Ratio Detection Sensitivity for varying Planet separation-------------------
    
    nsigmasens_finite = N_sigmaSens/uc.ppb
    nsigmasens_inf = infiniteTimeSens/uc.ppb
    separation = planetWA*lamD/uc.mas
    
    return separation, nsigmasens_finite, nsigmasens_inf

def getQE(scenarioData, QE_Data):
    
    lam = float(scenarioData.at['CenterLambda_nm','Latest'])
    indQE = QE_Data.df.loc[QE_Data.df['lambda_nm']<=(lam),'lambda_nm'].idxmax()
    det_QE  = QE_Data.df.at[indQE,'QE_at_neg100degC']
    
    return det_QE

def getNoiseVarianceRatesII(planetRate_proc, speckleRate_proc, zodiRate_proc,\
            ezo_bkgRate, lzo_bkgRate, residSpecRate, rate_planet_imgArea, ENF,\
                DarkCur_epoch_per_s, readNoise_leakage_in_current_units,\
                mpix, det_CIC_in_DC_units, readNoise, frameTime, isPhotonCounting,\
                    k_sp, k_lzo, k_ezo, k_det, dQE, strayLight_ph_pix_s):    
    
    planet_shot = planetRate_proc * ENF**2
    speckle_shot = speckleRate_proc * ENF**2
    zodi_shot = zodiRate_proc * ENF**2
    darkNoiseRate = mpix * DarkCur_epoch_per_s * ENF**2


    if isPhotonCounting:
        CIC_RNLK_noiseRate = ENF**2*mpix*(det_CIC_in_DC_units +\
                                          readNoise_leakage_in_current_units)
    else:
        CIC_RNLK_noiseRate = ENF**2*mpix*(det_CIC_in_DC_units)
        
        
    noiseVarRate_perSNRregion = strayLight_ph_pix_s * mpix * dQE    
    luminesRate = noiseVarRate_perSNRregion * ENF**2
    
    readNoiseRate = (mpix/frameTime) * readNoise**2   
    
    totNoiseVarRate = ENF**2 * (rate_planet_imgArea +\
                              k_sp * speckleRate_proc +\
                              k_lzo * lzo_bkgRate +\
                              k_ezo * ezo_bkgRate ) +\
                        k_det * (darkNoiseRate + CIC_RNLK_noiseRate + luminesRate) +\
                        k_det * readNoiseRate
    
    return totNoiseVarRate, readNoiseRate, luminesRate, noiseVarRate_perSNRregion,\
        CIC_RNLK_noiseRate, darkNoiseRate, zodi_shot, speckle_shot, planet_shot

def throughput_pars(perfLevel, THPT_Data, scenarioData, CG_occulter_transmission,\
                    CGcoreThruput, PSFpeakI, CGtauPol, omegaPSF, DPM, lamD):
    
    thpt_t_obsc = THPT_Data.df.at[0,'Pupil_Transmission']
    
    colArea =((math.pi/4)*DPM**2)*thpt_t_obsc
    tau_core_REQ = scenarioData.at['REQ_t_core','Latest']
    
    try:
        if perfLevel=="REQ":
            thpt_OTA_TCA = THPT_Data.df.at[0,"REQ_OTAplusTCA"]
            thpt_CGI = THPT_Data.df.at[0,"REQ_CGI"]
            thpt_t_core = scenarioData.at['REQ_t_core','Latest']
            thpt_tau_pk = tau_core_REQ/CGcoreThruput*PSFpeakI
        elif perfLevel=="CBE":
            thpt_OTA_TCA = THPT_Data.df.at[0,"CBE_OTAplusTCA"]
            thpt_CGI = THPT_Data.df.at[0,"CBE_CGI"]
            thpt_t_core = CGcoreThruput
            thpt_tau_pk = PSFpeakI
    except:
        print("Specify perfLevel as REQ or CBE")
    
    thpt_t_occ = CG_occulter_transmission/CGtauPol
    
    try:
        thpt_t_PSF =  thpt_t_core/thpt_t_occ
    except:
        thpt_t_PSF = 0
    
    thpt_t_refl = thpt_OTA_TCA * thpt_CGI;
    
    thpt_t_pol  = CGtauPol  # Same for REQ & CBE
    
    thpt_t_unif = thpt_t_occ * thpt_t_refl * thpt_t_pol
    
    #thpt_t_unif * thpt_t_PSF --t point source used for planet    
    thpt_t_pnt = thpt_t_pol * thpt_t_refl * thpt_t_core
   
    thpt_t_speckle = thpt_t_refl * thpt_t_pol  
       
    # omega_PSF is in as2
    thpt_coreArea_LamD = omegaPSF/(lamD*0.001/uc.mas)**2 
    
    thpt_bkgLimitedCore = thpt_t_core/ math.sqrt( thpt_coreArea_LamD )
    
    
    return  thpt_t_obsc, colArea, thpt_OTA_TCA, thpt_CGI, thpt_bkgLimitedCore,\
            thpt_coreArea_LamD, thpt_t_PSF, thpt_t_core,\
                 tau_core_REQ, thpt_t_occ, thpt_t_pnt, thpt_t_pol,\
                    thpt_t_refl, thpt_t_speckle, thpt_t_unif, thpt_tau_pk
    


