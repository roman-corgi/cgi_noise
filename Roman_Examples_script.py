#!/usr/bin/env python
# coding: utf-8

# # Example Notebook

# This page allows you to do three study cases
# 
# A) Time to SNR or SNR in allowed time
# 
# B) Planet Sensitivity
# 
# C) Dust Sensitivity
# 
# User may define or select values specifically for each one; open the table of contents and look for "USER"
# Specific variables, such as integration time, can optionally be set independently for each case

# In[1]:


import os
import sys
current_dir = os.getcwd()
sys.path.append(current_dir)
from loadXLcol import loadXLcol
import pandas as pd
import math
import func_library as fl
from matplotlib import pyplot as plt
import matplotlib.patches as patches
import pdb
from pathlib import Path
from IPython.display import Image, display
from datetime import datetime
import unitsConstants as uc
import numpy as np


# In[2]:


filenamedir = current_dir
print(filenamedir)

current_datetime = datetime.now()


# ## Calculation Parameters Specification

# In[3]:


# Primary Mirror Diameter for Roman
DPM = 2.363 * uc.meter


# Scenario specification is done via an excel "scenario file"

# In[4]:


scenario_filename = 'SCEN_IMG_NFOV_B2_SPC.xlsx' #'SCEN_DRM_IMG_WF_B4.xlsx' # 'SCEN_DRM_SP_B2_TVAC.xlsx' # 'SCEN_DRM_IMG_WF_B3_TVAC.xlsx' # 'SCEN_DRM_IMG_NF_B1_HLC.xlsx' # 

print('Select a scenario from the list of available options.')
SCEN_Dict = fl.open_folder('EBcsvData','Scenarios')
for key in SCEN_Dict.keys():
    if key.find('.xlsx') != -1:
        print(key)


# ### Scenario Parameters


# ### Scenario Parameters

# In[5]:


#-------------------------------------------------------
# SET THE SCENARIO PARAMETERS
#-------------------------------------------------------
# Set Study Parameters
SNRtarget  = 5.0 # desired signal to noise ratio
maxTinteg  = 100 # integration time in hours
pp_factor  = 2   # k_pp, post-processing factor
monthsAtL2 = 21  # months at L2

perfLevel = 'CBE'

# Contrast Stability options (choose 'MCBE_'):
CSprefix = 'MCBE_' # 'ICBE_'(instrument CBE), 'MCBE_'(mission CBE)
CBEtype = "mission", # "instrument", "mission"

# Detector Options (boolean):
isPhotonCounting = True 

# # Set Study Parameter values as authoritative
# keepConsistent = 1 # 0 to redefine Study Parameters for each study case


# In[6]:


scenFolder = fl.open_folder("EBcsvData","Scenarios")

# Specify here how many rows of the Scenario Excel to read:
scenarioDF = loadXLcol(scenFolder[scenario_filename],30).df

scenario = scenarioDF.at['Scenario','Latest']
scenarioDF


# In[7]:


lam = scenarioDF.at['CenterLambda_nm','Latest'] * uc.nm
print(lam)


# In[8]:


lamD = lam / DPM

print(f'Lambda/D = {lamD:5.3e} rad')
print(f'Lambda/D = {lamD/uc.mas:5.3f} mas')


# ### Host Star & Planet Specs

# In[9]:


# User set host star and planet specifications

hostStar_vmag = 5.05
hostStar_dist = 13.8 # parsec
hostStar_type = 'g0v'

planet_phaseAngle = 65   # degrees
planet_sma_AU     = 4    # AU
planet_radius     = 1 # R_Jupiter
planet_GeomAlbedo = 0.3

sep_mas    = fl.Target.phaseAng_to_sep(planet_sma_AU, hostStar_dist, planet_phaseAngle)

print(f'Separation = {sep_mas:5.0f} mas')

exoZodi = 1

detPCthreshold = 5.0


# In[10]:


# Assign core throughput values from file to local variables
thpt_t_pnt = scenarioDF.at['t_core_CBEmeas','Latest']
print(thpt_t_pnt)
RefStarSpecType = scenarioDF.at['RefStar_SpectralType','Latest']
print(RefStarSpecType)


# In[11]:


#-------------------------------------------------
# Derive Delta Mag from Desired Flux Ratio
#------------------------------------------------
# Set Flux Ratio and derive the deltaMag
fluxRatio0 = 5e-9
planet_deltaMag = fl.Target.fluxRatio_to_deltaMag(fluxRatio0)

print(f'Planet delta-mag = {planet_deltaMag:5.8e}')

#-------------------------------------------------
# Set Delta Mag and derive the Flux ratio
#-------------------------------------------------

# Uncomment to set the deltaMag and derive the fluxRatio
#planet_deltaMag = 20.752574989159953
#planet_fluxRatio = 10**(-0.4*planet_deltaMag)

planet_fluxRatio = fl.Target.deltaMag_to_fluxRatio(planet_deltaMag)

print(f'Flux Ratio       = {planet_fluxRatio:5.8e}')


# In[12]:


# set up target specs
target = fl.Target( v_mag = hostStar_vmag,\
                   phaseAng_deg = planet_phaseAngle,\
                   geomAlb_ag = planet_GeomAlbedo, \
                   radius_Rjup = planet_radius, \
                   dist_pc = hostStar_dist,\
                   exoZodi = exoZodi,\
                   specType = hostStar_type,\
                   sma_AU = planet_sma_AU
                   )


sma_AU = target.sep_and_dist_to_SMA(sep_mas, hostStar_dist, planet_phaseAngle)
target.sma_AU = sma_AU
target.albedo = target.fluxRatio_SMA_rad_to_albedo( planet_fluxRatio, sma_AU, planet_radius)

print('Target Specs:')
print(f'V:                 {target.v_mag}')
print(f'DIST:              {target.dist_pc}')
print(f'Geometric Albedo:  {target.geomAlb_ag}')
print(f'ExoZodi xSolar:    {target.exoZodi}')
print(f'SMA au:            {target.sma_AU}')
print(f'Spec type:         {target.specType}')
print(f'Phase Angle (deg): {target.phaseAng_deg}')
print(f'Albedo:            {target.albedo}')
print(f'Radius Rjup        {target.radius_Rjup}')
print(f'Separation_mas     {sep_mas}')


# ## Planet Sensitivity

# In[13]:


# get list of specific csv files for selected mode
filenameList = fl.getScenFileNames_DRM(scenarioDF)

# get the data tables of the mode-specific csv files
[CG_Data, QE_Data, DET_CBE_Data, STRAY_FRN_Data, THPT_Data, CAL_Data, CS_Data] \
    = fl.loadCSVs(filenameList)

# calculate planet working angle from given target specifications
IWA, OWA = fl.workingAnglePars(CG_Data, CS_Data)

planetWA = fl.Target.separation_to_planetWA(sep_mas, lamD, IWA, OWA)


# Constrast Stability for Mission CBE
selDeltaC, selContrast, SystematicCont, initStatRawContrast,\
    rawContrast, IntContStab, ExtContStab\
        = fl.contrastStabilityPars( CSprefix, planetWA, CS_Data)

det_QE = fl.getQE(scenarioDF, QE_Data)
detDarkBOM, detDarkEOM, DarkCur_epoch_per_s, DarkCur_epoch_per_hr,\
    missionFraction, detEOL_mos\
    = fl.getDarkCurrent(DET_CBE_Data, monthsAtL2)

CGtauPol, indWA, CGcoreThruput, PSFpeakI, omegaPSF, CGintSamp,\
    CGradius_arcsec, CGdesignWL, CGintmpix, CG_PSFarea_sqlamD, CGintensity,\
        CG_occulter_transmission, CGcontrast\
            = fl.coronagraph_pars(CG_Data, planetWA, IWA, DPM, lamD)


# In[14]:


thpt_t_obsc = THPT_Data.df.at[0,'Pupil_Transmission']

colArea =((math.pi/4)*DPM**2)*thpt_t_obsc

thpt_OTA_TCA = THPT_Data.df.at[0,"CBE_OTAplusTCA"]
thpt_CGI = THPT_Data.df.at[0,"CBE_CGI"]
thpt_t_core = CGcoreThruput
thpt_tau_pk = PSFpeakI

thpt_t_occ = CG_occulter_transmission/CGtauPol

try:
    thpt_t_PSF =  thpt_t_core/thpt_t_occ
except:
    thpt_t_PSF = 0

thpt_t_refl = thpt_OTA_TCA * thpt_CGI;
thpt_t_pol  = CGtauPol
thpt_t_unif = thpt_t_occ * thpt_t_refl * thpt_t_pol
thpt_t_speckle = thpt_t_refl * thpt_t_pol  

# omega_PSF is in as2
thpt_coreArea_LamD = omegaPSF/(lamD*0.001/uc.mas)**2 
thpt_bkgLimitedCore = thpt_t_core/ math.sqrt( thpt_coreArea_LamD )


# In[15]:


# Assign values from scenario file to local variables
opMode = scenarioDF.at['OperatingMode','Latest']
bandWidth = scenarioDF.at['BW','Latest']
RefStarVmag_CBE = scenarioDF.at['RefStar_V_mag_CBE','Latest']
RefStarDist = scenarioDF.at['RefStar_Distance_pc','Latest']
RefStarExoZodi = scenarioDF.at['RefStar_ExoZodi_Xsolar','Latest']
timeOnRef = scenarioDF.at['TimeonRefStar_tRef_per_tTar','Latest']
magLocalZodi = scenarioDF.at['LocZodi_magas2','Latest']
magExoZodi_1AU = scenarioDF.at['ExoZodi_magas2','Latest']
k_comp = 0
FWC_gr = 90000


# In[ ]:





# In[19]:


f_SR, CritLam, detPixSize_m, mpix, pixPlateSc\
    = fl.getFocalPlaneAttributes(opMode, scenarioDF, DET_CBE_Data, lam, bandWidth,\
                        DPM, CGdesignWL, omegaPSF)

inBandFlux0_sum, inBandZeroMagFlux, starFlux\
    = fl.getSpectra(target, lam, bandWidth)

print(starFlux)


# In[17]:


fluxRatio, planetFlux = fl.getFluxRatio(target, starFlux)

RefStarVmag, RefStarAbsMag, RefStarinBandZeroMagFlux, RefStarDeltaMag,\
     RefStarFlux, BrightnessRatio, betaRDI, k_sp, k_det, k_lzo, k_ezo,\
         v_sp, v_det, v_lzo, v_ezo \
             = fl.getrefStarRDI(target, inBandFlux0_sum, starFlux,\
              RefStarSpecType, RefStarVmag_CBE,\
                  RefStarDist, RefStarExoZodi, timeOnRef)

ZodiFlux, exoZodiAngFlux, loZodiAngFlux, loZodiFlux, exoZodiFlux, absMag \
    = fl.getZodi(magLocalZodi, magExoZodi_1AU, target, inBandZeroMagFlux, omegaPSF)

rate_planet_imgArea, rate_Zodi_imgArea, rate_exoZodi_incPht,\
    rate_loZodi_incPht, rate_speckleBkg, rate_photoConverted,\
        rate_totalwithoutplanet \
    = fl.getNoiseRates(f_SR, starFlux, fluxRatio, colArea, thpt_t_pnt, thpt_tau_pk,\
              thpt_t_speckle, det_QE, exoZodiFlux, loZodiFlux, thpt_t_unif, k_comp,\
              rawContrast,CGintmpix, mpix, DarkCur_epoch_per_s)

frameTime, frameTime_ANLG,maxANLGt_fr,maxPCt_fr, detEMgain\
    = fl.getFrameExposureTime(DET_CBE_Data, FWC_gr, rate_totalwithoutplanet,\
                     rate_photoConverted, isPhotonCounting)


# In[18]:


det_CIC_in_DC_units, det_CIC_at_epoch, det_CIC_gain\
    = fl.getDetectorCIC(DET_CBE_Data, CIC_adjust, detEMgain, missionFraction, frameTime)

CRhitsPerFrame,detPixAcross,detPixSize,CRrate,CRtailLen\
    = fl.getDetectorCosmicRays(perfLevel, DET_CBE_Data, detEMgain, frameTime )

CTE_clocking_efficiency, CTE_traps, signalPerPixPerFrame\
    = fl.getCTE(DET_CBE_Data, rate_photoConverted,frameTime,missionFraction)

hotPixFrac, hotPix = fl.gethotPixels(DET_CBE_Data, missionFraction)

ENF = fl.getENF(isPhotonCounting)
detCamRead = DET_CBE_Data.df.at[0,'ReadNoise_e']

readNoise, readNoise_leakage, readNoise_leakage_in_current_units,\
     PCeffloss, readNoise_w_gain\
   = fl.getReadNoiseandPCeffloss(detCamRead, detPCthreshold, isPhotonCounting, frameTime, detEMgain)
det_CTE = CTE_clocking_efficiency * CTE_traps

signal_region_electron_rate, det_PC_threshold_efficiency,\
    det_PC_coincid_effic, det_hotPix, det_cosmicRays, dQE\
    = fl.getdetdQE( det_CTE, PCeffloss, hotPix,\
            signalPerPixPerFrame, detPixAcross, CRtailLen, CRhitsPerFrame,\
                QE_adjust, det_QE)

planetRate_proc, speckleRate_proc, zodiRate_proc,\
            ezo_bkgRate, lzo_bkgRate, residSpecRate\
    = fl.getNoiseVarianceRatesI( f_SR, starFlux, planetFlux, colArea, rawContrast, \
          thpt_t_pnt, thpt_t_speckle, thpt_tau_pk, CGintmpix, k_pp_CBE,\
                 dQE, rate_exoZodi_incPht, rate_loZodi_incPht,\
                     selDeltaC)


# In[ ]:


intTimeDutyFactor, allocTinteg, usableTinteg = fl.intTime(perfLevel,scenarioDF)   

Kappa = SNR_for_NEFR/(f_SR*starFlux*colArea*thpt_t_pnt*dQE\
                      *usableTinteg)/uc.ppb

strayLight, strayLight_FRN, strayLight_ph_pix_s,\
    strayLight_ph_pix_h, strayLight_e_pix_h =\
        fl.getStrayLightFRN(scenario, perfLevel, STRAY_FRN_Data, CG_Data, IWA, OWA,\
                 opMode, DET_CBE_Data, ENF, detPixSize, mpix, dQE, Kappa,\
                      usableTinteg, hostStar_type, inBandFlux0_sum,\
                          f_SR, CGintmpix, colArea, thpt_t_speckle)

totNoiseVarRate, readNoiseRate, luminesRate, noiseVarRate_perSNRregion,\
    CIC_RNLK_noiseRate, darkNoiseRate, zodi_shot, speckle_shot, planet_shot=\
        fl.getNoiseVarianceRatesII(planetRate_proc, speckleRate_proc, zodiRate_proc,\
        ezo_bkgRate, lzo_bkgRate, residSpecRate, rate_planet_imgArea, ENF,\
            DarkCur_epoch_per_s, readNoise_leakage_in_current_units,\
            mpix, det_CIC_in_DC_units, readNoise, frameTime, isPhotonCounting,\
                k_sp, k_lzo, k_ezo, k_det, dQE, strayLight_ph_pix_s)


# ## Calculations for Specified Planet
# Substitute for usableTinteg:
usableTmini = 0.71 * 60 * 60


# In[ ]:


instaPlanetDF1 = fl.DRMinstaplanet(planetRate_proc, usableTinteg, totNoiseVarRate,\
                             residSpecRate, SNRtarget)
timetoSNR = instaPlanetDF1.loc[0,'Value']
print(f'Time to {SNRtarget} SNR = {timetoSNR:5.2f} hrs')
round(instaPlanetDF1,2)


# ### Reverse Calculation -- SNR in allowed time


# ### Reverse Calculation -- SNR in allowed time

# In[ ]:


# 2 Knobs to Turn:
allocatedTinteg = 100 # allocated integration time in hours
usableFrac      = 1 # usable fraction

usableTinteg = allocatedTinteg * usableFrac *3600
print(f'Usable Integration Time (usable) = {usableTinteg:5.2f} seconds')
instaPlanetDF2 = fl.DRMinstaplanet(planetRate_proc, usableTinteg,\
                                                totNoiseVarRate, residSpecRate, SNRtarget)
SNRinAllowedTime = instaPlanetDF2.loc[5,'Value']
print(f'SNR in allowed time = {SNRinAllowedTime:5.2f}')
round(instaPlanetDF2,2)


# ## Planet Sensitivity

# ### N-Sigma Sensitivity


# ## Planet Sensitivity

# In[ ]:


# Constrast Stability for Mission CBE
selDeltaC, selContrast, SystematicCont, initStatRawContrast,\
    rawContrast, IntContStab, ExtContStab\
        = fl.contrastStabilityPars( CSprefix, planetWA, CS_Data)


# In[ ]:


thpt_t_obsc = THPT_Data.df.at[0,'Pupil_Transmission']

colArea =((math.pi/4)*DPM**2)*thpt_t_obsc

thpt_OTA_TCA = THPT_Data.df.at[0,"CBE_OTAplusTCA"]
thpt_CGI = THPT_Data.df.at[0,"CBE_CGI"]
thpt_t_core = CGcoreThruput
thpt_tau_pk = PSFpeakI

try:
    thpt_t_PSF =  thpt_t_core/thpt_t_occ
except:
    thpt_t_PSF = 0

thpt_t_refl = thpt_OTA_TCA * thpt_CGI;

thpt_t_unif = thpt_t_occ * thpt_t_refl * thpt_t_pol

thpt_t_speckle = thpt_t_refl * thpt_t_pol

thpt_bkgLimitedCore = thpt_t_core/ math.sqrt( thpt_coreArea_LamD )


# In[ ]:


# Back to DRMgetSNRvals
f_SR, CritLam, detPixSize_m, mpix, pixPlateSc\
    = fl.getFocalPlaneAttributes(opMode, scenarioDF, DET_CBE_Data, lam, bandWidth,\
                        DPM, CGdesignWL, omegaPSF)

inBandFlux0_sum, inBandZeroMagFlux, starFlux\
    = fl.getSpectra(target, lam, bandWidth)

fluxRatio, planetFlux = fl.getFluxRatio(target, starFlux)


# ### N-Sigma Sensitivity

# In[ ]:


# fl.getrefStarRDI
RefStarAbsMag = RefStarVmag_CBE - 5*math.log10(RefStarDist/10)
RefStarinBandZeroMagFlux = inBandFlux0_sum.at[RefStarSpecType]
RefStarDeltaMag = target.v_mag - RefStarVmag_CBE
RefStarFlux = RefStarinBandZeroMagFlux*(10**((-0.4)*RefStarVmag_CBE))
BrightnessRatio = RefStarFlux/starFlux

# percentage of time on reference star
betaRDI = 1 / (BrightnessRatio*timeOnRef)

k_sp = 1 + betaRDI
k_det = 1 + betaRDI**2 * timeOnRef
k_lzo = k_det
k_ezo = k_sp

v_sp = math.sqrt(betaRDI)
v_det = betaRDI * math.sqrt(timeOnRef)
v_lzo = betaRDI * math.sqrt(timeOnRef)
v_ezo = math.sqrt(betaRDI)


# In[ ]:


absMag = target.v_mag - 5 * math.log10(target.dist_pc / 10)

# zodi solid angular flux in units of ph/s/m**2/arcsec**2
loZodiAngFlux =  inBandZeroMagFlux *10**((-0.4)* magLocalZodi) 
exoZodiAngFlux = target.exoZodi * inBandZeroMagFlux * \
              10**(-0.4*(absMag-uc.sunAbsMag+magExoZodi_1AU)) / (target.sma_AU**2)

ZodiFlux    = (loZodiAngFlux + exoZodiAngFlux) * omegaPSF
loZodiFlux  = loZodiAngFlux * omegaPSF
exoZodiFlux = exoZodiAngFlux * omegaPSF


# #### Non-Stellar Light Rates


# In[ ]:


# getNoiseRates
rate_planet_imgArea = f_SR * starFlux * fluxRatio * colArea * thpt_t_pnt * det_QE
rate_exoZodi_incPht = f_SR * exoZodiFlux * colArea * thpt_t_unif
rate_loZodi_incPht =  f_SR * loZodiFlux * colArea * thpt_t_unif 
rate_Zodi_imgArea = (rate_exoZodi_incPht + rate_loZodi_incPht)*det_QE

# speckleRate_imgArea
rate_speckleBkg = f_SR * starFlux * rawContrast * uc.ppb * thpt_tau_pk * \
    CGintmpix * colArea * thpt_t_speckle  * det_QE 

# total rate pixel with planet
rate_photoConverted = DarkCur_epoch_per_s + \
    (rate_planet_imgArea + rate_Zodi_imgArea + rate_speckleBkg )/mpix
# Why divide other values by mpix and not Dark Current?
rate_totalwithoutplanet =  DarkCur_epoch_per_s + (rate_Zodi_imgArea +\
                                                  rate_speckleBkg )/mpix


# In[ ]:


noise_rates =\
pd.DataFrame({f'{scenario}':['Planet Signal rate',\
                             'Exo Zodi rate',\
                             'Local Zodi rate',\
                             'Combined Zodi rate',\
                             'Speckle Background rate',\
                             'Average PSF rate w Planet',\
                             'Average PSF rate w/o Planet'],\
                    'Per core':[rate_planet_imgArea,\
                                rate_exoZodi_incPht,\
                                rate_loZodi_incPht,\
                                rate_Zodi_imgArea,\
                                rate_speckleBkg,\
                                rate_photoConverted*mpix,\
                                rate_totalwithoutplanet*mpix],\
                    'Per pixel':[rate_planet_imgArea/mpix,\
                                 rate_exoZodi_incPht/mpix,\
                                 rate_loZodi_incPht/mpix,\
                                 rate_Zodi_imgArea/mpix,\
                                 rate_speckleBkg/mpix,\
                                 rate_photoConverted,\
                                 rate_totalwithoutplanet],\
              'Units':['e-/s','ph/s','ph/s','e-/s','e-/s','e-/s','e-/s']})
noise_rates = noise_rates.style.format({"Value": "{:.2e}"})
noise_rates


# In[ ]:


# getFrameExposureTime
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


# In[ ]:


# getDetectorCIC
detCICdegradationEOM = DET_CBE_Data.df.at[0,'CICdegradationEOM']
detCIC1 = DET_CBE_Data.df.at[0,'CICatGain1BOM_e_per_pix_per_fr']
detCIC2 = DET_CBE_Data.df.at[0,'CICatGain2BOM_e_per_pix_per_fr']
detCICgain1 = DET_CBE_Data.df.at[0,'Gain1BOM']
detCICgain2 = DET_CBE_Data.df.at[0,'Gain2BOM']

det_CIC_gain = ((detCIC2 - detCIC1)/(detCICgain2 - detCICgain1))*detEMgain\
    +(detCIC1-((detCIC2 - detCIC1)/(detCICgain2 - detCICgain1))*detCICgain1)\
        *(1 + missionFraction*(detCICdegradationEOM-1))
det_CIC_at_epoch = CIC_adjust*det_CIC_gain
det_CIC_in_DC_units = det_CIC_at_epoch/frameTime


# In[ ]:


# getDetectorCosmicRays
# Cosmic Ray Tail Length
detCRtailGain1 = DET_CBE_Data.df.at[0,'CRtailGain1']
detCRtailGain2 = DET_CBE_Data.df.at[0,'CRtailGain2']
detCRtailLen1 = DET_CBE_Data.df.at[0,'CRtailLen1_pix']
detCRtailLen2 = DET_CBE_Data.df.at[0,'CRtailLen2_pix']
CRtailLen = ((detCRtailLen2-detCRtailLen1)/(detCRtailGain2-detCRtailGain1))*detEMgain+(detCRtailLen2- ((detCRtailLen2-detCRtailLen1)/(detCRtailGain2-detCRtailGain1))*detCRtailGain2)
CRrate = 5*(100**2) # rays/s/m^2  
detPixSize = DET_CBE_Data.df.at[0,'PixelSize_m']
detPixAcross = DET_CBE_Data.df.at[0,'PixelsAcross_pix']
CRhitsPerFrame = CRrate * frameTime * (detPixSize * detPixAcross)**2 # Cosmic Rays per frame


# In[ ]:


# getCTE
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


# In[ ]:


# getHotPixels
hotPixFrac = DET_CBE_Data.df.at[0,'HotPixFrac']
hotPix = hotPixFrac*missionFraction


# In[ ]:


# getENF
if isPhotonCounting:
    ENF = 1
else:
    ENF = math.sqrt(2)


# In[ ]:


# getReadNoiseandPCeffloss
detCamRead = DET_CBE_Data.df.at[0,'ReadNoise_e']
readNoise_w_gain = detCamRead/detEMgain # read noise with gain if analog
if isPhotonCounting:
    readNoise = 0 # Minimal read noise with photon counting
    PCeffloss = 1 - math.exp( -detPCthreshold*detCamRead/detEMgain)
else:
    readNoise = readNoise_w_gain # read noise with gain if analog
    PCeffloss = 0
readNoise_leakage = 0.5*math.erfc((detPCthreshold/math.sqrt(2)))
readNoise_leakage_in_current_units = readNoise_leakage/frameTime


# In[ ]:


# getdetCTE
det_CTE = CTE_clocking_efficiency * CTE_traps


# In[ ]:


# getdetdQE
det_PC_threshold_efficiency = 1 - PCeffloss
signal_region_electron_rate = signalPerPixPerFrame * det_CTE

# Photon-counting coincidence efficiency 
det_PC_coincid_effic = (1 - math.exp(-signal_region_electron_rate))\
    / (signal_region_electron_rate)

# Efficiency after subtracting fraction lost to hot pixels
det_hotPix = 1 - hotPix

det_cosmicRays = 1 - CRhitsPerFrame * CRtailLen/(detPixAcross**2)
dQE = det_QE * det_CTE * det_PC_threshold_efficiency * \
        det_PC_coincid_effic * det_hotPix * det_cosmicRays * QE_adjust


# In[ ]:


# getNoiseVarianceRatesI
planetRate_proc = f_SR*planetFlux*colArea*thpt_t_pnt*dQE

speckleRate_proc = f_SR * starFlux * rawContrast * thpt_tau_pk\
     * CGintmpix * thpt_t_speckle * colArea * dQE * uc.ppb

ezo_bkgRate = rate_exoZodi_incPht * dQE
lzo_bkgRate = rate_loZodi_incPht * dQE
zodiRate_proc = ezo_bkgRate + lzo_bkgRate

residSpecRate =  f_SR * starFlux * (selDeltaC/k_pp_CBE) *\
     thpt_tau_pk * CGintmpix * thpt_t_speckle * colArea * dQE 

# Next - display & describe values from cell above; initially copied from display of first cell in 1.4.1.1
noiseVar_rates =\
pd.DataFrame({f'{scenario}':['Planet rate',\
                             'Exo Zodi rate',\
                             'Local Zodi rate',\
                             'Combined Zodi rate',\
                             'Speckle Background rate',\
                             'Average PSF rate w Planet',\
                             'Average PSF rate w/o Planet'],\
                    'Value':[planetRate_proc,\
                             rate_exoZodi_incPht,\
                             rate_loZodi_incPht,\
                             rate_Zodi_imgArea,\
                             rate_speckleBkg,\
                             rate_photoConverted,\
                             rate_totalwithoutplanet],\
              'Units':['e-/s','ph/s','ph/s','e-/s','e-/s','e-/s','e-/s']})
noiseVar_rates = noiseVar_rates.style.format({"Value": "{:.2e}"})
noiseVar_rates


# In[ ]:


# intTime
intTimeDutyFactor = scenarioDF.at['DutyFactor_CBE','Latest']
allocTinteg = scenarioDF.at['t_max_hrs','Latest']   
# Integration Time in seconds
usableTinteg = allocTinteg * intTimeDutyFactor * uc.hour


# In[ ]:


Kappa = SNR_for_NEFR/(f_SR*starFlux*colArea*thpt_t_pnt*dQE\
                      *usableTinteg)/uc.ppb


# In[ ]:


# getStrayLightFRN
try:
    strayLight = fl.getStrayLightfromfile(scenario, perfLevel, STRAY_FRN_Data)
except:
    photonRateFlux = fl.getStrayLight_luminescenceBackground(perfLevel,opMode,\
                                                          DET_CBE_Data,\
                                                              mpix, dQE, False)
    leakageRateSL, leakageRate, strayLightBackground\
        = fl.getStrayLight_companionStarBackground(hostStar_type,\
                                                          IWA, OWA, CG_Data,\
                                      detPixSize, inBandFlux0_sum, f_SR,\
                                          CGintmpix, usableTinteg, mpix,\
                                          colArea, thpt_t_speckle, dQE, False)
    strayLight,_,_,_,_,_,_ = fl.getStrayLightfromcalc(CG_Data, photonRateFlux,\
                                                     leakageRateSL, IWA, OWA)

# This is the flux ratio noise from stray light:
strayLight_FRN = Kappa * math.sqrt(ENF**2 * strayLight * 1000000 \
                               * detPixSize**2 * mpix * dQE * usableTinteg) 
strayLight_ph_pix_s = strayLight * 1000000 * detPixSize**2
strayLight_ph_pix_h = strayLight_ph_pix_s * 3600
strayLight_e_pix_h = strayLight_ph_pix_h * dQE


# In[ ]:


# getNoiseVarianceRatesII
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

totNoiseVarRate = ENF**2 * (rate_planet_imgArea + k_sp * speckleRate_proc +\
                             k_lzo * lzo_bkgRate + k_ezo * ezo_bkgRate ) +\
                k_det * (darkNoiseRate + CIC_RNLK_noiseRate + luminesRate) +\
                k_det * readNoiseRate


# In[ ]:


# Noise rates for planet sensitivity

nonPlanetVarRate, nonpl_random, tot_nonpl_noise, N_sigmaSens =\
fl.DRM_planetSens_Nsigma(residSpecRate, usableTinteg, ENF, k_sp, speckleRate_proc,\
                         k_lzo, lzo_bkgRate, k_ezo, ezo_bkgRate, k_det,\
                             darkNoiseRate, CIC_RNLK_noiseRate, readNoiseRate,\
                                SNRtarget, Kappa)

print(f'At planet WA = {planetWA} lam/D,',end='')
print(f' and star v mag = {target.v_mag}, N-sigma Sensitivity = {N_sigmaSens}')     


# ### Contributions to Planet Sensitivity and FRN


# ### Contributions to Planet Sensitivity and FRN

# In[ ]:


# Calculates the FRN sensitivity and contributors to the sensitivity
nonplanetRandom, SNRcheck, planet_implicit, speckle, zodi, dark, cicrnleakage, read\
= fl.DRM_planetSens_contribtoNsigmaSens(N_sigmaSens, tot_nonpl_noise, ENF,\
                                       k_sp, speckleRate_proc, zodiRate_proc,\
                                k_det, darkNoiseRate, CIC_RNLK_noiseRate,\
                                   readNoiseRate, usableTinteg, f_SR,\
                                     starFlux, colArea, thpt_t_pnt, dQE)


# ### N-sigma Sensitivity with Infinite Integration Time


# ### N-sigma Sensitivity with Infinite Integration Time

# In[ ]:


tot_nonpl_noise_inftime, infiniteTimeSens =\
fl.DRM_planetSens_NsigmaSensInfiniteTime(k_det, readNoiseRate, usableTinteg,\
                                          nonpl_random, SNRtarget,\
                                         residSpecRate, Kappa)

print(f'N-sigma Sensitivity with Infinite time = {infiniteTimeSens:11.10e}') # NEED TO CHECK


# In[ ]:


res_Speckle, randomNoiseRate, nonplrandom, totnonplnoise,\
        cstab, bde, totalFRN =\
fl.DRM_planetSens_calcFRN(ENF, k_sp, k_lzo, k_ezo, k_det, lzo_bkgRate, \
                           ezo_bkgRate, darkNoiseRate, CIC_RNLK_noiseRate,\
                               readNoiseRate, residSpecRate, speckleRate_proc,\
                                   usableTinteg, Kappa)


# #### Non-Planetary Noise


# In[ ]:


planetSensdf1 =\
pd.DataFrame({f'{scenario}':['Star Vmag',\
                             'Pl. Separation',\
                             'Pl. Wrk. Angle','SNR',\
                             'Integ Time', 'post proc gain',\
                             'res. Speckle','Planet implicit',\
                             'speckle','zodi','Dark',\
                             'CIC + RN Leakage', 'Read',\
                             'non-pl random','tot nonpl noise',\
                             'N-sigma Sens'],\
                             'Value':[target.v_mag,\
                                       sep_mas,\
                                       planetWA,\
                                       SNRtarget, allocTinteg,\
                                       k_pp_CBE,\
                                       res_Speckle,\
                                       planet_implicit, speckle,\
                                       zodi, dark,\
                                       cicrnleakage, read,\
                                       nonpl_random,\
                                       tot_nonpl_noise,\
                                     f'{N_sigmaSens:3.2e}'],\
              'Units':['mag','mas','lam/D','SNR','hrs','k_pp',\
                      'e-','e-','e-','e-','e-','e-','e-','e-',\
                       'e- ($\\sigma$_o)','minimum']})
planetSensdf1
#pd.set_option('display.float_format', lambda x: '%3.2f' % x)


# In[ ]:


planetSensdf2 =\
pd.DataFrame({f'{scenario}':['SNR check',\
                             'nonPlanetVarRate',\
                             'non-pl random above',\
                             'non-pl random',\
                             'res. Speckle',\
                             'tot nonpl noise',\
                             'Infinite time sensitivity'],\
                             'Value':[SNRcheck,\
                                       nonPlanetVarRate,\
                                       nonplanetRandom,\
                                       nonpl_random,\
                                       res_Speckle,\
                                       tot_nonpl_noise_inftime,\
                                      f'{infiniteTimeSens:3.2e}'],\
              'Units':['','e/SR/s','for checking','e-',\
                       'e-', 'e- ($\\sigma$_o)','']})
planetSensdf2 
#pd.set_option('display.float_format', lambda x: '%3.2f' % x)


# In[ ]:


planetSensdf3 =\
pd.DataFrame({f'{scenario}':['Star Vmag',\
                             'Pl. Separation',\
                             'Pl. Wrk. Angle','SNR',\
                             'Integ Time', 'resid. Spckle',\
                             'res. Speckle',\
                             'non-pl random',\
                             'tot nonpl noise',\
                             'C Stab.','BDE','Total'],\
                             'Value':[ target.v_mag,\
                                       sep_mas,\
                                       planetWA,\
                                       SNRtarget, allocTinteg,\
                                       k_pp_CBE,\
                                       res_Speckle,\
                                       nonpl_random,\
                                       tot_nonpl_noise,\
                                       cstab, bde, totalFRN],\
              'Units':['mag','mas','lam/D','SNR','hrs','k_pp',\
                      'e-','e-','e- ($\\sigma$_o)','ppb','ppb','ppb']})
planetSensdf3


# ### Plot of N-sigma Sensitivities for Finite and Infinite Time

# In[ ]:


# TRY TO CONSOLIDATE ALL DATA into one dataframe

planetSensdf4 =\
pd.DataFrame({f'{scenario}':['Star Vmag',\
                             'Pl. Separation',\
                             'Pl. Wrk. Angle','SNR',\
                             'Integ Time', 'post proc gain',\
                             'res. Speckle','Planet implicit',\
                             'speckle','zodi','Dark',\
                             'CIC + RN Leakage', 'Read',\
                             'non-pl random','tot nonpl noise',\
                             'N-sigma Sens','SNR check',\
                             'nonPlanetVarRate',\
                             'non-pl random above',\
                             'Infinite time sensitivity',\
                            'C Stab.','BDE','Total'],\
                             'Value':[target.v_mag,\
                                       sep_mas,\
                                       planetWA,\
                                       SNRtarget, allocTinteg,\
                                       k_pp_CBE,\
                                       res_Speckle,\
                                       planet_implicit, speckle,\
                                       zodi, dark,\
                                       cicrnleakage, read,\
                                       nonpl_random,\
                                       tot_nonpl_noise,\
                                       f'{N_sigmaSens:3.2e}',\
                                       SNRcheck,\
                                       nonPlanetVarRate,\
                                       nonplanetRandom,\
                                      f'{infiniteTimeSens:3.2e}',\
                                     cstab, bde, totalFRN],\
              'Units':['mag','mas','lam/D','SNR','hrs','k_pp',\
                      'e-','e-','e-','e-','e-','e-','e-','e-',\
                       'e- ($\\sigma$_o)','minimum','','e/SR/s',\
                       'for checking','e- ($\\sigma$_o)',\
                       'ppb','ppb','ppb']})
#pd.set_option('display.float_format', lambda x: '%3.2f' % x)
planetSensdf4 


# ### Plot of N-sigma Sensitivities for Finite and Infinite Time 
# User may select the start and end radial slices (inclusive) in 
# lambda/D which by default are the inner and outer working angles of the dark hole.  


# In[ ]:


numWAs = 10 # number of radial slices including the IWA and OWA
planetWAs_1 = np.linspace(IWA,OWA, numWAs)

nsigmasens_finites_1 = []
nsigmasens_infs = []
separations_1 = []
for planetWA_1 in planetWAs_1:

    nsigmasens_finite, nsigmasens_inf, separation,\
    planet_implicit, speckle, zodi, cstab, bde, total =\
    fl.DRM_planetSens_vals(scenarioDF, perfLevel, CSprefix, isPhotonCounting, target, planetWA, lamD, detPCthreshold)

    nsigmasens_finites_1.append(nsigmasens_finite)
    nsigmasens_infs.append(nsigmasens_inf)
    separations_1.append(separation)


NsigmaSensDF = pd.DataFrame(list(zip(planetWAs_1, nsigmasens_finites_1, nsigmasens_infs, separations_1)),\
                   columns = ['Working Angle','N-sigma Sens','N-sigma Sens Inf Time','sep (mas)'])
round(NsigmaSensDF,4)


# In[ ]:


fig, ax = plt.subplots(figsize=(9, 6))

im1 = ax.plot(separations_1, nsigmasens_finites_1,'--bo')
im2 = ax.plot(separations_1, nsigmasens_infs,'--ro')
ax.set_title(f'{SNRtarget}-sigma {scenario} [{lam}nm, {usableTinteg/3600} hrs, PP = {k_pp_CBE}x]')
ax.set_ylabel('Flux Ratio Detection Sensitivity,   ppb')
ax.set_xlabel('Planet Separation,    milli-arseconds')
plt.legend([f'V = {hostStar_vmag}','t = inf'], loc = 'upper right')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)  # Customize grid style

# the time quoted is the usable target integration time (regardless of duty factor) 
datetime_string = current_datetime.strftime("%Y-%m-%d %H:%M:%S")
righttext = datetime_string + ' -- ' + scenario_filename 
# Add custom text
plt.text(
    1.02,           # x position (slightly right of axes)
    0.5,            # y position (middle of figure)
    righttext,  # your custom text (\n for line break if needed)
    transform=ax.transAxes,   # use axes coordinates
    rotation=90,              # vertical text
    verticalalignment='center',
    horizontalalignment='center',
    fontsize=9,
    color = 'gray'
)
plt.show()
# Adjust layout to prevent clipping
#plt.tight_layout()


# Caveat: Planet sensitivity calculation is approximate and same albedo is used to get same frame time for all working angles

# ### Flux Ratio Noise
# Calculating flux ratio noise in case question comes up whether the annular zone we have selected is adequately representative of behavior near IWA


# ### Flux Ratio Noise

# In[ ]:


planetWAs_2 = [3.1,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0]
separations_2 = []
cstabs = []
bdes = []
totals = []
for m in range(11):

    planetWA_2 = planetWAs_2[m]
    planetWA_2 = max(IWA, min(planetWA_2, OWA))
    planetWAs_2[m] = planetWA_2
    nsigmasens_finite, nsigmasens_inf, separation,\
    planet_implicit, speckle, zodi, cstab, bde, total =\
    fl.DRM_planetSens_vals(scenarioDF, perfLevel, CSprefix,\
                       isPhotonCounting, target, planetWA, lamD, detPCthreshold)

    separations_2.append(separation)
    cstabs.append(cstab)
    bdes.append(bde)
    totals.append(total)


FRN_df = pd.DataFrame(list(zip(separations_2,planetWAs_2,\
                               cstabs,bdes,totals)),\
                      columns = ['sep(mas)','WA','C Stab.','BDE','Total'])
FRN_df


# In[ ]:


fig, ax = plt.subplots(figsize=(8,4))

im1 = ax.plot(separations_2, cstabs,'--o')
im2 = ax.plot(separations_2, bdes,'--o')
im3 = ax.plot(separations_2, totals,'--o')
ax.set_title(f'FRN {scenario} [{lam}nm, {usableTinteg/3600} hrs, PP = {k_pp_CBE}x]')
ax.set_ylabel('Flux Ratio Noise,   ppb')
ax.set_xlabel('Planet Separation,    milli-arseconds')
plt.legend(['C Stab.','BDE','Total'], loc = 'upper right')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)  # Customize grid style
plt.show()

scenarioDF


# ### Contributions to N-sig Sensitivity
# 
# Note: Zodi contribution peaks at mid-was because of increasing occ trans with dec exo zodi flux as WA gets larger


# ### Contributions to N-sig Sensitivity

# In[ ]:


# User may select the upper and lower limits of range, but may default to 
# the inner and outer working angles (inclusive) of the range of WA 


numWAs = 7 # number of radial slices including the IWA and OWA
planetWAs_3 = np.linspace(IWA,OWA, numWAs)

nsigmasens_finites_3 = []

separations_3 = []
planet_implicits = []
speckles = []
zodis = []

for planetWA in planetWAs_3:

    print(planetWA)

    # Calculates the FRN sensitivity and contributors to the sensitivity

    nsigmasens_finite, nsigmasens_inf, separation,\
    planet_implicit, speckle, zodi, cstab, bde, total =\
    fl.DRM_planetSens_vals(scenarioDF, perfLevel, CSprefix,\
                           isPhotonCounting, target,planetWA, lamD, detPCthreshold)

    # For this radial slice, store the values for table

    # planet working angle (lam/D)
    # planetWAs_3.append(planetWA) 

     # separation (mas)
    separations_3.append(separation)

    # sensitivity
    nsigmasens_finites_3.append(nsigmasens_finite) 

    planet_implicits.append(planet_implicit)

    speckles.append(speckle)

    zodis.append(zodi)

Contributors_Sensitivity = pd.DataFrame(list(zip(separations_3,planetWAs_3,\
                                                 planet_implicits,speckles,\
                                                 zodis,nsigmasens_finites_3)),\
                   columns = ['sep(mas)','WA','planet','Spec','Zodi','~NEFR'])
Contributors_Sensitivity.round(2)

