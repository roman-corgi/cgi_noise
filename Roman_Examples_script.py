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

# In[6]:


scenFolder = fl.open_folder("EBcsvData","Scenarios")

# Specify here how many rows of the Scenario Excel to read:
scenarioDF = loadXLcol(scenFolder[scenario_filename],30).df

scenario = scenarioDF.at['Scenario','Latest']
scenarioDF

# In[7]:


lam = scenarioDF.at['CenterLambda_nm','Latest'] * uc.nm
print(lam)


lamD = lam / DPM

print(f'Lambda/D = {lamD:5.3e} rad')
print(f'Lambda/D = {lamD/uc.mas:5.3f} mas')


# In[5]:

#-------------------------------------------------------
# SET THE SCENARIO PARAMETERS
#-------------------------------------------------------

# User set host star and planet specifications
hostStar_vmag = 5.05
hostStar_dist = 13.8 # parsec
hostStar_type = 'g0v'

planet_phaseAngle = 65   # degrees
planet_sma_AU     = 4    # AU
planet_radius     = 1 # R_Jupiter
planet_GeomAlbedo = 0.3
exoZodi = 1  # X solar Zodi

sep_mas    = fl.Target.phaseAng_to_sep(planet_sma_AU, hostStar_dist, planet_phaseAngle)
print(f'Separation = {sep_mas:5.0f} mas')

# Set Study Parameters
SNRtarget  = 4.0 # desired signal to noise ratio
allocTinteg = 100 # integration time in hours
monthsAtL2 = 21  # months at L2
isPhotonCounting = True 
detPCthreshold = 5.0

pp_factor  = 2   # k_pp, post-processing factor

perfLevel = 'CBE'
# Contrast Stability options (choose 'MCBE_'):
CSprefix = 'MCBE_' # 'ICBE_'(instrument CBE), 'MCBE_'(mission CBE)
CBEtype = "mission", # "instrument", "mission"


intTimeDutyFactor = scenarioDF.at['DutyFactor_CBE','Latest']
 
# Integration Time in seconds
usableTinteg = allocTinteg * intTimeDutyFactor * uc.hour


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
opMode = scenarioDF.at['OPMODE_IMG_SPEC','Latest']
bandWidth = scenarioDF.at['BW','Latest']
RefStarVmag_CBE = scenarioDF.at['RefStar_V_mag_CBE','Latest']
RefStarDist = scenarioDF.at['RefStar_Distance_pc','Latest']
RefStarExoZodi = scenarioDF.at['RefStar_ExoZodi_Xsolar','Latest']
timeOnRef = scenarioDF.at['TimeonRefStar_tRef_per_tTar','Latest']
magLocalZodi = scenarioDF.at['LocZodi_magas2','Latest']
magExoZodi_1AU = scenarioDF.at['ExoZodi_magas2','Latest']
k_pp_CBE = scenarioDF.at['pp_Factor_CBE','Latest']
k_comp = 0
FWC_gr = 90000





# In[ ]:





# In[19]:


f_SR, CritLam, detPixSize_m, mpix, pixPlateSc\
    = fl.getFocalPlaneAttributes(opMode, scenarioDF, DET_CBE_Data, lam, bandWidth,\
                        DPM, CGdesignWL, omegaPSF)

inBandFlux0_sum, inBandZeroMagFlux, starFlux = fl.getSpectra(target, lam, bandWidth)

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


det_CIC_in_DC_units, det_CIC_gain  = fl.getDetectorCIC(DET_CBE_Data, detEMgain, missionFraction, frameTime)

CRhitsPerFrame,detPixAcross,detPixSize,CRrate,CRtailLen\
    = fl.getDetectorCosmicRays(perfLevel, DET_CBE_Data, detEMgain, frameTime )

CTE_clocking_efficiency, CTE_traps, signalPerPixPerFrame = fl.getCTE(DET_CBE_Data, rate_photoConverted,frameTime,missionFraction)

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
            signalPerPixPerFrame, detPixAcross, CRtailLen, CRhitsPerFrame, det_QE)

planetRate_proc, speckleRate_proc, zodiRate_proc,\
            ezo_bkgRate, lzo_bkgRate, residSpecRate\
    = fl.getNoiseVarianceRatesI( f_SR, starFlux, planetFlux, colArea, rawContrast, \
          thpt_t_pnt, thpt_t_speckle, thpt_tau_pk, CGintmpix, k_pp_CBE,\
                 dQE, rate_exoZodi_incPht, rate_loZodi_incPht,\
                     selDeltaC)


# In[ ]:


Kappa = SNRtarget / (f_SR*starFlux*colArea*thpt_t_pnt*dQE\
                      *usableTinteg)/uc.ppb

# strayLight, strayLight_FRN, strayLight_ph_pix_s,\
#     strayLight_ph_pix_h, strayLight_e_pix_h =\
#         fl.getStrayLightFRN(scenario, perfLevel, STRAY_FRN_Data, CG_Data, IWA, OWA,\
#                  opMode, DET_CBE_Data, ENF, detPixSize, mpix, dQE, Kappa,\
#                       usableTinteg, hostStar_type, inBandFlux0_sum,\
#                           f_SR, CGintmpix, colArea, thpt_t_speckle)

totNoiseVarRate, readNoiseRate, CIC_RNLK_noiseRate, darkNoiseRate, zodi_shot, speckle_shot, planet_shot=\
        fl.getNoiseVarianceRatesII(planetRate_proc, speckleRate_proc, zodiRate_proc,\
        ezo_bkgRate, lzo_bkgRate, residSpecRate, rate_planet_imgArea, ENF,\
            DarkCur_epoch_per_s, readNoise_leakage_in_current_units,\
            mpix, det_CIC_in_DC_units, readNoise, frameTime, isPhotonCounting,\
                k_sp, k_lzo, k_ezo, k_det, dQE)


# ## Calculations for Specified Planet
# Substitute for usableTinteg:
usableTmini = 0.71 * 60 * 60


# In[ ]:


instaPlanetDF1 = fl.DRMinstaplanet(planetRate_proc, usableTinteg, totNoiseVarRate, residSpecRate, SNRtarget)
timetoSNR = instaPlanetDF1.loc[0,'Value']
print(f'Time to {SNRtarget} SNR = {timetoSNR:5.2f} hrs')
round(instaPlanetDF1,2)
