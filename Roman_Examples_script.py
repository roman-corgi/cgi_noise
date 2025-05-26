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

# In[4]:

# Scenario specification is done via an excel "scenario file"

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

# Primary Mirror Diameter for Roman
DPM = 2.363 * uc.meter

lam = scenarioDF.at['CenterLambda_nm','Latest'] * uc.nm
print(f"lam = {lam / uc.nm} nm")


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

indQE = QE_Data.df.loc[QE_Data.df['lambda_nm']<=(lam/uc.nm),'lambda_nm'].idxmax()
det_QE  = QE_Data.df.at[indQE,'QE_at_neg100degC']


# mission End of Life
detEOL_mos = DET_CBE_Data.df.at[0,'DetEOL_mos']
missionFraction = monthsAtL2 /detEOL_mos
    
# Detector dark current
detDarkBOM = DET_CBE_Data.df.at[0,'DarkBOM_e_per_pix_per_hr']
detDarkEOM = DET_CBE_Data.df.at[0,'DarkEOM_e_per_pix_per_hr']
DarkCur_epoch_per_hr = detDarkBOM+missionFraction*(detDarkEOM-detDarkBOM)
DarkCur_epoch_per_s = DarkCur_epoch_per_hr/3600


# In[14]:


# CGtauPol, indWA, CGcoreThruput, PSFpeakI, omegaPSF, CGintSamp,\
#     CGradius_arcsec, CGdesignWL, CGintmpix, CG_PSFarea_sqlamD, CGintensity,\
#         CG_occulter_transmission, CGcontrast\
#             = fl.coronagraph_pars(CG_Data, planetWA, IWA, DPM, lamD)


# #Coronagraph Parameters ---need planet WA and IWA and OWA

CGtauPol = 1
indWA = CG_Data.df[(CG_Data.df.rlamD <=planetWA)]['rlamD'].idxmax(axis=0)# index of last radial slice < planet's radius
 
# Frac light incident on primary in FWHM of PSF centered at r
CGcoreThruput = CG_Data.df.loc[indWA,'coreThruput']*CGtauPol

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

FocalPlaneAtt = fl.loadCSVrow(Path(filenamedir,'EBcsvData','CONST_SNR_FPattributes.csv'))
    
AmiciPar = fl.loadCSVrow(Path(filenamedir,'EBcsvData','CONST_Amici_parameters.csv'))
    
detPixSize_m = DET_CBE_Data.df.at[0,'PixelSize_m']

## Get Focal Plane Attributes for Selected Operating Mode

if opMode == "SPEC":   
    try:
        resolution = scenarioDF.at['R_required','Latest']
        print(f"resolution = {resolution}")
        f_SR = 1/(resolution*bandWidth)        
    except:
        print("R_required is not specified in scenario-- set to default 0")
        resolution = 0.0001
        f_SR = -1 

            
    CritLam = FocalPlaneAtt.df.at[1,'Critical_Lambda_m'] 
    mpix = fl.mpix_Amici(AmiciPar, lam, DPM,detPixSize_m,resolution)
    pixPlateSc = CritLam/DPM/2/uc.mas 
    
elif opMode == "IMG":
    f_SR = 1 
    CritLam = FocalPlaneAtt.df.at[0,'Critical_Lambda_m']
    mpix = omegaPSF * uc.arcsec**2 * (lam/CGdesignWL)**2*(2*DPM/CritLam)**2 
    pixPlateSc = CritLam/DPM/2/uc.mas 
else:
    raise Exception("Valid Operational Modes are IMG and SPEC")



# In[19]:


f_SR, CritLam, detPixSize_m, mpix, pixPlateSc\
    = fl.getFocalPlaneAttributes(opMode, scenarioDF, DET_CBE_Data, lam, bandWidth,\
                        DPM, CGdesignWL, omegaPSF)

# inBandFlux0_sum, inBandZeroMagFlux, starFlux = fl.getSpectra(target, lam, bandWidth)

SpectraFolder = os.path.join(filenamedir,'EBcsvData','Spectra')

#### Using hardcoded csv filename
SPECTRA_file = os.path.abspath(os.path.join(SpectraFolder,"SPECTRA_ALL_BPGS.csv"))
SPECTRA_Data = fl.loadCSVrow(SPECTRA_file)# Spectra data csv to dataframe

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

print(f"starFlux = {starFlux:3e}")


# In[17]:

# fluxRatio = fl.Target.alb_rad_sma_to_fluxRatio(target.albedo,\
#                                             target.radius_Rjup,\
#                                                 target.sma_AU)
fluxRatio = target.albedo * (target.radius_Rjup * uc.jupiterRadius / (target.sma_AU * uc.AU) )**2
planetFlux = fluxRatio * starFlux

# fluxRatio, planetFlux = fl.getFluxRatio(target, starFlux)

# Observing scenario is RDI with one observation of a brighter reference
# and one observation of the target.

RefStarAbsMag = RefStarVmag_CBE - 5*math.log10(RefStarDist/10)
RefStarinBandZeroMagFlux = inBandFlux0_sum.at[RefStarSpecType]
RefStarDeltaMag = target.v_mag - RefStarVmag_CBE
RefStarFlux = RefStarinBandZeroMagFlux*(10**((-0.4)*RefStarVmag_CBE))
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

RefStarAbsMag = RefStarVmag_CBE - 5*math.log10(RefStarDist/10)
RefStarinBandZeroMagFlux = inBandFlux0_sum.at[RefStarSpecType]
RefStarDeltaMag = target.v_mag - RefStarVmag_CBE
RefStarFlux = RefStarinBandZeroMagFlux*(10**((-0.4)*RefStarVmag_CBE))
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

# In[18]:

rate_planet_imgArea = f_SR * starFlux * fluxRatio * colArea * thpt_t_pnt * det_QE

rate_exoZodi_incPht = f_SR * exoZodiFlux * colArea * thpt_t_unif

rate_loZodi_incPht =  f_SR * loZodiFlux * colArea * thpt_t_unif 
  
rate_Zodi_imgArea = (rate_exoZodi_incPht + rate_loZodi_incPht)*det_QE

# speckleRate_imgArea
rate_speckleBkg = f_SR * starFlux * rawContrast * uc.ppb * thpt_tau_pk * \
    CGintmpix * colArea * thpt_t_speckle  * det_QE 
    
# total rate pixel with planet
rate_photoConverted = DarkCur_epoch_per_s +  (rate_planet_imgArea + rate_Zodi_imgArea + rate_speckleBkg )/mpix

rate_totalwithoutplanet =  DarkCur_epoch_per_s + (rate_Zodi_imgArea + rate_speckleBkg )/mpix



# rate_planet_imgArea, rate_Zodi_imgArea, rate_exoZodi_incPht,\
#     rate_loZodi_incPht, rate_speckleBkg, rate_photoConverted,\
#         rate_totalwithoutplanet \
#     = fl.getNoiseRates(f_SR, starFlux, fluxRatio, colArea, thpt_t_pnt, thpt_tau_pk,\
#               thpt_t_speckle, det_QE, exoZodiFlux, loZodiFlux, thpt_t_unif, k_comp,\
#               rawContrast,CGintmpix, mpix, DarkCur_epoch_per_s)

# frameTime, frameTime_ANLG,maxANLGt_fr,maxPCt_fr, detEMgain\
#     = fl.getFrameExposureTime(DET_CBE_Data, FWC_gr, rate_totalwithoutplanet,\
#                      rate_photoConverted, isPhotonCounting)

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

# In[18]:


det_CIC_in_DC_units, det_CIC_gain  = fl.getDetectorCIC(DET_CBE_Data, detEMgain, missionFraction, frameTime)

CRhitsPerFrame,detPixAcross,detPixSize,CRrate,CRtailLen\
    = fl.getDetectorCosmicRays(perfLevel, DET_CBE_Data, detEMgain, frameTime )

CTE_clocking_efficiency, CTE_traps, signalPerPixPerFrame = fl.getCTE(DET_CBE_Data, rate_photoConverted,frameTime,missionFraction)

hotPixFrac, hotPix = fl.gethotPixels(DET_CBE_Data, missionFraction)

ENF = fl.getENF(isPhotonCounting)

detCamRead = DET_CBE_Data.df.at[0,'ReadNoise_e']

# readNoise, readNoise_leakage, readNoise_leakage_in_current_units,\
#      PCeffloss, readNoise_w_gain\
#    = fl.getReadNoiseandPCeffloss(detCamRead, detPCthreshold, isPhotonCounting, frameTime, detEMgain)

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

# In[18]:

det_CTE = CTE_clocking_efficiency * CTE_traps

det_PC_threshold_efficiency = 1 - PCeffloss

signal_region_electron_rate = signalPerPixPerFrame * det_CTE

# Photon-counting coincidence efficiency 
det_PC_coincid_effic = (1 - math.exp(-signal_region_electron_rate))\
    / (signal_region_electron_rate)

# Efficiency after subtracting fraction lost to hot pixels
det_hotPix = 1 - hotPix

det_cosmicRays = 1 - CRhitsPerFrame * CRtailLen/(detPixAcross**2)

dQE = det_QE * det_CTE * det_PC_threshold_efficiency * det_PC_coincid_effic * det_hotPix * det_cosmicRays


# In[ ]:

planetRate_proc = f_SR*planetFlux*colArea*thpt_t_pnt*dQE

speckleRate_proc = f_SR * starFlux * rawContrast * thpt_tau_pk * CGintmpix * thpt_t_speckle * colArea * dQE * uc.ppb

residSpecRate =  f_SR * starFlux * (selDeltaC/k_pp_CBE) * thpt_tau_pk * CGintmpix * thpt_t_speckle * colArea * dQE

Kappa = SNRtarget / (f_SR*starFlux*colArea*thpt_t_pnt*dQE\
                      *usableTinteg)/uc.ppb

darkNoiseRate = mpix * DarkCur_epoch_per_s * ENF**2


if isPhotonCounting:
    CIC_RNLK_noiseRate = ENF**2*mpix*(det_CIC_in_DC_units +\
                                      readNoise_leakage_in_current_units)
else:
    CIC_RNLK_noiseRate = ENF**2*mpix*(det_CIC_in_DC_units)

readNoiseRate = (mpix/frameTime) * readNoise**2   

totNoiseVarRate = ENF**2 * (rate_planet_imgArea +\
                          k_sp * speckleRate_proc +\
                          k_lzo * rate_loZodi_incPht * dQE +\
                          k_ezo * rate_exoZodi_incPht * dQE ) +\
                    k_det * (darkNoiseRate + CIC_RNLK_noiseRate) +\
                    k_det * readNoiseRate

# ## Calculations for Specified Planet
# Substitute for usableTinteg:
usableTmini = 0.71 * 60 * 60


# In[ ]:


instaPlanetDF1 = fl.DRMinstaplanet(planetRate_proc, usableTinteg, totNoiseVarRate, residSpecRate, SNRtarget)
timetoSNR = instaPlanetDF1.loc[0,'Value']
print(f'Time to {SNRtarget} SNR = {timetoSNR:5.2f} hrs')
round(instaPlanetDF1,2)
