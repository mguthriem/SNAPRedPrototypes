from mantid.simpleapi import *
from mantid import config
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import json
import time 
import importlib
sys.path.append('/SNS/SNAP/shared/Malcolm/code/SimpleCalibrationScripts/')
import SNAPTools as snp
importlib.reload(snp)

#SNAPReduce


#################################################################################
#20230929 - This is an attempt to modifiy an initial prototype reduction script,
# which incorporates containers into the workflow and uses SNAPRed functionality

#################################################################################
#BASIC SET UP
#################################################################################

runList =[59058]#,59056,59058,59061]#Comma separated list of run numbers. All must be in same state
isLite = True
downSampleFactor = 4

#container parameters
fullMaskPath = '/SNS/SNAP/IPTS-30264/shared/' #relative path from run IPTS typically will be shared directory but doesn't need to be
mskFiles = ['snap59053_dSpacing.json',
'snap59053_wavelength.json',
'SNAP59053_sides.xml',
'snap59058_dSpacing.json',
'snap59061_dSpacing.json']

##################################################################################
#ADVANCED SET UP
#################################################################################

#Flags
abbreviate = False #use apply a time filter when loading events 
subtractBackground = False #subtract an existing workspace
applyAbsorption = False
applyDiamWindow = False
saveGSAS = True
saveFullprof = False
DiagnosticMode = False
switchOffVan = False

applyDiamAtten = True
getMonitors = False

#dip stuff

L2_mon = 1.92 #m downstream distance of monitor
monitorNormWS = 'SNAP58813_monitors_scl_sm' #this is an empty inst data set multiplied by 0.0075 binned exactly as data.
upstream = 2 #which diamond is upstream
transDataPath = '/SNS/SNAP/IPTS-30264/shared/'

#configure
timeLimit=12 #opnly process this many hours if abbreviate = True

#subtractBackground
RunToSubtract='59014' #must have already been reduced

#Grouping
groupingList = ['Bank','Column']#,'Bank','All'] # Pixel grouping schemes to use

#Masking. A list of mask files can be provided and all will be applied
#The list can contain xml files (pixel masks, that discard entire pixels)
# or json files (Current implementation of swiss cheese masks)
# the json files are generated separately using 
# /SNS/SNAP/shared/Malcolm/devel/SandPitSNAP/ExtractMaskBinHistory.py (instructions in script)
fullMaskPath = '/SNS/SNAP/IPTS-30264/shared/' #relative path from run IPTS typically will be shared directory but doesn't need to be
mskFiles = ['snap59053_dSpacing.json',
'snap59053_wavelength.json',
'SNAP59053_sides.xml',
'snap59058_dSpacing.json',
'snap59061_dSpacing.json']

# Attenuation: Wavelength window correction (for dips), using "Sachith" method
# set parameters here. UB matrix file must exist. 
# The same Window is applied to all runs in runList 

windowWidthPars = [0.0165,0.0205]
UBPath = '/SNS/SNAP/IPTS-30264/shared/SNAP59053UB2.mat' #full file path to UB matrices output 

# Activate and et up sample and container absorption
sampleMaterial={'ChemicalFormula':'(Li7)6-P-S5-Br', 
    'PackingFraction':1.0, 
    'MassDensity':2.08}
sampleGeometry={'Shape':'Cylinder',
        'Height':0.3,
        'Radius':0.1,
        'Center':[0.,0.,0.]}

containerMaterial={'ChemicalFormula': 'Ti0.6765-Zr0.324',
    'NumberDensity': 0.100}
containerGeometry={'Shape':'HollowCylinder',
        'Height':0.3,
        'InnerRadius':0.1,
        'OuterRadius':0.11,
        'Center':[0.,0.,0.]}
        
# Output: Some parameters here controlling output files for Refinement

gsaTag = '' #just a label appended to GSAS powder data file name
if applyDiamWindow:
            gsaTag += '_LamWin'
customOutputBinning = False
customBinningParamsFile = '/SNS/SNAP/IPTS-30614/shared/SNAP58307_customBinningParams.json'



#############################################################################
# Data reduction workflow begins here. DON'T EDIT BELOW
#############################################################################
isLite = True
#attempt to auto assign state and calibration
stateID,stateDict,errorState=snp.StateFromRunFunction(runList[0])
if stateDict['wav']==1.8:
    
    stateCalibrationFolder = '/SNS/SNAP/shared/Calibration_Prototype/Powder/b810d6da5d4af06e/'
    geomCalibFile = f'{stateCalibrationFolder}SNAP058882_calib_geom_20230630.lite.h5'
    rawVMBFile = f'{stateCalibrationFolder}RVMB59000.lite.nxs'
elif stateDict['wav']==2.1:
    # print(f'Wavelength is: {stateDict["wav"]}')
    stateCalibrationFolder = '/SNS/SNAP/shared/Calibration_Prototype/Powder/04bd2c53f6bf6754/'
    geomCalibFile = f'{stateCalibrationFolder}SNAP058882_calib_geom_20230630.lite.h5'
    rawVMBFile = f'{stateCalibrationFolder}RVMB57473.lite.nxs'

if not DiagnosticMode:
    config.setLogLevel(0, quiet=True)

# Announce starting work
print('/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/')
print('/_/ SNAPRed Prototype v1.0 Starting...  /_/')
print('/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/')

print('\nINITIAL SET UP\n')
print(f'Configured for wavelength: {stateDict["wav"]}')
#STEP 0: Gather build all reduction parameter dictionaries

# File defining reduction parameters for this instrument state (v 1.4)
if extend:
    sPrmFile = f'{stateCalibrationFolder}CalibrationParameters_ext.json' 
else:
    sPrmFile = f'{stateCalibrationFolder}CalibrationParameters.json'
    
#Dict 0: sPRM (State-specific parameters)
with open(sPrmFile, "r") as json_file:
  sPrm = json.load(json_file) 

#STEP 1:  
# Create workspace names to hold raw vanadium data.
TOFrawVMB = f'_TOFrawVMB_{stateID}'
DSPrawVMB = f'_DSPrawVMB_{stateID}'
LAMrawVMB = f'_LAMrawVMB_{stateID}'

# Create lists to group workspaces
calibrationWorkspaces = []
if applyDiamWindow:
    diamondWindowWorkspaces = []

time0 = time.time() 
fSize = float(os.path.getsize(rawVMBFile))/1e9 
print(f'Loading raw vanadium correction (file {fSize:.3f} GB)')

#if not already loaded, load raw vanadium
try:
    tmp = mtd[DSPrawVMB]
    print(f'Vanadium for state: {stateID} already loaded')
except:
    LoadNexus(Filename=rawVMBFile,
        OutputWorkspace=TOFrawVMB)
    print(f' Raw vanadium loaded in: {time.time()-time0:.3f} s')

#apply diffraction calibration to raw vanadium ws 
ApplyDiffCal(InstrumentWorkspace=TOFrawVMB,
   CalibrationFile=geomCalibFile)    
# LoadInstrument(Workspace=wsVan,Filename=SNAPLite.xml
       
    
############################################################################################################################
# calc and apply wavelength windowing correction to VANADIUM: non-standard approach to remove diamond 
# dips without a d/stream detector#
############################################################################################################################      

#wavelength binning parameters (happily, group independent)
lamCen = sPrm["instrumentState"]["detectorState"]["wav"]
lamBand = sPrm["instrumentState"]["instrumentConfig"]["bandwidth"]
lamMin = lamCen - lamBand/2.0
lamMax = lamCen + lamBand/2.0
lamBin = -0.001 #TODO: check this is a good value

if applyDiamWindow:    
    time0 = time.time()
    print('Calculating diamond window')
    ConvertUnits(InputWorkspace=DSPrawVMB,
        Target='Wavelength',
        OutputWorkspace=LAMrawVMB)
    
    Rebin(InputWorkspace=LAMrawVMB,
        PreserveEvents=False,
        Params=f'{lamMin},{lamBin},{lamMax}', #binning approximately highest our of all PGS subgroups
        OutputWorkspace='tmp')
        
    snp.diamondWavelengthWindow('tmp',
        UBPath,
        True,
        windowWidthPars,
        'lambdaWindow') #lamdaWindow multiplies a ws (in wavelength units) to remove affected wavelengths.
        
    diamondWindowWorkspaces.extend(['lambdaWindow',
        'tmp_pks',
        'tmp_window_ticks',
        'tmp_window'])
        
    GroupWorkspaces(InputWorkspaces=diamondWindowWorkspaces,outputWorkspace='diamondWindowWorkspaces')    
    print(f' Calculated diamond window in: {time.time()-time0:.3f} s')
    
    #apply wavelength window to raw (unfocused) vanadium
    Multiply(LHSWorkspace=LAMrawVMB,
        RHSWorkspace='lambdaWindow',
        OutputWorkspace=LAMrawVMB)
        
    ConvertUnits(InputWorkspace=LAMrawVMB,
        Target='dSpacing',
        OutputWorkspace=DSPrawVMB)
               
else:
    ConvertUnits(InputWorkspace=TOFrawVMB,
        Target='dSpacing',
        OutputWorkspace=DSPrawVMB)
##############################################################

#Set up pixel grouping schemes.
#I got caught in mantid hell on this one. TODO: fix mantid generation/reading of pixel grouping schemes...

time0 = time.time()
print('Establishing pixel grouping schemes...')

overallDMin = []
overallDMax = []
overallDDoD = []
overallDBin = []
for igrp,group in enumerate(groupingList):
        
    if isLite:
        groupingFile=f'/SNS/SNAP/shared/Calibration/Powder/PixelGroupingDefinitions/SNAPFocGroup_{group}.lite.nxs'
        LoadNexusProcessed(Filename=groupingFile,
            OutputWorkspace=f'FocGrp_{group}_lite')
        ws = mtd[f'FocGrp_{group}_lite']
        calibrationWorkspaces.append(f'FocGrp_{group}_lite')
    else:
        groupingFile=f'/SNS/SNAP/shared/Calibration/Powder/PixelGroupingDefinitions/SNAPFocGroup_{group}.xml'
        LoadDetectorsGroupingFile(InputFile = groupingFile,
            OutputWorkspace=f'FocGrp_{group}')
        ws = mtd[f'FocGrp_{group}']
        calibrationWorkspaces.append(f'FocGrp_{group}')
    #collect info related to grouping workspaces.
    GroupIDs = ws.getGroupIDs()

    #instantiate grouping workspace
    snp.instantiateGroupingWS(sPrm,ws,isLite)

    #get groupingspecific parameters
    gPrm = snp.initGroupingParams(sPrm,ws,isLite)

    #Here, extract the "worst case" binning parameters for the entire group 
    overallDMin.append(min(gPrm["dMin"])) #smalled dMin of all subgroups
    overallDMax.append(max(gPrm["dMax"])) #largest dMax of all subgroups
    overallDDoD.append(abs(min(gPrm["delDOverD"]))) #smallest delDOverD 
    overallDBin.append( -1*abs(min(gPrm["delDOverD"]))/10 ) #finest log binning parameter of all subgroups

try:
    calibrationWorkspaces.extend([TOFrawVMB,DSPrawVMB,LAMrawVMB])
except:
    calibrationWorkspaces.extend([TOFrawVMB,DSPrawVMB])

# GroupWorkspaces(InputWorkspaces=calibrationWorkspaces,outputWorkspace='calibrationWorkspaces')
            
print(f' {len(groupingList)} pixel groups defined in: {time.time()-time0:.3f} s')
# for i,group in enumerate(groupingList):
#     print(f'{group} min: {overallDMin[i]} max: {overallDMax[i]}')

print('\nREDUCING DATA\n')
# Load and process all runs in list
for runNo in runList:
    
    time0=time.time()
    print(f'working on run {runNo}')
    
    #Dict 2: rPRM (Run-specific parameters #currently set to be the same for all runs)
    IPTSLoc = GetIPTS(RunNumber=runNo,Instrument='SNAP')
    print('creating mask location')
    print('iptsloc:',IPTSLoc)

    maskFileLoc = fullMaskPath
    print('mask location: ',maskFileLoc)
    
    rPrm={'maskFileName':mskFiles,
        'maskFileLoc': maskFileLoc,
        'saveGSAS':saveGSAS,
        'IPTSLoc':IPTSLoc,
        'upStreamDiamondUBPath':UBPath
        }
    
    runStr = str(runNo)
    wsName = f'SNAP{runStr}'
    wsName_lam = f'LAM{runStr}'
    wsName_dsp = f'DSP{runStr}'

    if isLite:
        NxsDataPath = rPrm['IPTSLoc'] +\
        sPrm['instrumentState']['instrumentConfig']['sharedDirectory'] +\
        'lite/' +\
        sPrm['instrumentState']['instrumentConfig']['nexusFilePrefix'] +\
        f'{runStr}.lite' +\
        sPrm['instrumentState']['instrumentConfig']['nexusFileExtension']
        
        if not os.path.exists(NxsDataPath):
            
            #need to create lite file
            time1 = time.time()
            print('    Lite version of file doesn\'t exist! Creating...') 
            #path to original file:
            origFile = rPrm['IPTSLoc'] +\
            sPrm['instrumentState']['instrumentConfig']['nexusDirectory'] +\
            sPrm['instrumentState']['instrumentConfig']['nexusFilePrefix'] +\
            runStr +\
            sPrm['instrumentState']['instrumentConfig']['nexusFileExtension']
            
            snp.makeLite(origFile,NxsDataPath)
            
            print(f'Created Lite file: {NxsDataPath} in: {time.time()-time1:.3f} s')
    else:
        NxsDataPath = rPrm['IPTSLoc'] +\
        sPrm['instrumentState']['instrumentConfig']['nexusDirectory'] +\
        sPrm['instrumentState']['instrumentConfig']['nexusFilePrefix'] +\
        runStr +\
        sPrm['instrumentState']['instrumentConfig']['nexusFileExtension']
        
    fSize = float(os.path.getsize(NxsDataPath))/1e9 
    
    print(f'    Loading raw event data (file {fSize:.3} GB)')   

    if abbreviate:
        LoadEventNexus(Filename=NxsDataPath, 
            OutputWorkspace=wsName, 
            FilterByTofMin=sPrm['instrumentState']['particleBounds']['tof']['minimum'], 
            FilterByTofMax=sPrm['instrumentState']['particleBounds']['tof']['maximum'],
            FilterByTimeStop=timeLimit*3600,
            NumberOfBins=500,
            LoadMonitors=getMonitors)
    else:
        print(NxsDataPath)
        LoadEventNexus(Filename=NxsDataPath, 
        OutputWorkspace=wsName, 
        FilterByTofMin=sPrm['instrumentState']['particleBounds']['tof']['minimum'], 
        FilterByTofMax=sPrm['instrumentState']['particleBounds']['tof']['maximum'], 
        NumberOfBins=500,
        LoadMonitors=getMonitors)

    NormaliseByCurrent(InputWorkspace=wsName,
        OutputWorkspace=wsName)
        
    CloneWorkspace(InputWorkspace=wsName,
        OutputWorkspace=f'TOF{runStr}_raw')


    if getMonitors:
        NormaliseByCurrent(InputWorkspace=wsName+'_monitors',
        OutputWorkspace=wsName+'_monitors')
        
        L_mon = sPrm['instrumentState']['instrumentConfig']['L1']+L2_mon
        
        TOFMin_mon = lamMin*L_mon/3.9561e-3
        TOFMax_mon = lamMax*L_mon/3.9561e-3
        
        LoadNexus(Filename=transDataPath+monitorNormWS+'.nxs',
            OutputWorkspace=monitorNormWS)
        
        Rebin(InputWorkspace=wsName+'_monitors',
            OutputWorkspace=wsName+'_monitors',
            Params=f'{TOFMin_mon},{lamBin},{TOFMax_mon}')        
            
        Divide(LHSWorkspace=wsName+'_monitors',
            RHSWorkspace=monitorNormWS,
            OutputWorkspace=wsName+'_monitors_ioio')
            
        ExtractSpectra(InputWorkspace=wsName+'_monitors_ioio',
            WorkspaceIndexList=1,
            OutputWorkspace=wsName+'_monitors_ioio')
                  

# apply diffraction calibration

    ApplyDiffCal(InstrumentWorkspace=wsName,
       CalibrationFile=geomCalibFile)


    
# apply absorption correction
################################################################################################################
#  Calculate absorption                                                                                                            #
################################################################################################################
    if applyAbsorption:
        gsaTag = '_A' #A label for this attenuation model
        
        time_0 = time.time()
        logger.notice('Starting to calculate absorption')
        ConvertUnits(InputWorkspace=wsName,
            OutputWorkspace=wsName_lam,
            Target='Wavelength')
            
        SetSample(wsName_lam,
            Geometry=sampleGeometry,
            Material=sampleMaterial,
            ContainerGeometry=containerGeometry,
            ContainerMaterial=containerMaterial)
        
        MonteCarloAbsorption(InputWorkspace=wsName_lam,
            OutputWorkspace='__AbsCorr',
            SimulateScatteringPointIn='SampleOnly')
        
        Divide(LHSWorkspace=wsName_lam,
            RHSWorkspace='__AbsCorr',
            OutputWorkspace=wsName)
            
        DeleteWorkspace(Workspace='__AbsCorr')
        logger.notice(f'Total time to execute atten corr: {time.time() - time_0:.3f}s')

################################################################################################################
#  Apply Diamond Attenuation                                                                                                            #
################################################################################################################
    if applyDiamAtten:

        time_0 = time.time()
        logger.notice('Starting to calculate absorption')
        
        LoadNexus(Filename=transDataPath+f'SNAP{runStr}_trns_diam{upstream}.nxs',
            OutputWorkspace='tmpDiamAtten')
        
        gsaTag = f'_DA{upstream}' #A label for this attenuation model
        
        ConvertUnits(InputWorkspace=wsName,
            OutputWorkspace=wsName_lam,
            Target='Wavelength')
        
        Rebin(InputWorkspace=wsName_lam,
            OutputWorkspace=f'{wsName_lam}_r',
            PreserveEvents=True,
            Params = f'{lamMin},{lamBin},{lamMax}')            
            
        Scale(InputWorkspace='tmpDiamAtten',
            OutputWorkspace='tmpDiamAtten',
            Factor=-1.0,
            Operation='Add')
                        
        Rebin(InputWorkspace='tmpDiamAtten',
            OutputWorkspace='tmpDiamAtten',
            Params = f'{lamMin},{lamBin},{lamMax}')
            
        Scale(InputWorkspace='tmpDiamAtten',
            OutputWorkspace='tmpDiamAtten',
            Factor=1.0,
            Operation='Add')
            
        ConvertToHistogram(InputWorkspace='tmpDiamAtten',
            OutputWorkspace='tmpDiamAtten')

        Divide(LHSWorkspace=f'{wsName_lam}_r',
            RHSWorkspace='tmpDiamAtten',
            OutputWorkspace=wsName)

################################################################################################################
#  Apply Masks to data and vanadium                                                                                                            #
################################################################################################################

# apply all masks to unfocused data and to vanadium
    print('    Masking')

    snp.SNAPMsk(wsName,rPrm)

    DSPlocalVMB='DSPlocalVMB'
    CloneWorkspace(InputWorkspace=DSPrawVMB,
        OutputWorkspace=DSPlocalVMB)
    
    snp.SNAPMsk(DSPlocalVMB,rPrm)    
    
################################################################################################################
#  Reduce data for each group in list of pixel grouping schemes                                                                                                            #
################################################################################################################
    
    for iGrp,focGrp in enumerate(groupingList):
        
        print(f'    working on group {focGrp}')
        #get groupingspecific parameters
        if isLite:
            gpws = mtd[f'FocGrp_{focGrp}_lite']
        else:
            gpws = mtd[f'FocGrp_{focGrp}']
            
        gPrm = snp.initGroupingParams(sPrm,gpws,isLite)


        if applyDiamWindow:
            
            ConvertUnits(InputWorkspace=wsName,
                Target='Wavelength',
                OutputWorkspace=wsName)
            
            Rebin(InputWorkspace=wsName,
                PreserveEvents=True,
                Params=f'{lamMin},{lamBin},{lamMax}', #binning 
                OutputWorkspace=wsName)
            
            Multiply(LHSWorkspace=wsName,
                RHSWorkspace='lambdaWindow',
                OutputWorkspace=wsName)
        
        # if applyDiamAtten:    
        #     # ConvertUnits(InputWorkspace=wsName,        
        #     ConvertUnits(InputWorkspace=wsName_lam+'_corr',
        #             Target='dSpacing',
        #             OutputWorkspace=wsName)
        # else:
        ConvertUnits(InputWorkspace=wsName,        
                Target='dSpacing',
                OutputWorkspace=wsName)  
        Rebin(InputWorkspace=wsName,
                PreserveEvents=True,
                Params=f'{overallDMin[iGrp]},{overallDBin[iGrp]},{overallDMax[iGrp]}', 
                OutputWorkspace=wsName)
        
        #first focus and process vanadium for this group
        # print(f'{iPrm["calibrationDirectory"]}{iPrm["pixelGroupingDirectory"]}{sPrm["focGroupDefinition"][iGrp]}')
        # print('groupingWOrkspace = ',f'SNAPGrp_{focGrp}')
        
        #First create Focused vanadium rebin and smooth

        if isLite:
            DiffractionFocussing(InputWorkspace=DSPlocalVMB,
                GroupingWorkspace=f'FocGrp_{focGrp}_lite',
                PreserveEvents = True,
                OutputWorkspace=f'{DSPlocalVMB}_{focGrp}')
        else:
            DiffractionFocussing(InputWorkspace=DSPlocalVMB,
                GroupingWorkspace=f'FocGrp_{focGrp}',
                PreserveEvents = True,
                OutputWorkspace=f'{DSPlocalVMB}_{focGrp}')
        
        SmoothData(Inputworkspace=f'{DSPlocalVMB}_{focGrp}',
            NPoints=20, #Temp fix
            OutputWorkspace=f'{DSPlocalVMB}_{focGrp}_sm')
        
        #Rebin Vanadium data to match group: need equal binning in all spectra
        Rebin(InputWorkspace=f'{DSPlocalVMB}_{focGrp}_sm', 
            PreserveEvents=True,
            Params=f'{overallDMin[iGrp]},{overallDBin[iGrp]},{overallDMax[iGrp]}', #common bin pars for group
            OutputWorkspace=f'{DSPlocalVMB}_{focGrp}_sm')
            

#         if len(sPrm["VPeaks"])!=0:
#             StripPeaks(InputWorkspace=f'{DSPlocalVMB}_{focGrp}', 
#                 OutputWorkspace=f'{DSPlocalVMB}_{focGrp}', 
#                 FWHM=5, 
#                 PeakPositions='1.2356,1.5133,2.1401', 
#                 BackgroundType='Quadratic')
#         
        # sPrm["VSmoothPoints"]='20' #over-ride
        
            
        #second focus data and divide by vanadium
        
        if isLite:
            DiffractionFocussing(InputWorkspace=wsName,
                GroupingWorkspace=f'FocGrp_{focGrp}_lite',
                PreserveEvents = False,
                OutputWorkspace=f'{wsName}_{focGrp}')
        else:
            DiffractionFocussing(InputWorkspace=wsName,
                GroupingWorkspace=f'FocGrp_{focGrp}',
                PreserveEvents = False,
                OutputWorkspace=f'{wsName}_{focGrp}')
                        
        if applyDiamWindow:
            outputWS = f'{wsName}_{focGrp}_nrm_wind'
        else:
            outputWS = f'{wsName}_{focGrp}_nrm'
            
        if applyDiamAtten:
            outputWS = outputWS + '_diamCorr'
            print('outputWS =',outputWS) 
            
        if extend:
            outputWS = outputWS + '_ext'
            
        if abbreviate:
            outputWS = outputWS + f'_{timeLimit}hrs'
            
        if len(mskFiles) == 0:
            outputWS = outputWS + 'noMsk'

        Rebin(InputWorkspace=f'{wsName}_{focGrp}', 
            PreserveEvents=True,
            Params=f'{overallDMin[iGrp]},{overallDBin[iGrp]},{overallDMax[iGrp]}', #common bin pars for group
            OutputWorkspace=f'{wsName}_{focGrp}')

        
        if switchOffVan:
            CloneWorkspace(InputWorkspace=f'{wsName}_{focGrp}',
                OutputWorkspace=outputWS)
        else:
            Divide(LHSWorkspace=f'{wsName}_{focGrp}',
                RHSWorkspace=f'{DSPlocalVMB}_{focGrp}_sm',
                OutputWorkspace=outputWS)      

        #subtract background if requested
        
        if subtractBackground:
            
            Rebin(InputWorkspace=f'SNAP{RunToSubtract}_{focGrp}_nrm', 
            PreserveEvents=True,
            Params=f'{overallDMin[iGrp]},{overallDBin[iGrp]},{overallDMax[iGrp]}', #common bin pars for group
            OutputWorkspace='Bgnd')
            
            
            Minus(LHSWorkspace=outputWS,
            RHSWorkspace='Bgnd',
            OutputWorkspace=outputWS + '_bs')
            
            outputWS = outputWS + '_bs'
            
        #finally rebin ragged   
        XBin = [-x*downSampleFactor/10 for x in gPrm["delDOverD"] ]
        RebinRagged(InputWorkspace=outputWS, 
            PreserveEvents=True,
            XMin=gPrm["dMin"],
            XMax=gPrm["dMax"],
            Delta=XBin,
            OutputWorkspace=outputWS) 
            
        if customOutputBinning:
            #Load custom binning parameters from file
            with open(customBinningParamsFile, "r") as json_file:
                custBin = json.load(json_file)
            gPrm["dMin"][iGrp] = custBin["dMin"][iGrp]
            gPrm["dMax"][iGrp] = custBin["dMax"][iGrp]  
            gPrm["delDOverD"][iGrp] = custBin["delDOverD"][iGrp]
             
        # print(f'group: {iGrp} {focGrp}')
        # print(gPrm)

        
        #If not in diagnostic mode, do a bunch of clean up...
        
        if not DiagnosticMode:
            DeleteWorkspaces([f'{DSPlocalVMB}_{focGrp}_sm',
                f'{DSPlocalVMB}_{focGrp}',
                f'{wsName}_{focGrp}'])
            try:
                DeleteWorkspace('pixMsk')
            except:
                pass
                
    #If requested save reduced data as GSAS-II files
        if saveGSAS:
            # print(f'    saving gsas-II file')
            ConvertUnits(InputWorkspace=outputWS,
                OutputWorkspace='gsasTOFdata',
                Target='TOF')
        
        
            gsasFileName = rPrm['IPTSLoc'] +\
                sPrm['instrumentState']['instrumentConfig']['reducedDataDirectory'] +\
                str(runNo)
                
            gsasFileName += focGrp
            gsasFileName += gsaTag
            gsasFileName += '.gsa'
            
            ws = mtd['gsasTOFdata']
            runInfo=ws.getRun()
            gsasHeader=runInfo.getLogData('run_title').value
            # rPrm.update({'runTitle':gsasFileName})
                            
            SaveGSS(InputWorkspace='gsasTOFdata',
                Filename=gsasFileName,
                SplitFiles=False,
                Format='SLOG',
                MultiplyByBinWidth=True,
                UseSpectrumNumberAsBankID=True,
                Append=False,
                UserSpecifiedGSASHeader=gsasHeader)
             
            # DeleteWorkspace(Workspace='gsasTOFdata')
            
        if saveFullprof:
            
            FPFileName = rPrm['IPTSLoc'] +\
                sPrm['instrumentState']['instrumentConfig']['reducedDataDirectory'] +\
                str(runNo)
                
            FPFileName += focGrp
            FPFileName += gsaTag
            FPFileName += '.fxye'
            
            ws = mtd[outputWS]
            runInfo=ws.getRun()
            # FPHeader=runInfo.getLogData('run_title').value
            # rPrm.update({'runTitle':gsasFileName})
                            
            SaveFocusedXYE(InputWorkspace=outputWS,
                Filename=FPFileName,
                SplitFiles=False,
                # Format='SLOG',
                # MultiplyByBinWidth=True,
                # UseSpectrumNumberAsBankID=True,
                Append=False,
                Format='XYE')
    
    RenameWorkspace(InputWorkspace=wsName,
        OutputWorkspace=wsName+'_unFocussed')
    
    if saveGSAS or saveFullprof:    
        print(f'Run {runNo} reduced in: {time.time()-time0:.3f} s (output files saved)')
    else:
        print(f'Run {runNo} reduced in: {time.time()-time0:.3f} s')
    
    if not DiagnosticMode:
        DeleteWorkspaces([f'{wsName}_unFocussed','DSPlocalVMB'])
    
    # Save Record of all reduction parameters

    redParRecordPath = IPTSLoc +\
                sPrm["instrumentState"]['instrumentConfig']["reductionRecordDirectory"]
    
    isExist = os.path.exists(redParRecordPath)
    if not isExist:
        os.makedirs(redParRecordPath)
    redParRecordPath += sPrm['instrumentState']['instrumentConfig']['name']
    redParRecordPath += str(runNo)
    redParRecordPath += gsaTag
    redParRecordPath += '.json'
    
    redParRecord = {}

    redParRecord.update(rPrm)
    redParRecord.update(sPrm)
    # if applyAttenuation:
    #     redParRecord.update(sampleMaterial)
    #     redParRecord.update(sampleGeometry)
    #     redParRecord.update(containerMaterial)
    #     redParRecord.update(containerGeometry)

    with open(redParRecordPath,'w+') as file:
        json.dump(redParRecord, file, indent=2)
        
    print('Reduction record written to file\n')
        
#Lastly tidy up reduced data
for focGrp in groupingList:
    GroupWorkspaces(GlobExpression=f'*{focGrp}*',
        OutputWorkspace=f'Reduced_{focGrp}')
        
print('\nEXECUTION COMPLETE\n')
DeleteWorkspace(Workspace='groupedWS')
config.setLogLevel(3, quiet=True)
