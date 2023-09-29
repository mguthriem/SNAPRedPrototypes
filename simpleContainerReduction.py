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

#################################################################################s
#HACKS:
#
# 20230709: 1) adding modifications to support dip fitting/correction.
#           2) modified to retain raw vanadium, which is slow to load
#           3) debugging event/histogram vanadium normalisation issues
#
# 20230927: 1) attempt to standardize for "simple" containers
 
#################################################################################


#################################################################################
#BASIC SET UP
#################################################################################

runList =[59041]#59042 #Comma separated list of run numbers. All must be in same state

downSampleFactor = 4
#Grouping
groupingList = ['All','Bank','Column']#,'Bank','All'] # Pixel grouping schemes to use
#Container
container = 'quartz-capillary' 	#Options:
			       	#'kapton-capillary'
                               	#'single-toroid-TiZr-gasket'
  			       	#'quartz-capillary'
				#'none' (if none, must specify variable sampleRadius

#sampleRadius = 0.1 #cm

sampleMaterial={'ChemicalFormula':'C5-N-H4-D0.47-Mo-O3', 
	    'PackingFraction':0.67,#TODO 
	    'MassDensity':0.0} #TODO

##################################################################################
#ADVANCED SET UP
#################################################################################

#############
#Flags      #
#############
abbreviate = False
subtractBackground = False

# absorption/attenuation settings.
applyAbsorption = False #mass absorption for cylindrical geometry
applyDiamWindow = False #removes wavelengths affected by dips
applyDiamAtten = False #full diamond dip correction
getMonitors = False #used during set up of full diamond dip correction

#output options
saveGSAS = True
saveFullprof = False

#diagnostic/testing options
DiagnosticMode = True
switchOffVan = False

#############
#Parameters #
#############

#dip stuff

L2_mon = 1.92 #m downstream distance of monitor
monitorNormWS = 'SNAP58813_monitors_scl_sm' #this is an empty inst data set multiplied by 0.0075 binned exactly as data.
upstream = 2 #which diamond is upstream
transDataPath = '/SNS/SNAP/IPTS-30264/shared/'

#configure
timeLimit=12 #opnly process this many hours if abbreviate = True

#subtractBackground
RunToSubtract='59014' #must have already been reduced



#Masking. A list of mask files can be provided and all will be applied
#The list can contain xml files (pixel masks, that discard entire pixels)
# or json files (Current implementation of swiss cheese masks)
# the json files are generated separately using 
# /SNS/SNAP/shared/Malcolm/devel/SandPitSNAP/ExtractMaskBinHistory.py (instructions in script)
fullMaskPath = '/SNS/SNAP/IPTS-31732/shared/' #relative path from run IPTS typically will be shared directory but doesn't need to be
mskFiles = []

# Attenuation: Wavelength window correction (for dips), using "Sachith" method
# set parameters here. UB matrix file must exist. 
# The same Window is applied to all runs in runList
# These are  ignored if applyDiamWindow is false 

windowWidthPars = [0.0165,0.0205]
UBPath = '/SNS/SNAP/IPTS-30264/shared/SNAP59053UB2.mat' #full file path to UB matrices output 

# Set up sample and container absorption

# 
if container == 'kapton-capillary':
	containerMaterial={'ChemicalFormula': 'C22-H10-N2-O5',
	    'MassDensity': 1.42}
	containerGeometry={'Shape':'HollowCylinder',
		'Height':0.3,
		'InnerRadius':0.1, #2mm diameter
		'OuterRadius':0.1013, #13 um wall thickness
		'Center':[0.,0.,0.]}
elif container == 'single-toroid-TiZr-gasket':
#PE Gasket
	containerMaterial={'ChemicalFormula': 'Ti0.6765-Zr0.324',
	     'NumberDensity': 0.100}
	containerGeometry={'Shape':'HollowCylinder',
		 'Height':0.3,
		 'InnerRadius':0.3,
		 'OuterRadius':0.7,
		 'Center':[0.,0.,0.]}
elif container == 'quartz-capillary':
	containerMaterial={'ChemicalFormula': 'SiO2',
	     'MassDensity': 2.65}
	containerGeometry={'Shape':'HollowCylinder',
		 'Height':0.3,
		 'InnerRadius':0.1,
		 'OuterRadius':0.12, #2mm 200um thick container
		 'Center':[0.,0.,0.]}
elif container == 'none':
	continue	

if container =='none:
	sampleGeometry={'Shape':'Cylinder',
		'Height':0.3,
		'Radius':sampleRadius,
		'Center':[0.,0.,0.]}
else:
	sampleGeometry={'Shape':'Cylinder',
		'Height':0.3,
		'Radius':containerGeometry["InnerRadius"],
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
    # RVMB59000_3_noAbs.lite.nxs
    stateCalibrationFolder = '/SNS/SNAP/shared/Calibration_Prototype/Powder/b810d6da5d4af06e/'
    geomCalibFile = f'{stateCalibrationFolder}SNAP058882_calib_geom_20230630.lite.h5'
    rawVMBFile = f'{stateCalibrationFolder}RVMB59000_dsp.lite.nxs'
elif stateDict['wav']==2.1:
    
    stateCalibrationFolder = '/SNS/SNAP/shared/Calibration_Prototype/Powder/04bd2c53f6bf6754/'
    geomCalibFile = f'{stateCalibrationFolder}SNAP058882_calib_geom_20230630.lite.h5'
    rawVMBFile = f'{stateCalibrationFolder}RVMB58810_dsp.lite.nxs'
    print(f'Wavelength is: {stateDict["wav"]}')
    print(f'DiffractionCalibration: {geomCalibFile}')
    print(f'Raw Vanadium: {rawVMBFile}')
elif stateDict['wav']==6.4:
    # print(f'Wavelength is: {stateDict["wav"]}')
    stateCalibrationFolder = '/SNS/SNAP/shared/Calibration_Prototype/Powder/c073719d9101e8f2/'
    geomCalibFile = f'{stateCalibrationFolder}SNAP058882_calib_geom_20230630.lite.h5'
    rawVMBFile = f'{stateCalibrationFolder}RVMB59016_dsp.lite.nxs'

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
        OutputWorkspace=DSPrawVMB)
    print(f' Raw vanadium loaded in: {time.time()-time0:.3f} s')
     
    
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
    
    # Rebin(InputWorkspace=LAMrawVMB,
    #     PreserveEvents=False,
    #     Params=f'{lamMin},{lamBin},{lamMax}', #binning approximately highest our of all PGS subgroups
    #     OutputWorkspace='tmp')
        
    snp.diamondWavelengthWindow(LAMrawVMB,
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
               

##############################################################

#Set up pixel grouping schemes.
#I got caught in mantid hell on this one. TODO: fix mantid generation/reading of pixel grouping schemes...

time0 = time.time()
print('Establishing pixel grouping schemes...')

# overallDMin = []
# overallDMax = []
# overallDDoD = []
# overallDBin = []
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


try:
    calibrationWorkspaces.extend([DSPrawVMB])
except:
    calibrationWorkspaces.extend([DSPrawVMB,LAMrawVMB])

GroupWorkspaces(InputWorkspaces=calibrationWorkspaces,outputWorkspace='calibrationWorkspaces')
            
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
    # wsName = f'SNAP{runStr}'
    wsName_lam = f'LAM{runStr}'
    wsName_dsp = f'DSP{runStr}'

    snp.loadAndPreProcess(runNo,sPrm,isLite,geomCalibFile,wsName_dsp,getMonitors)
        
    CloneWorkspace(InputWorkspace=wsName_dsp,
        OutputWorkspace=f'DSP{runStr}_unfocussed') #

    # if getMonitors:
    #     NormaliseByCurrent(InputWorkspace=f'SNAP{str(runNo)}_monitors',
    #     OutputWorkspace=f'SNAP{str(runNo)}_monitors')
    #     
    #     L_mon = sPrm['instrumentState']['instrumentConfig']['L1']+L2_mon
    #     
    #     TOFMin_mon = lamMin*L_mon/3.9561e-3
    #     TOFMax_mon = lamMax*L_mon/3.9561e-3
    #     
    #     LoadNexus(Filename=transDataPath+monitorNormWS+'.nxs',
    #         OutputWorkspace=monitorNormWS)
    #     
    #     Rebin(InputWorkspace=wsName+'_monitors',
    #         OutputWorkspace=wsName+'_monitors',
    #         Params=f'{TOFMin_mon},{lamBin},{TOFMax_mon}')        
    #         
    #     Divide(LHSWorkspace=wsName+'_monitors',
    #         RHSWorkspace=monitorNormWS,
    #         OutputWorkspace=wsName+'_monitors_ioio')
    #         
    #     ExtractSpectra(InputWorkspace=wsName+'_monitors_ioio',
    #         WorkspaceIndexList=1,
    #         OutputWorkspace=wsName+'_monitors_ioio')
                  


    
# apply absorption correction
################################################################################################################
#  Calculate absorption                                                                                                            #
################################################################################################################
    if applyAbsorption:
        gsaTag = '_A' #A label for this attenuation model
        
        time_0 = time.time()
        logger.notice('Starting to calculate absorption')
        ConvertUnits(InputWorkspace=wsName_dsp,
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
            
        ConvertUnits(InputWorkspace=wsName_lam,
            OutputWorkspace=wsName_dsp,
            Target='dSpacing')
            
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
        
        ConvertUnits(InputWorkspace=wsName_dsp,
            OutputWorkspace=wsName_lam,
            Target='Wavelength')
        
        # Rebin(InputWorkspace=wsName_lam,
        #     OutputWorkspace=f'{wsName_lam}_r',
        #     PreserveEvents=False,
        #     Params = f'{lamMin},{lamBin},{lamMax}')            
            
        Scale(InputWorkspace='tmpDiamAtten',
            OutputWorkspace='tmpDiamAtten',
            Factor=-1.0,
            Operation='Add')
                        
        RebinToWorkspace(InputWorkspace='tmpDiamAtten',
            OutputWorkspace='tmpDiamAtten',
            WorkspaceToMatch=wsName_lam)
            
        Scale(InputWorkspace='tmpDiamAtten',
            OutputWorkspace='tmpDiamAtten',
            Factor=1.0,
            Operation='Add')
            
        ConvertToHistogram(InputWorkspace='tmpDiamAtten',
            OutputWorkspace='tmpDiamAtten')

        Divide(LHSWorkspace=wsName_lam,
            RHSWorkspace='tmpDiamAtten',
            OutputWorkspace=wsName_lam)
            
        ConvertUnits(InputWorkspace=wsName_lam,
            OutputWorkspace=wsname_dsp,
            Target='dSpacing')

################################################################################################################
#  Apply Masks to data and vanadium                                                                                                            #
################################################################################################################

# apply all masks to unfocused data and to vanadium
    print('    Masking')

    snp.SNAPMsk(wsName_dsp,rPrm,sPrm)

    DSPlocalVMB='DSPlocalVMB'
    CloneWorkspace(InputWorkspace=DSPrawVMB,
        OutputWorkspace=DSPlocalVMB)
    
    snp.SNAPMsk(DSPlocalVMB,rPrm,sPrm)    
    
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
            
            ConvertUnits(InputWorkspace=wsName_dsp,
                Target='Wavelength',
                OutputWorkspace=wsName_lam)
            
            
            Multiply(LHSWorkspace=wsName_lam,
                RHSWorkspace='lambdaWindow',
                OutputWorkspace=wsName)
                
            ConvertUnits(InputWorkspace=wsName_lam,
                Target='dSpacing',
                OutputWorkspace=wsName_dsp)
        
       
        #first focus and process vanadium for this group

        if isLite:
            DiffractionFocussing(InputWorkspace=DSPlocalVMB,
                GroupingWorkspace=f'FocGrp_{focGrp}_lite',
                PreserveEvents = False,
                OutputWorkspace=f'{DSPlocalVMB}_{focGrp}')
        else:
            DiffractionFocussing(InputWorkspace=DSPlocalVMB,
                GroupingWorkspace=f'FocGrp_{focGrp}',
                PreserveEvents = False,
                OutputWorkspace=f'{DSPlocalVMB}_{focGrp}')
        
        SmoothData(Inputworkspace=f'{DSPlocalVMB}_{focGrp}',
            NPoints=20, #Temp fix
            OutputWorkspace=f'{DSPlocalVMB}_{focGrp}_sm')

        #second focus data and divide by vanadium
        
        if isLite:
            DiffractionFocussing(InputWorkspace=wsName_dsp,
                GroupingWorkspace=f'FocGrp_{focGrp}_lite',
                PreserveEvents = False,
                OutputWorkspace=f'{wsName_dsp}_{focGrp}')
        else:
            DiffractionFocussing(InputWorkspace=wsName_dsp,
                GroupingWorkspace=f'FocGrp_{focGrp}',
                PreserveEvents = False,
                OutputWorkspace=f'{wsName_dsp}_{focGrp}')
                        
        if applyDiamWindow:
            outputWS = f'{wsName_dsp}_{focGrp}_nrm_wind'
        else:
            outputWS = f'{wsName_dsp}_{focGrp}_nrm'
            
        if applyDiamAtten:
            outputWS = outputWS + '_diamCorr'
            print('outputWS =',outputWS) 

            
        if abbreviate:
            outputWS = outputWS + f'_{timeLimit}hrs'
            
        if len(mskFiles) == 0:
            outputWS = outputWS + 'noMsk'
       
        if switchOffVan:
            CloneWorkspace(InputWorkspace=f'{wsName_dsp}_{focGrp}',
                OutputWorkspace=outputWS)
        else:
            print(f'trying to divide {wsName_dsp}_{focGrp} by {DSPlocalVMB}_{focGrp}_sm')
            Divide(LHSWorkspace=f'{wsName_dsp}_{focGrp}',
                RHSWorkspace=f'{DSPlocalVMB}_{focGrp}_sm',
                OutputWorkspace=outputWS)      

        #subtract background if requested
        
        if subtractBackground:
            
            
            Minus(LHSWorkspace=outputWS,
            RHSWorkspace=f'SNAP{RunToSubtract}_{focGrp}_nrm',
            OutputWorkspace=outputWS + '_bs')
            
            outputWS = outputWS + '_bs'
            
        #finally rebin ragged   
        XBin = [-x*downSampleFactor/10 for x in gPrm["delDOverD"] ]
        RebinRagged(InputWorkspace=outputWS, 
            PreserveEvents=False,
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
                f'{wsName_dsp}_{focGrp}'])
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
            rPrm.update({'runTitle':gsasFileName})
                            
            SaveGSS(InputWorkspace='gsasTOFdata',
                Filename=gsasFileName,
                SplitFiles=False,
                Format='SLOG',
                MultiplyByBinWidth=True,
                UseSpectrumNumberAsBankID=True,
                Append=False,
                UserSpecifiedGSASHeader=gsasHeader)
             
            DeleteWorkspace(Workspace='gsasTOFdata')
            print(f'reduced data written to: {gsasFileName}')
            
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
    
    RenameWorkspace(InputWorkspace=wsName_dsp,
        OutputWorkspace=wsName_dsp+'_unFocussed')
    
    if saveGSAS or saveFullprof:    
        print(f'Run {runNo} reduced in: {time.time()-time0:.3f} s (output files saved)')
    else:
        print(f'Run {runNo} reduced in: {time.time()-time0:.3f} s')
    
    if not DiagnosticMode:
        DeleteWorkspaces([f'{wsName_dsp}_unFocussed','DSPlocalVMB'])
    
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
