# todo: taper off correction array

tileNum = 0 # 0 to ignore. ONLY WORKS ON SERVERS. For running tiled bands to be stitched together later.
numTiles = 8 # number of horizontal tiles to break into if tileNum > 0.
#---------------------------------------------------------------------
# INPUTS #
#---------------------------------------------------------------------
options = {}

# MOVING WINDOW AND TARGET BLOCK SIZES ####
options['radius'] = 100# In PIXELS. Search radius, with sources activated within the radius and outside of the center (target) block.
options['blockSize'] = 15 # Odd number. Targets will be square blocks of pixels with this number of pixels on a side.

options['projectDir'] = r'C:\DATADRIVE\DUKE_PNW_DATA\PNW_Omniscape_1080m'#r'D:\Users\bmcrae\Duke_PNW_Omniscape\d8_omniscape' #r'C:\DATADRIVE\DUKE_PNW_DATA\PNW_OmniScape'# this is where all the input data are, and where out directory will be created.
options['resisRasterBase'] = 'rDraft7_1080mclip.tif'#'rNullFadeTestClip2.tif'#'ones.tif'#
options['outputDirBase'] = 'correctionTestBS15'

options['useSourceRaster'] = False #Use a layer specifying source/target pixels to connect. 0 = no source, anything positive will indicate strength of a source at that pixel. If not using a separate source raster, will use r Cutoff and consider anything with resistance lower than this to be a source/target
options['sourceRasterBase'] = 'te9x_sel82_sm_ones.tif'#'DUKE_CIRCUITSCAPE_resistances_draft4_Habitat_min_540m_clip.tif'#'srcInvR.tif'# Name of source raster, if using one
options['sourceInverseR'] = False # Source strengths are inverse of input resistances (before squaring)

# RESISTANCES AND CUTOFF VALUES ####
options['rCutoff'] = 100# If NOT using a source raster, everything <= this value in the ORIGINAL resistance raster will be a source (before squaring if squareResistances is True)
options['squareResistances'] = False # Will square resistance raster before doing any calculations. sourceInverseR and rCutoff applied before squaring.

options['squareRootInputs'] = False


# CALCULATE MAXIMUM INSTEAD OF SUM OF CURRENT ACROSS ITERATIONS
options['calcMaxCur'] = False #Note: this has more tiling/fading issues.

# CLIMATE ####
options['useClimate'] = False
options['matchClimatePCs'] = True  # climate method- principle components if True, present-day temperature differences if False
# parameters to match temperature differences
options['tDiff'] = 4 # Match pixels that differ in temperature by this value
options['tWindow'] = 1 #How close does target need to be to the tDiff value
options['climateRasterBase'] = 'tmean_adks.asc'#'6x6t1pc2.asc'#'6x6clim.asc'#'TMEAN_NE_clip.tif'
options['absClim'] = False # connect if ABS VAL of climate differs by tcutoff. Meant to help with fade.
# parameters to match current and future PCs
options['t1PC1RasterBase'] = 'NORM_6190_PC1_PNWmockup.tif'
options['t1PC2RasterBase'] = 'NORM_6190_PC2_PNWmockup.tif'
options['t2PC1RasterBase'] = 'MIROC5_2080s_RCP85_PC1_PNWmockup.tif'
options['t2PC2RasterBase'] = 'MIROC5_2080s_RCP85_PC2_PNWmockup.tif'
options['PCWindow'] = 0.9 # Euclidean distance between PCs to determine match

# DISTANCE FUNCTION#### 
options['useDistanceFunction'] = False # If true can set minimum and maximum distances or function between source and target pixels. Sources outside this range won't be activated.
options['distEq'] = "(options['radius']-dist)/options['radius']" # Distance equation for source strengths. Use NONE to ignore

options['minDist'] = None # In Pixels. Sources closer than this will have no current injected. Use None to ignore.
options['maxDist'] = None # Sources farther than this will have no current injected. Use None to ignore.
# Note: current can actually be higher in some pixels when min distances are used. Will be common with targetOnly mode because same amount
# of current is injected farther from target. But can also happen when opposing (canceling) currents occur without min dist.

# FLOW ACCUMULATION CALCULATIONS ####
options['calcFA'] = False
options['addRandomResistances'] = True #add random values to FA resis raster
options['limitCalcsByCWD'] = False
options['cwdLimit'] = options['radius'] 

options['useCwDistanceFunction'] = False # NOT COMPLETED/TESTED YET
options['cwDistEq'] = "(cwdLimit-cwDist)/cwdLimit" # Distance equation for source strengths. Use NONE to ignore

# CURRENT CALCULATIONS ####
options['calcCurrent'] = True

# OPTION TO LIMIT ANALYSIS EXTENT ####  
# Bands are horizontal, with width equal to blockSize. There are approx nrows/blockSize bands in a raster. Stripes are vertical, with width equal to blocksize.
# These options allow you to only process a subset of bands and stripes. 
# options['startBand'] = 2675# First band to process. Use 0 to ignore.
# options['endBand'] = 2675# # Stops after processing this band. Use 0 to ignore.
# options['startStripe'] = 2227 # First stripe to process. Use 0 to ignore.
# options['endStripe'] = 2227 # Stops after processing this stripe. 0 to ignore.
options['startBand'] = 0#4#16#4# First band to process. Use 0 to ignore.
options['endBand'] = 0#15# # Stops after processing this band. Use 0 to ignore.
options['startStripe'] = 0 # First stripe to process. Use 0 to ignore.
options['endStripe'] = 0
 # Stops after processing this stripe. 0 to ignore.

                                                              
# VOLTAGES #### Note: code not complete yet!
options['calcVoltages'] = True
options['adjustVoltages'] = False

# TARGETS AND WEIGHTING ####
options['weightTargOnly'] = True # Total current flow is ~ ntargs. May make sense if # dispersers limited, or number of dispersers a pixel can accept is limited. If false and noweight is false, current flow will be ~ntargs*nsources
# Recommend not changing following
options['noWeight'] = False # Recommend False. 
options['centerGround'] = True# Recommend True.
options['negTargets'] = False # Recommend False. negative sources at targets- blocks can work better with centerground, fade out, and no neg targs.


# FADE CURRENTS #FIXME: check out divide by zero in fade calc
options['fadeIn'] = False# Decreases current from each solve linearly with distance to center. Helps with affects of blocks. Best with center ground.
from math import *
options['fadeInDist'] = (options['blockSize'])#(options['blockSize'])#-1)/2#/4.0#sqrt(2)*(options['blockSize']/2.0)#options['radius']/2 #
options['fadeConstant'] = 0.5 # Constant added to numerator and denominator in fade equation. Avoids zero current at center ground. Higher values mean more current at ground 
options['setKernelMean'] = True

options['removeArtifacts'] = True

# Quantilize maps ####
options['quantilizeMaps'] = True # Additional current map will be written with values binned by quantile (e.g., 0-100 with 100 quantiles). Makes display easier (just use min-max stretch).
options['numQuantiles'] = 100

# SPECIAL FUNCTIONS
options['calcNull'] = False # Calculates a 'null' result with all resistances = 1. 
options['calcFANull'] = False
options['saveSourceAndTargetCounts'] = False # #Save source and target counts when doing climate analyses
options['doSourceAndTargetCountsOnly'] = False # if saving counts, DON'T calculate current
options['printTimings'] = True

options['cleanUpBandFiles'] = True
    
    
#---------------------------------------------------------------------
# END INPUTS #
#---------------------------------------------------------------------
if tileNum > 0:
    print 'Tiling with ' + str(numTiles) + ' tiles. Running tile #' + str(tileNum) + ' in this instance.'
    print 'Any manually set start and end bands will be ignored.'


from functools import wraps
PROF_DATA = {}
import arcgisscripting
import os 
import arcpy
import math
import ConfigParser
import string
import os.path as path
import glob
import time
import datetime
import subprocess
import gc
import sys
import numpy as npy
import shutil
import glob
# from lm_retry_decorator import retry

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

def omniscape(options):
    try:
        theStart = datetime.datetime.now()  
        # TILING
        options['tileNum'] = tileNum
        options['numTiles'] = numTiles

        options = set_options_and_dirs(options)   
        copy_this_file(options)
        
        arcpy.env.scratchWorkspace = options['scratchDir']
        arcpy.env.Workspace = options['scratchDir']
        os.environ["TEMP"] = options['scratchDir']
        os.environ["TMP"] = options['scratchDir']    
        
        #FIXME: check correction arrays, especially positioning
        if options['calcCurrent'] and options['removeArtifacts']:
            print 'Calculating null results to remove artifacts'
            artifactArray = npy.ones((2*options['radius']+1,2*options['radius']+1),dtype = 'int32')
            artifactHeader = createArtifactHeader(options)

            subsetCenterRow = centerRow =subsetCenterCol=centerCol= options['radius']
            artifactBS1Result = get_null_result(artifactArray, artifactHeader, subsetCenterRow, subsetCenterCol, blockSize = 1)
            artifactBSXResult = get_null_result(artifactArray, artifactHeader, subsetCenterRow, subsetCenterCol, blockSize = options['blockSize'])
            # then create larger grid, add in XxX nullBS1Result results
#FIXME CHECK:            
            artifactBS1SumArray = npy.zeros((2*options['radius']+options['blockSize']+1,2*options['radius']+options['blockSize']+1),dtype = 'float64')
            for centerRow in range(options['radius'], options['radius']+options['blockSize']): 
                for centerCol in range(options['radius'], options['radius']+options['blockSize']):
                    artifactBS1SumArray = addData(artifactBS1SumArray, artifactBS1Result, centerRow, centerCol, options)

            centerArtifactBS1SumArray = artifactBS1SumArray[(options['blockSize']-1)/2:2*options['radius']+(options['blockSize']-1)/2+1,(options['blockSize']-1)/2:2*options['radius']+(options['blockSize']-1)/2+1]
            # print centerArtifactBS1SumArray
            correctionArray = npy.where(artifactBSXResult>0,npy.divide(centerArtifactBS1SumArray,artifactBSXResult),0)
            # print correctionArray

            grid = npy.indices((correctionArray.shape))
            subsetRowArray = grid[0]
            subsetColArray = grid[1]
            del grid
            subsetDistArray = npy.sqrt(npy.multiply(subsetRowArray - subsetCenterRow, subsetRowArray- subsetCenterRow) + npy.multiply(subsetColArray-subsetCenterCol, subsetColArray-subsetCenterCol))           

            correctionArray = npy.where(subsetDistArray > 9*options['radius']/10,1,correctionArray)
            del subsetDistArray
            
            correctionFile = path.join(options['scratchDir'], 'correction.asc')
            ascii_grid_writer(correctionFile, correctionArray, artifactHeader, options)



        # if options['removeArtifacts']:
            # correctionArray = arcpy.RasterToNumPyArray(os.path.join(options['projectDir'],'div3.tif'),"#","#","#",-9999)
            # correctionArray = npy.where(correctionArray<0,0.0,correctionArray)

            # Set raster paths and export to ascii if needed. FIXME: don't really need ascii, just convenient for header code for now
        resisRaster = path.join(options['projectDir'],options['resisRasterBase'])
        descData = arcpy.Describe(resisRaster)
        spatialReference = descData.spatialReference
        options['cellSize'] = descData.MeanCellHeight
        
        if options['useClimate']:
            if options['matchClimatePCs']:
                t1PC1Raster = path.join(options['projectDir'],options['t1PC1RasterBase'])
                t1PC2Raster = path.join(options['projectDir'],options['t1PC2RasterBase'])
                t2PC1Raster = path.join(options['projectDir'],options['t2PC1RasterBase'])
                t2PC2Raster = path.join(options['projectDir'],options['t2PC2RasterBase']) 
            else:    
                climateRaster = path.join(options['projectDir'],options['climateRasterBase'])
        
        if options['useSourceRaster']:
            sourceRaster = path.join(options['projectDir'],options['sourceRasterBase'])
        elif options['sourceInverseR']:
            sourceRaster = 1 / arcpy.Raster(resisRaster)
        else:
            sourceRaster = arcpy.sa.Con((arcpy.Raster(resisRaster) <= options['rCutoff']),1,0)
        resisRaster = square_resistances(resisRaster,options)
        
        if options['squareRootInputs']:
            resisRaster, sourceRaster = sqr_root_inputs(resisRaster, sourceRaster, options)
            
        # Next line needed to avoid shift of one pixel in band function using rasterToNumpyArray. Don't ask me why.               
        sourceRaster = arcpy.sa.Con(arcpy.sa.Raster(resisRaster)>0,sourceRaster) 
        
        if options['calcCurrent']:
            cumCurrentRaster = arcpy.sa.Con((arcpy.Raster(resisRaster) > 0),0,0)
        else:
            cumCurrentRaster = None
        
        if options['calcFA']:    
            resisRasterFA = add_random_resistances(resisRaster,options)
        
        header = get_header(resisRaster)
          
        descData = arcpy.Describe(resisRaster)
        arcpy.env.extent = descData.Extent    
        
        if options['calcVoltages']:
            cumVdiffRaster = arcpy.sa.Con((arcpy.Raster(resisRaster) > 0),0,0)   
        else: cumVdiffArray = cumVdiffRaster = None       
        if options['calcFA']:
            cumFlowRaster = arcpy.sa.Con((arcpy.Raster(resisRaster) > 0),0,0)
        else: cumFlowRaster = None
        if options['saveSourceAndTargetCounts']:# and options['useClimate']: 
            cumTargetRaster = arcpy.sa.Con((arcpy.Raster(resisRaster) > 0),0,0)
            cumSourceRaster = arcpy.sa.Con((arcpy.Raster(resisRaster) > 0),0,0)
        else:
            cumTargetRaster = cumSourceRaster = None
        
        if options['endBand'] == 0:
            approxEndBand = str(int(header['nrows']/options['blockSize'])+1)
        else: approxEndBand = str(options['endBand'])
        if options['endStripe'] == 0:
            approxEndStripe = str(int(header['ncols']/options['blockSize'])+1)
        else: approxEndStripe = str(options['endStripe'])
        maxNumSolvesToDo = (int(approxEndBand)-options['startBand'] + 1)*(int(approxEndStripe)-options['startStripe']+1) + 1
        iter = 0
        bandNum = 0
        pctDone=0
        pctDone2=0
        for centerRow in range((options['blockSize']-1)/2,header['nrows'],options['blockSize']):
            solveInBand = False
            bandNum += 1 
            if options['startBand'] > 0 and bandNum < options['startBand']: continue
            if options['endBand'] > 0 and bandNum >= options['endBand']+1: break
            print 'Starting band #',str(bandNum)+'/'+approxEndBand,' centered on row '+str(centerRow)

            bandArray = band(resisRaster,header,centerRow, options)
            if options['calcFA']:
                bandArrayFA = band(resisRasterFA,header,centerRow, options)
            if npy.max(bandArray) == -9999: continue
            sourceBandArray = band(sourceRaster,header,centerRow, options)
            
            if npy.max(sourceBandArray) <= 0:       
                print 'No sources in band; continuing'
                continue

            if options['useClimate']:
                if options['matchClimatePCs']:
                    t1PC1BandArray = band(t1PC1Raster, header, centerRow, options)
                    t1PC2BandArray = band(t1PC2Raster, header, centerRow, options)
                    t2PC1BandArray = band(t2PC1Raster, header, centerRow, options)
                    t2PC2BandArray = band(t2PC2Raster, header, centerRow, options)
                else:
                    climateBandArray = band(climateRaster, header, centerRow, options)
                    
    #fixme: getting extra nodata row on top of climate band?
                #xprint 'climband'
                #xprint climateBandArray
            if options['calcCurrent']:        
                cumCurrentArray = npy.zeros(bandArray.shape, dtype = 'float64') 
            else:
                cumCurrentArray = None
            if options['saveSourceAndTargetCounts']:# and options['useClimate']:
                cumTargetArray = cumCurrentArray.copy() 
                cumSourceArray = cumCurrentArray.copy()
                
            if options['calcVoltages']: cumVdiffArray = cumCurrentArray.copy()
            subsetCenterRow = min(centerRow,options['radius'])
            
            # Check for all nodata in center band of resistance raster
            if options['blockSize']>1:
                bandCenterArray = bandArray[(subsetCenterRow-(options['blockSize']-1)/2):(subsetCenterRow+(options['blockSize']-1)/2),:]
            else:
                bandCenterArray = bandArray[subsetCenterRow,:]
            
             
            if npy.max(bandCenterArray) == -9999:       
                del bandCenterArray
                continue
            del bandCenterArray
            
            if options['blockSize']>1:
                sourceCenterArray = sourceBandArray[(subsetCenterRow-(options['blockSize']-1)/2):(subsetCenterRow+(options['blockSize']-1)/2+1),:]
            else:
                sourceCenterArray = sourceBandArray[subsetCenterRow,:]
                
            if npy.max(sourceCenterArray) <= 0:            
                print 'no targets; continuing' #fixme: need to differently handle sources and targets
                del sourceCenterArray
                continue
            # del sourceCenterArray fixme: saving for now to quickly look for valid target areas 
            
            sourceCenterArraySum0 = npy.sum(npy.where(sourceCenterArray > 0, sourceCenterArray, 0), axis = 0)
            
            grid = npy.indices((bandArray.shape))
            rowArray = grid[0]
            colArray = grid[1]
            del grid
            stripeNum = 0
            for centerCol in range((options['blockSize']-1)/2, header['ncols'],options['blockSize']):
                start_time0 = time.clock()
                stripeNum += 1
                if options['startStripe'] > 0 and stripeNum < options['startStripe']: continue
                if options['endStripe'] > 0 and stripeNum >= options['endStripe']+1: break
                iter += 1
                print 'Band #',bandNum,'Stripe #',stripeNum
    #fixme- check this:
                subsetCenterCol = min(centerCol,options['radius']) 
                subsetCenterRow = min(centerRow,options['radius'])
                # time consuming- add smart search. 
                if options['blockSize']>1:
                    if npy.max(sourceCenterArraySum0[centerCol-(options['blockSize']-1)/2:centerCol+(options['blockSize']-1)/2+1]) <= 0: 
                        print'No sources in radius'
                        continue
                else:
                    if sourceCenterArray[centerCol] <= 0: #BHM ADDED INDEX 10/12/15 
                        print'No sources in radius'
                        continue 

                # time consuming- add smart search. 
                circleResisArray = circ(bandArray, rowArray, colArray, subsetCenterRow, centerCol, options['radius']) #fixme: move grid/rowarray calc into circ fn? 
                if options['calcFA']:            
                    circleResisArrayFA = circ(bandArrayFA, rowArray, colArray, subsetCenterRow, centerCol, options['radius'])

                sourceArray = circ(sourceBandArray, rowArray, colArray, subsetCenterRow, centerCol, options['radius'])
                sourceArray = npy.where(sourceArray < 0, 0, sourceArray)    #fixme- just set nodata sources to zero earlier    

                targetArray = center_block(sourceArray, options['blockSize'], subsetCenterRow, subsetCenterCol)                             
                targetSum = targetArray.sum()
                if targetSum == 0:
                    continue
                if options['useDistanceFunction'] and (options['maxDist'] or options['minDist'] or options['distEq'] is not None):
                    sourceArray = modify_source_strengths_by_distance(sourceArray, subsetCenterRow, subsetCenterCol, options)          

                sourceArray = npy.where(targetArray > 0, 0, sourceArray)    
                sourceSum = sourceArray.sum()
                if sourceSum == 0:
                    continue
                
                if options['centerGround']:
                    groundArray = npy.where(targetArray > 0, -9999, -9999) # fixme- replaces with zeros-9999 to speed up?
                    groundArray[subsetCenterRow, subsetCenterCol] = 0
                else:
                    groundArray = npy.where(targetArray > 0, 10000000, -9999)
                if len(sourceArray) < 10:
                    print 'raw sources and targets:'
                    print sourceArray
                    print targetArray

                if options['useClimate']:
                    if options['matchClimatePCs']:               
                        sourceArray, targetArray, sourceSum, targetSum, targetTally = match_climate_pcs(sourceArray, targetArray, t1PC1BandArray, t1PC2BandArray, t2PC1BandArray, t2PC2BandArray, rowArray, colArray, subsetCenterRow, centerCol, options)                
                    else:
                        sourceArray, targetArray, sourceSum, targetSum, targetTally = match_temperature_diffs(sourceArray, targetArray, climateBandArray, rowArray, colArray, subsetCenterRow, centerCol, options)                
                if sourceSum == 0 or targetSum == 0:
                    continue 
                        
                circleHeader = get_subset_header(sourceArray, header, options, centerRow, centerCol)
                yMin = circleHeader['yllcorner'] #fixme- check. coudl be something like: max(circleHeader['yllcorner'],circleHeader['yllcorner'] + ((circleHeader['nrows'] - centerRow - options['radius'] - 1) * circleHeader['cellsize']))
                if options['saveSourceAndTargetCounts']:# and options['useClimate']:
                    cumSourceArray = addData(cumSourceArray, sourceArray, subsetCenterRow, centerCol, options)
                    cumTargetArray = addData(cumTargetArray, targetArray, subsetCenterRow, centerCol, options)
                    if options['doSourceAndTargetCountsOnly']:
                        continue
                #then normalizing. 1 amp injected, 1 amp taken out
                sourceArray = sourceArray/(sourceSum+0.0)
                
                if options['negTargets']:
                    targetArray = -targetArray/(targetSum+0.0) 
                    sourceArray += targetArray
                print 'Source and target sums:',sourceSum,targetSum
                del targetArray
                
                if options['calcCurrent']:
                    sourceAsciiFile = path.join(options['scratchDir'], 'source_r'+str(centerRow) + 'c' +str(centerCol)+'iter'+str(iter)+'.asc') 
                    
                
                groundAsciiFile = path.join(options['scratchDir'], 'ground_r'+str(centerRow) + 'c' +str(centerCol)+'iter'+str(iter)+'.asc')            
                ascii_grid_writer(groundAsciiFile, groundArray, circleHeader, options)
                # save_numpy(groundArray) #FIXME: shows that could save some time, not a huge amount.
                del groundArray 
                
                outputFN = 'target_r'+str(centerRow) + 'c' +str(centerCol)+'.out'
                outputFile = path.join(options['scratchDir'], outputFN)
                resisAsciiFile = path.join(options['scratchDir'], 'resis_r'+str(centerRow) + 'c' +str(centerCol)+'iter'+str(iter)+'.asc') 

                if options['calcCurrent']:
                    csOptions = setCircuitscapeOptions(resisAsciiFile,sourceAsciiFile,groundAsciiFile,outputFile)

                if options['calcVoltages']:
                    csOptions['write_volt_maps'] = True
                
                if options['calcNull']:
                    circleResisArray = npy.where(circleResisArray>0,1,circleResisArray)
                    if options['calcFA'] and  options['calcFANull']:
                        circleResisArrayFA = npy.where(circleResisArrayFA>0,1,circleResisArrayFA)

                print '\nDone with prep'
                start_time0 = elapsed_time(start_time0)           

                # print 'Solving stripe #' + str(stripeNum) + '/' + approxEndStripe,' and band #' + str(bandNum) + '/' + approxEndBand+'\n'

                # maxSolvesDone=(bandNum-options['startBand']+1)*numstripes + (stripeNum-options['startStripe']+1)
                pctDone = report_pct_done(iter, maxNumSolvesToDo, -1)

                print 'Elapsed time so far: {0}'.format(datetime.datetime.now()-theStart)  
                
                solveInBand = True
                if options['calcFA']:
                    rFA2, options, flowSourceCorrection, maxFlow, sourceArray = calc_fa(options,groundAsciiFile,circleHeader,yMin,sourceArray,circleResisArrayFA,iter,resisRaster)
                    
                    if maxFlow <= 0:
                        print 'NO FLOW, continuing'
                        continue
                        
                if options['calcCurrent']:
                    if options['limitCalcsByCWD']:          
                        LLC = arcpy.Point(circleHeader['xllcorner'],circleHeader['yllcorner'])
                        # cwdLimit = (options['cwdLimit'] * options['cellSize'])
                        cwdArray = arcpy.RasterToNumPyArray(options['cwdRaster'],LLC,"#","#",-9999)  
                        
                        # Fixme: put in fn
                        if cwdArray.shape[1] > circleHeader['ncols']: #FIXME: temp fix because random raster yields bands that are 1 too long for some reason      
                            print 'Array is too long'
                            print 'ncols',circleHeader['ncols']
                            print 'shape',cwdArray.shape[1]
                            print 'last col min,max', npy.min(cwdArray[:,circleHeader['ncols']]), npy.max(cwdArray[:,circleHeader['ncols']])
                            print 'Col 0 min,max', npy.min(cwdArray[:,0]), npy.max(cwdArray[:,0])
                            cwdArray = cwdArray[:,1:circleHeader['ncols']]                       
    #fixme: should sources or resistances be changed based on cwd?
                        sourceArray =npy.where((cwdArray < 0), 0, sourceArray)
                        # circleResisArray = npy.where((cwdArray < 0), -9999, circleResisArray)
                        del cwdArray
                    ascii_grid_writer(sourceAsciiFile, sourceArray, circleHeader, options)
                    ascii_grid_writer(resisAsciiFile, circleResisArray, circleHeader, options)
                    curMapPath, outConfigFile = calc_current(options, csOptions, centerRow, centerCol) 
                    if not path.exists(curMapPath):
                        print "Can't find circuitscape output"
                        raw_input('\nPress Enter to continue.') 
                        exit(0)
                        continue
                    currentArray = ascii_grid_reader(curMapPath, header['nodata'], 'float64')
                    if options['centerGround']:
                        maxCur = currentArray[subsetCenterRow, subsetCenterCol]
                    else:
                        maxCur = npy.max(currentArray)
                    if maxCur <= 0:
                        print 'NO CURRENT, continuing'
                        continue
                    
                    currentSourceCorrection = 1
                    if (options['weightTargOnly']) or (options['noWeight']):
                        if 1-maxCur > 0.01:
                            currentSourceCorrection = maxCur

                    if options['fadeIn']:
                        currentArray = fade_currents(currentArray, subsetCenterRow, subsetCenterCol, options)  
                    elif options['removeArtifacts']:
                        correctionArray2 = trim(correctionArray, centerRow, centerCol, header['nrows'], header['ncols'], options)
                        currentArray = correctCurrents(currentArray, subsetCenterRow, subsetCenterCol, options, correctionArray2)
                        
                if options['weightTargOnly']:
                    if options['useClimate']:
                        multiplier = targetTally 
                    else:
                        multiplier = targetSum
                elif options['useClimate']:
                    multiplier = targetSum
                elif options['noWeight']:
                    multiplier = 1
                else:    
                    multiplier = sourceSum * targetSum #fixme: combos of large radii and large blocks can exceed the ~2.1 billion limit for int32. Even going to int64 seems to crash arc.sa.Times.
                
                del sourceArray
                if options['calcFA']:                     
                    flowMultiplier = multiplier/flowSourceCorrection
                    #fixme: next line uses int, does this negatively affect any cases?
                    rFA3 = arcpy.sa.Times(rFA2, int(flowMultiplier)) #fixme: may be faster to convert rfa or rfa2 into array and do numpy calcs for entire band, then add band in
                    cumFlowRaster = addData_arcpy(rFA3, cumFlowRaster)
                    if options['limitCalcsByCWD']:
                        delete_data(options['cwdRaster'])

                        
                        
    #fixme- may not be compatible with subsources, weights, polygons, etc:
                
    #put in fn once striping solved below
                if options['calcCurrent']:
                    currentMultiplier = multiplier/currentSourceCorrection  
                    currentArray = currentArray * currentMultiplier    
                    delete_data(outConfigFile)
                    delete_data(resisAsciiFile)
                    delete_data(sourceAsciiFile)
                    cumCurrentArray = addData(cumCurrentArray, currentArray, subsetCenterRow, centerCol, options)
                    del currentArray


                if options['calcVoltages']:

                    voltMapPath = path.join(options['scratchDir'], 'target_r'+str(centerRow) + 'c' +str(centerCol) + '_voltmap.asc')
                    voltMap = (ascii_grid_reader(voltMapPath, header['nodata'], 'float64'))
                    voltMap = npy.where(circleResisArray<0,npy.nan,voltMap)
                    #xprint 'voltmap'
                    #xprint voltMap
                    
                    max_vdiff = get_max_vdiff(voltMap,circleResisArray,csOptions)
                    max_vdiff = npy.where(npy.isnan(max_vdiff),0,max_vdiff)

    # fixme temp fix for striping at edges
                    max_vdiff[0,::] = 0
                    max_vdiff[:,0] = 0
                    max_vdiff[:,circleHeader['ncols']-1] = 0
                    max_vdiff[circleHeader['nrows']-1,:] = 0
                    # vDiffPath = path.join(options['scratchDir'], 'target_r'+str(centerRow) + 'c' +str(centerCol) + '_VDIFF.asc')
                    # ascii_grid_writer(vDiffPath, max_vdiff, circleHeader, options)

                    if options['fadeIn']: #Fixme: do we want to fade voltages?
                        max_vdiff = fade_currents(max_vdiff, subsetCenterRow, subsetCenterCol, options)  
                    elif options['removeArtifacts']:
                        correctionArray2 = trim(correctionArray, centerRow, centerCol, header['nrows'], header['ncols'], options)
                        max_vdiff = correctCurrents(max_vdiff, subsetCenterRow, subsetCenterCol, options, correctionArray2)

                        # npy.where(max_vdiff>0,npy.multiply(max_vdiff,fadeArray),max_vdiff)
                        # del fadeArray

                    max_vdiff = currentMultiplier * max_vdiff       
                    # if options['weightTargOnly']:
                        # # fixme: need to adjust multiplier for currentAdjustment?
                        # max_vdiff = targetSum * max_vdiff
                    # elif options['noWeight']:
                        # pass
                    # else:    
                        # max_vdiff = sourceSum * targetSum * max_vdiff
                    cumVdiffArray = addData(cumVdiffArray, max_vdiff, subsetCenterRow, centerCol, options)
                    delete_data(voltMapPath)
                delete_data(groundAsciiFile)
                del circleResisArray               
                start_time1 = time.clock()

    # FIXME: do all raster stuff for each band, not each solve
                # cumCurrentRaster = addData_arcpy(cumCurrentRaster, currentRaster)
    # FIXME: move file writing to band iteration

            if options['calcCurrent']:
                del rowArray,colArray
                del sourceCenterArray,bandArray,sourceBandArray
            if options['calcFA']:
                del bandArrayFA
                
            if options['saveSourceAndTargetCounts']:
                yMin = max(header['yllcorner'],header['yllcorner'] + ((header['nrows'] - centerRow - options['radius'] - 1) * header['cellsize']))
                LLC = arcpy.Point(header['xllcorner'],yMin)
                bandSourceRaster = arcpy.NumPyArrayToRaster(cumSourceArray,LLC, header['cellsize'],header['cellsize'],-9999)
                cumSourceRaster = addData_arcpy(cumSourceRaster, bandSourceRaster)       
            
                del cumSourceArray, bandSourceRaster
                bandTargetRaster = arcpy.NumPyArrayToRaster(cumTargetArray,LLC, header['cellsize'],header['cellsize'],-9999)
                cumTargetRaster = addData_arcpy(cumTargetRaster, bandTargetRaster)       
                del cumTargetArray, bandTargetRaster
            
            
            print 'Done with band#',str(bandNum)+'/'+approxEndBand
            # pctDone = report_pct_done((bandNum-options['startBand']+1)*(stripeNum-options['startStripe']+1), maxNumSolvesToDo, -1)

            pctDone = report_pct_done(iter+1, maxNumSolvesToDo, -1)
            if solveInBand: #
                if options['calcCurrent']:
                    yMin = max(header['yllcorner'],header['yllcorner'] + ((header['nrows'] - centerRow - options['radius'] - 1) * header['cellsize']))
                    LLC = arcpy.Point(header['xllcorner'],yMin)
                    bandCurrentRaster = arcpy.NumPyArrayToRaster(cumCurrentArray,LLC, header['cellsize'],header['cellsize'],-9999)                          

                    # SAVING BAND RASTER REMOVES OCCASIONAL HORIZONTAL STRIPING
                    tempBandFile = os.path.join(options['scratchDir'], 'justBAND'+str(bandNum)+'cur_' + options['outputFileText']+'.tif')
                    bandCurrentRaster.save(tempBandFile)
                    delete_data(tempBandFile)
                    
                    cumCurrentRaster = addData_arcpy(cumCurrentRaster, bandCurrentRaster)
                    del bandCurrentRaster,cumCurrentArray
                    if options['calcVoltages']:
                        bandVdiffRaster = arcpy.NumPyArrayToRaster(cumVdiffArray,LLC,header['cellsize'],header['cellsize'],-9999)
                        cumVdiffRaster = addData_arcpy(cumVdiffRaster, bandVdiffRaster)
                        del cumVdiffArray, bandVdiffRaster
                  
                options = write_temp_maps(options,bandNum,cumCurrentRaster,cumVdiffRaster,cumFlowRaster)
                                         
            # print 'Done with band #',bandNum,',
            print 'Elapsed time so far: {0}'.format(datetime.datetime.now()-theStart)  
        print 'Done with solves.'  
        print_prof_data()
        write_final_maps(options,cumCurrentRaster,cumVdiffRaster,cumFlowRaster,cumSourceRaster,cumTargetRaster) 
        #xprint locals()
        if options['calcFA']: delete_data(resisRasterFA)
        if options['squareResistances']: delete_data(resisRaster)
        return options
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()
    
def profile(fn):
    @wraps(fn)
    def with_profiling(*args, **kwargs):
        start_time = time.time()

        ret = fn(*args, **kwargs)

        elapsed_time = time.time() - start_time

        if fn.__name__ not in PROF_DATA:
            PROF_DATA[fn.__name__] = [0, []]
        PROF_DATA[fn.__name__][0] += 1
        PROF_DATA[fn.__name__][1].append(elapsed_time)

        return ret

    return with_profiling

def print_prof_data():
    try:
        if not options['printTimings']:
            return
    except:
        return
    print "*** TIMING DATA ****"
    for fname, data in PROF_DATA.items():
        max_time = max(data[1])
        avg_time = sum(data[1]) / len(data[1])
        
        print "Function %s called %d times. " % (fname, data[0]),
        print 'Execution time max: %.3f, average: %.3f' % (max_time, avg_time),
        print '. Total time = ' + str(data[0]*avg_time) 
    print '\n'    
        
def clear_prof_data():
    global PROF_DATA
    PROF_DATA = {}

    
def get_null_result(artifactArray, artifactHeader, centerRow, centerCol, blockSize):
    try:
        grid = npy.indices((artifactArray.shape))
        rowArray = grid[0]
        colArray = grid[1]
        del grid
        sourceArray = circ(artifactArray, rowArray, colArray, centerRow, centerRow, options['radius'])
           
        resisArray = sourceArray.copy()
        sourceArray = npy.where(sourceArray < 0, 0, sourceArray) 
        targetArray = center_block(sourceArray, blockSize, centerRow, centerCol)                             
        targetSum = targetArray.sum()
        if options['useDistanceFunction'] and (options['maxDist'] or options['minDist'] or options['distEq'] is not None):
            sourceArray = modify_source_strengths_by_distance(sourceArray, centerRow, centerCol, options)          

        sourceArray = npy.where(targetArray > 0, 0, sourceArray)    
        sourceSum = sourceArray.sum()
        if options['centerGround']:
            groundArray = npy.where(targetArray > 0, -9999, -9999) # fixme- replaces with zeros-9999 to speed up?
            groundArray[centerRow, centerCol] = 0
        else:
            groundArray = npy.where(targetArray > 0, 10000000, -9999)

        sourceAsciiFile = path.join(options['scratchDir'], 'source_artifact_BS'+str(blockSize)+'.asc') 
        resisAsciiFile = path.join(options['scratchDir'], 'resis_artifact_BS'+str(blockSize)+'.asc') 
        groundAsciiFile = path.join(options['scratchDir'], 'ground_artifact_BS'+str(blockSize)+'.asc') 

        ascii_grid_writer(groundAsciiFile, groundArray, artifactHeader, options)
        ascii_grid_writer(resisAsciiFile, resisArray, artifactHeader, options)
        ascii_grid_writer(sourceAsciiFile, sourceArray, artifactHeader, options)
        outputFN = 'target_r'+str(centerRow) + 'c' +str(centerCol) + '.out'#_curmap.asc'

        outputFile = path.join(options['scratchDir'], outputFN)
        csOptions = setCircuitscapeOptions(resisAsciiFile,sourceAsciiFile,groundAsciiFile,outputFile)
        curMapPath, outConfigFile = calc_current(options, csOptions, centerRow, centerCol) 
        if not path.exists(curMapPath):
            print "Can't find circuitscape output"
            raw_input('\nPress Enter to continue.') 
            exit(0)
            
        currentArray = ascii_grid_reader(curMapPath, artifactHeader['nodata'], 'float64')

        if options['weightTargOnly']:
            multiplier = targetSum
        elif options['noWeight']:
            multiplier = 1
        else:    
            multiplier = sourceSum * targetSum  
        print sourceSum, targetSum, multiplier
        currentArray = currentArray * multiplier
        
        return currentArray
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()
        
def createArtifactHeader(options): 
    nullHeader = {}
    nullHeader['ncols'] = 2*options['radius']+1
    nullHeader['nrows'] = 2*options['radius']+1
    nullHeader['xllcorner'] = 1
    nullHeader['yllcorner'] = 1
    nullHeader['cellsize'] = 1
    nullHeader['nodata'] = -9999 #if (nodata == False) else nodata 
    return nullHeader
    
    
def report_pct_done(current, goal, last):
    try:
        """Reports percent done"""
        goal = float(goal)
        pctDone = ((float(current) / goal) * 100)
        # pctDone = 10 * (npy.floor(pctDone/10))
        pctDone2 = (float(npy.floor(pctDone*10))/10)
        # print pctDone2,'CT'
        # if pctDone - last >= 1 or pctDone > 0.5:
            # print str(int(pctDone))+  "% done.",
            # # return 10*int((npy.floor(pctDone/10)))
            # return int(npy.floor(pctDone))
        # elif pctDone < 0.5:
        print "~" + str(pctDone2)+  "% done.",
        # return 10*int((npy.floor(pctDone/10)))
        return int(npy.floor(pctDone))
        # return last
    except:    
        print_python_error()

# To evaluate speed differences. definitely faster than saving and reading asciis. may speed up cs too?    
# @profile
# def save_numpy(array):
    # npyFile = r'c:\temp\temp.npy'
    # npy.save(npyFile, array)     
    # npyArray = npy.load(npyFile, mmap_mode = None)

    
@profile
def calc_current(options, csOptions, centerRow, centerCol): 
    try:
        configFN = 'config_target_b'+str(centerRow) + 'c' +str(centerCol)+'.ini'
        outConfigFile = path.join(options['scratchDir'], configFN)
        writeCircuitscapeConfigFile(outConfigFile, csOptions)  
        solveInBand = True
        CSPATH = get_cs_path()                            
        start_time1 = time.clock()
        call_circuitscape(CSPATH, outConfigFile)
        
        start_time1 = elapsed_time(start_time1)

        delete_data(configFN)
        configFN2 = 'target_r'+str(centerRow) + 'c' +str(centerCol)+'.ini' #fixme?
        delete_data(path.join(options['scratchDir'], configFN2))

        curMapPath = path.join(options['scratchDir'], 'target_r'+str(centerRow) + 'c' +str(centerCol) + '_curmap.asc')
        # Adjust currents if some sources in different components than target block
        return curMapPath, outConfigFile
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

@profile
def calc_fa(options,groundAsciiFile,circleHeader,yMin,sourceArray,circleResisArrayFA, iter, resisRaster):
    try:
        oldExtent = arcpy.env.extent
        arcpy.env.extent = "MAXOF"
        # fixme: convoluted way to get a shapefile via raster conversion. simplify
        outPoint = path.join(options['scratchDir'], 'point_iter'+str(iter)+'.shp')
        if arcpy.Exists(outPoint):
            arcpy.Delete_management(outPoint)
        field = "VALUE"
        
        arcpy.RasterToPoint_conversion(groundAsciiFile, outPoint, field)    
        LLC = arcpy.Point(circleHeader['xllcorner'],yMin)

        resisRasFA = arcpy.NumPyArrayToRaster(circleResisArrayFA,LLC, circleHeader['cellsize'],circleHeader['cellsize'],-9999)    
        del circleResisArrayFA
        
        print '\nStarting flow accumulation'
        start_timeFA = time.clock()
        if options['limitCalcsByCWD']:
            options['cwdRaster'] = path.join(options['scratchDir'],'cwd'+str(iter)+'.tif')
            cwdLimit = options['cwdLimit'] * options['cellSize']
            rCB = arcpy.sa.CostBackLink(outPoint, resisRasFA, cwdLimit, options['cwdRaster'])#, rOut + "d" + str(FID) ) #BHM can limit cwd here 999999999
            rFD = arcpy.sa.Con( rCB > 0, arcpy.sa.Power(2, (rCB - 1))) 
            rFD = arcpy.sa.Con( rCB == 0, 11, rFD) #Bogus direction ensures that flow at target point will be accumulative

            # outCon = arcpy.sa.Con(arcpy.sa.Raster(options['cwdRaster'])<cwdLimit,sourceRasFA,0)       
            if options['useCwDistanceFunction']:
                print 'CWDDISTANCE FUNCTION NOT TESTED YET'
                sourceArray = modify_source_strengths_by_cwd(sourceArray, circleHeader, LLC, yMin, options) #do this once, pass back
            sourceRasFA = arcpy.NumPyArrayToRaster(sourceArray,LLC, circleHeader['cellsize'],circleHeader['cellsize'],-9999)              
            rFA = arcpy.sa.FlowAccumulation( rFD, sourceRasFA ) #Note! removed source strengths again 

        else:
            sourceRasFA = arcpy.NumPyArrayToRaster(sourceArray,LLC, circleHeader['cellsize'],circleHeader['cellsize'],-9999)              
            rCB = arcpy.sa.CostBackLink(outPoint, resisRasFA, '#')
            rFD = arcpy.sa.Con( rCB > 0, arcpy.sa.Power(2, (rCB - 1))) 
            rFD = arcpy.sa.Con( rCB == 0, 11, rFD) #Bogus direction ensures that flow at target point will be accumulative
            rFA = arcpy.sa.FlowAccumulation( rFD, sourceRasFA ) 
    #temp
        # sourceRasFA.save('c:\\temp\\sources.tif')
        rFA2 = arcpy.sa.Con(arcpy.sa.IsNull(rFA),0,rFA)
        start_timeFA = elapsed_time(start_timeFA)
        arcpy.env.extent = oldExtent 
        arcpy.Delete_management(outPoint)
        del rCB, rFD, rFA
        # Adjust flow if some sources in different components than target block
        flowSourceCorrection = 1 
        maxFlow = get_raster_max(rFA2)
        if (options['weightTargOnly']) or (options['noWeight']):
            #Need to know if some sources were in different components than target block
            # if noweight, divide by sourceCorrection
            # if targonly, ""
            # else, do nothing (comes out in wash)  
            if 1-maxFlow > 0.01 and maxFlow > 0:
                flowSourceCorrection = maxFlow
        return rFA2, options, flowSourceCorrection, maxFlow, sourceArray
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()   

@profile
def get_raster_max(raster):
    try:
        maxObject = arcpy.GetRasterProperties_management(raster, "MAXIMUM") 
        max = float(str(maxObject.getOutput(0)))
        return max
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

    
def correctCurrents(currentArray, subsetCenterRow, subsetCenterCol, options, correctionArray):
    try:   
        #with just block:
        # currentArrayBlock = currentArray[(subsetCenterRow-(options['blockSize']-1)/2):(subsetCenterRow+(options['blockSize']-1)/2)+1,(subsetCenterCol-(options['blockSize']-1)/2):(subsetCenterCol+(options['blockSize']-1)/2+1)]
        # currentArrayBlock = npy.multiply(currentArrayBlock,correctionArray)
        # currentArray[(subsetCenterRow-(options['blockSize']-1)/2):(subsetCenterRow+(options['blockSize']-1)/2)+1,(subsetCenterCol-(options['blockSize']-1)/2):(subsetCenterCol+(options['blockSize']-1)/2+1)] = currentArrayBlock
        # return currentArray 
        # currentArrayBlock = currentArray[(subsetCenterRow-(options['blockSize']-1)/2):(subsetCenterRow+(options['blockSize']-1)/2)+1,(subsetCenterCol-(options['blockSize']-1)/2):(subsetCenterCol+(options['blockSize']-1)/2+1)]
        
        
        
        currentArray = npy.multiply(currentArray,correctionArray)
        # currentArray[(subsetCenterRow-(options['blockSize']-1)/2):(subsetCenterRow+(options['blockSize']-1)/2)+1,(subsetCenterCol-(options['blockSize']-1)/2):(subsetCenterCol+(options['blockSize']-1)/2+1)] = currentArrayBlock
        return currentArray 


    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()    
    
@profile
def fade_currents(currentArray, subsetCenterRow, subsetCenterCol, options):
    try:
        grid = npy.indices((currentArray.shape))
        subsetRowArray = grid[0]
        subsetColArray = grid[1]
        del grid
        subsetDistArray = npy.sqrt(npy.multiply(subsetRowArray - subsetCenterRow, subsetRowArray- subsetCenterRow) + npy.multiply(subsetColArray-subsetCenterCol, subsetColArray-subsetCenterCol))           
        fadeArray = npy.where(subsetDistArray<options['fadeInDist'],(subsetDistArray+options['fadeConstant'])/(options['fadeInDist']+options['fadeConstant']),1)

        # fadeArray = npy.where(subsetDistArray<options['fadeInDist'],0,1)
        currentArray = npy.where(currentArray>0,npy.multiply(currentArray,fadeArray),currentArray)

        # fixme: need to check for Nodata in kernel?
        if options['setKernelMean'] and options['blockSize'] > 1:
            if options['blockSize'] > 3:
                center16Currents = currentArray[subsetCenterRow-2:subsetCenterRow+3,subsetCenterCol-2:subsetCenterCol+3]
                currentArray[subsetCenterRow-2:subsetCenterRow+3,subsetCenterCol-2:subsetCenterCol+3]=npy.mean(center16Currents)
            else:
                center9Currents = currentArray[subsetCenterRow-1:subsetCenterRow+2,subsetCenterCol-1:subsetCenterCol+2]
                currentArray[subsetCenterRow-1:subsetCenterRow+2,subsetCenterCol-1:subsetCenterCol+2]=npy.mean(center9Currents)
        
        return currentArray 
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()
    
    
def square_resistances(inRaster,options):
    try:
        if options['squareResistances']: 
            descData = arcpy.Describe(inRaster)
            arcpy.env.extent = descData.Extent    

            fileBase, fileExtension = path.splitext(options['resisRasterBase'])
            resisRasterSQ = path.join(options['scratchDir'],'sq_'+fileBase+'.tif')
            outResistanceRaster = arcpy.sa.Times(inRaster, inRaster) 
            outResistanceRaster.save(resisRasterSQ)           
            inRaster = resisRasterSQ    
        return inRaster
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

        
def sqr_root_inputs(resisRaster, sourceRaster,options):
    try:
        descData = arcpy.Describe(resisRaster)
        arcpy.env.extent = descData.Extent    

        fileBase, fileExtension = path.splitext(options['resisRasterBase'])
        resisRasterSQRT = path.join(options['scratchDir'],'sqrt_'+fileBase+'.tif')
        outResistanceRaster = arcpy.sa.SquareRoot(resisRaster) 
        outResistanceRaster.save(resisRasterSQRT)           
        resisRaster = resisRasterSQRT    
            
        
        fileBase, fileExtension = path.splitext(options['sourceRasterBase'])
        sourceRasterSQRT = path.join(options['scratchDir'],'sqrt_'+fileBase+'.tif')
        outSourceRaster = arcpy.sa.SquareRoot(sourceRaster) 
        outSourceRaster.save(sourceRasterSQRT)           
        sourceRaster = sourceRasterSQRT            
        return resisRaster, sourceRaster
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()
        
            
@profile        
def match_climate_pcs(sourceArray, targetArray, t1PC1BandArray, t1PC2BandArray, t2PC1BandArray, t2PC2BandArray, rowArray, colArray, subsetCenterRow, centerCol, options):                
    try:
        print 'Starting climate using PCs'
        start_timeClimate = time.clock()
        t1PC1Array = circ(t1PC1BandArray, rowArray, colArray, subsetCenterRow, centerCol, options['radius'])              
        t1PC2Array = circ(t1PC2BandArray, rowArray, colArray, subsetCenterRow, centerCol, options['radius'])              
        t2PC1Array = circ(t2PC1BandArray, rowArray, colArray, subsetCenterRow, centerCol, options['radius'])              
        t2PC2Array = circ(t2PC2BandArray, rowArray, colArray, subsetCenterRow, centerCol, options['radius'])              
        sourceArrayNew = npy.zeros(sourceArray.shape,dtype = 'float64')
        # targetArrayNew = npy.zeros(sourceArray.shape,dtype = 'float64')
        
        # indices of valid targets
        tRows, tCols = npy.where(targetArray)
        targetTally = 0

        for i in range(0,len(tRows)):
            t2PC1Target = t2PC1Array[tRows[i],tCols[i]]     #Targets are future. Fixme: would need to calculate t1pcs for targets if doing absclim analogy
            t2PC2Target = t2PC2Array[tRows[i],tCols[i]] 
            if t2PC1Target == -9999 or t2PC2Target == -9999: #fixme: do nans
                targetArray[tRows[i],tCols[i]] = 0
                continue
            PC1DistArray = t1PC1Array - t2PC1Target #fixme: do nans here
            PC2DistArray = t1PC2Array - t2PC2Target
            if len(sourceArray) < 10:
                print 't2PC1Target,t2PC2Target'
                print t2PC1Target,t2PC2Target
                print 'PC1DistArray'
                print PC1DistArray
                print 'PC2DistArray'
                print PC2DistArray
                
            PCDistArray = npy.sqrt((PC1DistArray*PC1DistArray)+(PC2DistArray*PC2DistArray))
            PCDistArray = npy.where((t1PC1Array > -9999) & (t1PC2Array > -9999), PCDistArray, -9999) #fixme do nans or some other way to speed up. Only checking t1 here, assuming identical nodata areas
            if len(sourceArray) < 10:
                print 'PCDistArray '
                print PCDistArray 
            PCCutoffArray = npy.where((PCDistArray>= -options['PCWindow']) & (PCDistArray<= options['PCWindow']),1,0) 
            if len(sourceArray) < 10:
                print 'PCCutoffArray'
                print PCCutoffArray
            sourceArrayTarget_i = npy.multiply(sourceArray,PCCutoffArray)# sources for this target cell
            if npy.max(sourceArrayTarget_i)>0:
                targetTally+= 1 # There's a valid source for this target. Increment target count for scaling current in options['weightTargOnly'] setting
            
            # For climate, may have more or fewer targets. Need to scale sources by ntargs, and targets by nsources
            sourceArrayNew = sourceArrayNew+sourceArrayTarget_i
            
    #?/ not sure how to handle                                        
    #This is where sourcesum and targetsum become the same                    
    # putting in current equal to nsource*ntarg, taking same out. then normalizing to unit current. then mult by nsource*ntarg.
            #taking same out. 
            targetArray[tRows[i],tCols[i]] = sourceArrayTarget_i.sum()                   
        del tRows, tCols
               
        # take sum of sourceArrayTarget, set target strength to that...
        if len(sourceArray) < 10:
            print 'source new'
            print sourceArrayNew
            print 'target new'
            print targetArray
            
        sourceArray = sourceArrayNew
        del sourceArrayNew
        sourceSum = sourceArray.sum()
        targetSum = targetArray.sum() 
        #xprint 'targetTally',targetTally
        # Note that we now have identical source sums and targetsums
        start_timeClimate = elapsed_time(start_timeClimate) #Fixme: can climate be sped up? Taking up to 1 second with blocksize 25, radius 100           

        return sourceArray, targetArray, sourceSum, targetSum, targetTally
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

@profile
def match_temperature_diffs(sourceArray, targetArray, climateBandArray, rowArray, colArray, subsetCenterRow, centerCol, options):                
    try:
        # Fixme: need to add a quick search for overall potential for matches within bands
        print 'Starting climate using temperature differential'
        start_timeClimate = time.clock()
        climateArray = circ(climateBandArray, rowArray, colArray, subsetCenterRow, centerCol, options['radius'])              
        sourceArrayNew = npy.zeros(sourceArray.shape,dtype = 'float64')
        
        # indices of valid targets
        tRows, tCols = npy.where(targetArray)
        if len(sourceArray) < 10:
            print 'clim'
            print climateArray
        targetTally = 0
        for i in range(0,len(tRows)):
            tTarget = climateArray[tRows[i],tCols[i]] 
            if tTarget == -9999:
                targetArray[tRows[i],tCols[i]] = 0
                continue
            tDiffArray = npy.where(climateArray == -9999, 0, climateArray-tTarget) #target is cooler. # fixme- could save a bit of time by removing where- maybe could just mask outinvalid- or uses nans, or remove -999 values later
           
            if options['absClim']:
                tCutoffArray = npy.where((abs(tDiffArray)>= options['tDiff']-options['tWindow']) & (abs(tDiffArray)<= options['tDiff']+options['tWindow']),1,0)
            else:
                tCutoffArray = npy.where((tDiffArray>= options['tDiff']-options['tWindow']) & (tDiffArray<= options['tDiff']+options['tWindow']),1,0) 
            if len(sourceArray) < 10:
                print 'tCutoffArray'
                print tCutoffArray
            sourceArrayTarget_i = npy.multiply(sourceArray,tCutoffArray)# sources for this target cell
            if npy.max(sourceArrayTarget_i)>0:
                targetTally+= 1 # There's a valid source for this target. Increment target count for scaling current in options['weightTargOnly'] setting
            
            # For climate, may have more or fewer targets. Need to scale sources by ntargs, and targets by nsources
            sourceArrayNew = sourceArrayNew+sourceArrayTarget_i
            
    #?/ not sure how to handle                                        
    #This is where sourcesum and targetsum become the same                    
    # putting in current equal to nsource*ntarg, taking same out. then normalizing to unit current. then mult by nsource*ntarg.
            #taking same out. 
            targetArray[tRows[i],tCols[i]] = sourceArrayTarget_i.sum()                   
        del tRows, tCols
            
        # take sum of sourceArrayTarget, set target strength to that...
        if len(sourceArray) < 10:
            print 'source new'
            print sourceArrayNew
            print 'target new'
            print targetArray
        sourceArray = sourceArrayNew
        del sourceArrayNew
        sourceSum = sourceArray.sum()
        targetSum = targetArray.sum()
        print 'targetTally',targetTally
        # Note that we now have identical source sums and targetsums
        start_timeClimate = elapsed_time(start_timeClimate) #Fixme: can climate be sped up? Taking up to 1 second with blocksize 25, radius 100           
        return sourceArray, targetArray, sourceSum, targetSum, targetTally
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()
    
@profile
def modify_source_strengths_by_cwd(sourceArray, circleHeader, LLC, yMin, options): 
    try:
        if options['cwDistEq'] is None: return sourceArray
        
        cwdArray = arcpy.RasterToNumPyArray(options['cwdRaster'],LLC,"#","#",-9999)  
        # Fixme: put in fn
        if cwdArray.shape[1] > circleHeader['ncols']: #FIXME: temp fix because random raster yields bands that are 1 too long for some reason      
            print 'Array is too long'
            print 'ncols',circleHeader['ncols']
            print 'shape',cwdArray.shape[1]
            print 'last col min,max', npy.min(cwdArray[:,circleHeader['ncols']]), npy.max(cwdArray[:,circleHeader['ncols']])
            print 'Col 0 min,max', npy.min(cwdArray[:,0]), npy.max(cwdArray[:,0])
            cwdArray = cwdArray[:,1:circleHeader['ncols']]                       

        cwDist = cwdArray
        cwdLimit = options['cwdLimit'] * options['cellSize']
        cwDistanceFunctionArray = npy.where(cwdArray > 0, eval(options['cwDistEq']), 0)
        if len(sourceArray) < 10:
            print 'cwDist'
            print cwDist
            print'cw distance Function'
            print cwDistanceFunctionArray
        del cwDist

        sourceArrayNew = npy.where(sourceArray>0,npy.multiply(sourceArray,cwDistanceFunctionArray),sourceArray)
        # print npy.max(sourceArrayNew-sourceArray)
        # print npy.min(sourceArrayNew-sourceArray)
        if len(sourceArray) < 10:
            print 'source, dist, distfn, distmodifier,newsource:'
            print 'sourceArray\n',sourceArray
            print 'cwdArray\n',cwdArray
            print 'distanceFunctionArray\n',cwDistanceFunctionArray
            print 'sourceArrayNew\n',sourceArrayNew        
       
        return sourceArrayNew
        # del subsetRowArray, sourceArrayNew, subsetColArray, subsetDistArray, distanceModifierArray, distanceFunctionArray
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

        

@profile
def modify_source_strengths_by_distance(sourceArray, subsetCenterRow, subsetCenterCol, options):
    try:
        grid = npy.indices((sourceArray.shape))
        subsetRowArray = grid[0]
        subsetColArray = grid[1]
        del grid
        subsetDistArray = npy.sqrt(npy.multiply(subsetRowArray - subsetCenterRow, subsetRowArray- subsetCenterRow) + npy.multiply(subsetColArray-subsetCenterCol, subsetColArray-subsetCenterCol))           

        # Array with distance function- not yet implemented
        if options['distEq'] is not None:
            dist = subsetDistArray
            distanceFunctionArray = eval(options['distEq'])
            if len(sourceArray) < 10:
                print 'Dist'
                print dist
                print'distance Function'
                print distanceFunctionArray
            del dist
        else:
            distanceFunctionArray = npy.ones(sourceArray.shape,dtype = 'float64') #FIXME! Need distance function. Right now just doing min and max. 

        #Array to modify source array- includes distance function, max and min distances
        if not options['maxDist']:
            distanceModifierArray = npy.where(subsetDistArray>= options['minDist'],distanceFunctionArray,0) #x >= None returns X
        elif not options['minDist']:
            distanceModifierArray = npy.where(subsetDistArray<= options['maxDist'],distanceFunctionArray,0) 
        else:
            distanceModifierArray = npy.where((subsetDistArray<= options['maxDist']) & (subsetDistArray>= options['minDist']),distanceFunctionArray,0)
   
        sourceArrayNew = npy.where(sourceArray>0,npy.multiply(sourceArray,distanceModifierArray),sourceArray)
        if len(sourceArray) < 10:
            print 'source, dist, distfn, distmodifier,newsource:'
            print 'sourceArray\n',sourceArray
            print 'subsetDistArray\n',subsetDistArray
            print 'distanceFunctionArray\n',distanceFunctionArray
            print 'distanceModifierArray\n',distanceModifierArray
            print 'sourceArrayNew\n',sourceArrayNew        
       
        return sourceArrayNew
        # del subsetRowArray, sourceArrayNew, subsetColArray, subsetDistArray, distanceModifierArray, distanceFunctionArray
        return sourceArray
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()
    
    
def add_random_resistances(inRaster, options):
    try:
        if options['addRandomResistances']: 
            outRandomRaster = os.path.join(options['projectDir'],'randomRaster.tif') #FIXME: give option to use any random raster
            if not os.path.exists(outRandomRaster):
                print'Creating Random raster'
                descData = arcpy.Describe(inRaster)
                extent = descData.Extent
                cellSize = descData.MeanCellHeight
                seedValue = 1
                outRandomRaster = arcpy.sa.CreateRandomRaster(seedValue, cellSize, extent) 
            outResisRaster = arcpy.sa.Plus(outRandomRaster, inRaster) 
            fileBase, fileExtension = path.splitext(options['resisRasterBase'])
            resisRasterFA = path.join(options['scratchDir'],'rand_'+fileBase+'.tif')
            outResisRaster.save(resisRasterFA)
            inRaster = resisRasterFA #Fixme: can this just be in memory isntead?
        return inRaster 
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()
    
def quantilize(raster):
    try:
        if raster is not None:       
            inArray = arcpy.RasterToNumPyArray(raster,"#","#","#",-9999)
            interval = 100.0/options['numQuantiles']
            breakSequence = seq(interval, 100-interval, interval)
            quantileArray = npy.zeros_like(inArray, dtype = 'int32')        
            quantileArray = npy.where(inArray > 0,options['numQuantiles'],quantileArray)
                       
            inVector = inArray.flatten()
            ind = npy.where(inVector > 0) 
            inVector = inVector[ind] #Eliminate nodata and zeros

            # quantilize                      
            if len(inVector) == 0:
                print 'Current array is all zero. No results to process.'
                return
            print '\nPartitioning ' + str(len(inVector)) + ' non-zero values into ' + str(options['numQuantiles']) + ' quantiles.'
            quantileBreaks = npy.percentile(inVector, breakSequence)

            # fixme: consider replacing below with an arcpy reclassify command. Would be faster and take less memory.
            # Not critical and remap table is fussy- can't be a string
            st = time.clock()
            for i in range(options['numQuantiles']-2,-1,-1):
                quantileArray = npy.where (inArray<quantileBreaks[i],i+1,quantileArray)
            st = elapsed_time(st)                                    
            
            quantileArray = npy.where(inArray < 0,-9999,quantileArray)
            quantileArray = npy.where(inArray == 0,0,quantileArray)
            del inArray 
            
            descData = arcpy.Describe(raster)
            cellSize = descData.meanCellHeight
            extent = descData.Extent
            spatialReference = descData.spatialReference
            
            pnt = arcpy.Point(extent.XMin,extent.YMin)
            quantileRaster = arcpy.NumPyArrayToRaster(quantileArray,pnt,
                                                 cellSize,cellSize,-9999)
            arcpy.DefineProjection_management(quantileRaster,spatialReference)

            del quantileArray
            return quantileRaster
        else: return None
    except: 
        print 'Failed to quantilize'
        return None
      

        
def seq(start, stop, step = 1):
    n = int(round((stop - start)/float(step)))
    if n > 1:
        return([start + step*i for i in range(n+1)])

def delete_data(file):
    # FIXME if possible: still having file lock problems with Arc.
    if file is None:
        return
    try:
        if os.path.isfile(file):
            try:
                os.remove(file)
                gc.collect()
            except:
                try:
                    arcpy.delete_management(file)
                except:
                    pass
    except:
        pass
    return

def delete_dir(dir):
    try:
        if os.path.exists(dir):
            shutil.rmtree(dir)
        return
    except:
        return
    
def write_final_maps(options,cumCurrentRaster, cumVdiffRaster,cumFlowRaster,cumSourceRaster,cumTargetRaster):
    try:
        print 'Creating final output maps'  
        quantCurrentRaster = quantVdiffRaster = quantFlowRaster = None
        if options['quantilizeMaps']:
            quantCurrentRaster = quantilize(cumCurrentRaster)
            # quantVdiffRaster = quantilize(cumVdiffRaster) 
            # quantFlowRaster = quantilize(cumFlowRaster)

        if options['calcCurrent']:
            sumFile = path.join(options['outputDir'], 'cur_' + options['outputFileText']+'.tif')
            print 'Writing:\n',sumFile
            outRasterString = sumFile
            cumCurrentRaster.save(sumFile)
            del cumCurrentRaster
            if options['cleanUpBandFiles']:
                delete_data(options['prevCumCurrentFile'])
            if quantCurrentRaster is not None:
                sumFile = path.join(options['outputDir'], 'cur_PCT_' + options['outputFileText']+'.tif')
                outRasterString = outRasterString + '; ' + sumFile
                print 'Writing:\n',sumFile
                quantCurrentRaster.save(sumFile)
                del quantCurrentRaster
        
        if options['calcFA']:       
            sumFile = path.join(options['outputDir'], 'flow_' + options['outputFileText']+'.tif')
            if not options['calcCurrent']:
                outRasterString = sumFile
            else:
                outRasterString = outRasterString + '; ' + sumFile
            cumFlowRaster3 = arcpy.sa.Con(cumFlowRaster > 0,arcpy.sa.Int(cumFlowRaster))#(arcpy.Raster(cumFlowRaster))
            print 'Writing:\n',sumFile
            cumFlowRaster3.save(sumFile)
            delete_data(options['prevFlowFile'])
            if quantFlowRaster is not None:
                sumFile = path.join(options['outputDir'], 'flow_PCT_' + options['outputFileText']+'.tif')
                outRasterString = outRasterString + '; ' + sumFile
                print 'Writing:\n',sumFile
                quantFlowRaster.save(sumFile)
            
        if options['calcVoltages']:
            sumFile = path.join(options['outputDir'], 'vDiff_' + options['outputFileText']+'.tif')
            outRasterString = outRasterString + '; ' + sumFile
            print 'Writing:\n',sumFile
            cumVdiffRaster.save(sumFile)
            delete_data(options['prevVdiffFile'])    
            if quantVdiffRaster is not None:
                sumFile = path.join(options['outputDir'], 'vDiff_PCT_' + options['outputFileText']+'.tif')
                outRasterString = outRasterString + '; ' + sumFile
                print 'Writing:\n',sumFile
                quantVdiffRaster.save(sumFile)    
        if cumSourceRaster is not None:
            sumFile = path.join(options['outputDir'], 'sources_' + options['outputFileText']+'.tif')
            print 'Writing:\n',sumFile
            cumSourceRaster.save(sumFile)
        if cumTargetRaster is not None:    
            sumFile = path.join(options['outputDir'], 'targets_' + options['outputFileText']+'.tif')
            print 'Writing:\n',sumFile
            cumTargetRaster.save(sumFile)
        print '\nSaved output files (you can just copy and paste entire path into ArcMap add data dialog):'
        print outRasterString,'\n'
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()
    
    
@profile    
def write_temp_maps(options,bandNum,cumCurrentRaster,cumVdiffRaster,cumFlowRaster): 
    try:
        print 'Writing temporary grids...\n'
        cumCurrentFile = cumVdiffFile = cumFlowFile = None
        if options['calcCurrent']:
            cumCurrentFile = os.path.join(options['outputDir'], 'BAND'+str(bandNum)+'cur_' + options['outputFileText']+'.tif')
            try:
                cumCurrentRaster.save(cumCurrentFile)
            except:
                try:
                    print 'Error writing'
                    time.sleep(5)                
                    cumCurrentFile = os.path.join(options['outputDir'], 'BAND'+str(bandNum)+'cur_' + options['outputFileText']+'2.tif')
                    cumCurrentRaster.save(cumCurrentFile)

                except:
                    print 'Second error writing- may need to reinstate creation of new scratch directory if fails again'
                    time.sleep(5)
                    cumCurrentFile = os.path.join(options['outputDir'], 'BAND'+str(bandNum)+'cur_' + options['outputFileText']+'3.tif')
                    cumCurrentRaster.save(cumCurrentFile)                       
                    
        if options['calcVoltages']:
            cumVdiffFile = path.join(options['outputDir'], 'BAND'+str(bandNum)+'vDiff_' + options['outputFileText']+'.tif')
            try:
                cumVdiffRaster.save(cumVdiffFile)
            except:    
                try:
                    print 'Error writing'
                    time.sleep(5)
                    cumVdiffFile = os.path.join(options['outputDir'],'BAND'+str(bandNum)+options['resisRasterBase']+'_vdiffSumtemp2.tif')
                    cumVdiffRaster.save(cumVdiffFile)
                except:
                    print 'Second error writing'
                    time.sleep(5)
                    cumVdiffFile = os.path.join(options['outputDir'],'BAND'+str(bandNum)+options['resisRasterBase']+'_vdiffSumtemp3.tif')
                    cumVdiffRaster.save(cumVdiffFile)          

        if options['calcFA']:
            cumFlowFile = os.path.join(options['outputDir'], 'BAND'+str(bandNum)+'flow_' + options['outputFileText']+'.tif')
            try:
                cumFlowRaster.save(cumFlowFile)
            except:
                try:
                    time.sleep(5)
                    cumFlowFile = os.path.join(options['outputDir'], 'BAND'+str(bandNum)+'flow_' + options['outputFileText']+'a.tif')
                    cumFlowRaster.save(cumFlowFile)
                except:
                    time.sleep(5)
                    cumFlowFile = os.path.join(options['outputDir'], 'BAND'+str(bandNum)+'flow_' + options['outputFileText']+'b.tif')
                    cumFlowRaster.save(cumFlowFile)
        if options['cleanUpBandFiles']:
            delete_data(options['prevCumCurrentFile'])
        options['prevCumCurrentFile'] = cumCurrentFile
        delete_data(options['prevVdiffFile'])
        options['prevVdiffFile'] = cumVdiffFile
        delete_data(options['prevFlowFile'])
        options['prevFlowFile'] = cumFlowFile        
                
        return options
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()
    
            
def clean_up(options):
    if not options['cleanUpBandFiles']:
        return
    print 'Cleaning up...'
    # arcpy.RefreshCatalog(options['outputDir'])
    # arcpy.RefreshCatalog(options['scratchDir'])
    # arcpy.env.Workspace = 'c:\\temp'
    # arcpy.env.scratchWorkspace = 'c:\\temp'
    try:
        if os.path.exists(options['scratchDir']):
            shutil.rmtree(dir)
    except:
        for filename in glob.glob(os.path.join(options['scratchDir'],'*.*')):
            delete_data( filename )
        try:
            if os.path.exists(options['scratchDir']):
                shutil.rmtree(options['scratchDir'])
        except:
            pass
    for filename in glob.glob(os.path.join(options['outputDir'],'band*.*')):
        try:    
            delete_data( filename )
        except:
            pass


def set_options_and_dirs(options):
    try:
        """derives output base filename by combining input options. Also handles filenames for tiling. Checks and corrects options    """    

        if options['tileNum'] > 0:
            if options['tileNum'] > options['numTiles']:
                print 'Error: tile number (tileNum) cannot be higher than number of tiles (numTiles).'
                raw_input('\nPress Enter to continue.') 
                exit(0)    

            resisRaster = path.join(options['projectDir'],options['resisRasterBase'])
            header = get_header(resisRaster)
            approxEndBand = int(header['nrows']/options['blockSize'])+1   
            tileSize = int(approxEndBand/options['numTiles'])+1
            print 'There are approximately',approxEndBand,'bands in the dataset.'
            print 'Tile size is',tileSize,'bands.'
            
            options['startBand'] = (tileNum-1) * tileSize + 1
            options['endBand'] = (tileNum) * tileSize 
            if options['endBand']>approxEndBand:
                options['endBand'] = 0
            print 'startBand is',options['startBand']
            print 'endBand is',options['endBand']
            
        
        if options['limitCalcsByCWD'] and not options['calcFA']:
            print 'Error: cannot limit calculations by CWD unless calcFA is set to True.'
            raw_input('\nPress Enter to continue.') 
            exit(0)    
        if not options['calcCurrent']:
            options['calcVoltages'] = False
        # if not options['useClimate']:
            # options['saveSourceAndTargetCounts'] = False
        if not options['calcCurrent'] and not options['calcFA']:
            print 'Error: must either calculate current or flow accumulation. Both are set to False.'
            raw_input('\nPress Enter to continue.') 
            exit(0)
        
        if options['calcNull']:
            options['outputDirBase'] = options['outputDirBase']+'_NULL'

        options['scratchDirBase'] = 'scratch'+options['outputDirBase'] 
        if options['tileNum'] > 0:
            options['scratchDirBase'] = options['scratchDirBase'] + str(options['tileNum'])

        options['compress'] = False
        
        # if tile >= 0: 
            # options['outputDirBase'] = options['outputDirBase']+str(tile)
            # options['scratchDirBase'] = options['scratchDirBase']+str(tile)
            # fileBase,ext = os.path.splitext(options['resisRasterBase'])
            # options['resisRasterBase'] = fileBase+str(tile)+ext
            # if options['useClimate']:
                # if not options['matchClimatePCs']:
                    # fileBase,ext = os.path.splitext(options['climateRasterBase'])
                    # options['climateRasterBase'] = fileBase+str(tile)+ext
                # else:
                    # fileBase,ext = os.path.splitext(options['t1PC1RasterBase'])
                    # options['t1PC1RasterBase'] = fileBase+str(tile)+ext
            
                    # fileBase,ext = os.path.splitext(options['t2PC1RasterBase'])
                    # options['t2PC1RasterBase'] = fileBase+str(tile)+ext

                    # fileBase,ext = os.path.splitext(options['t1PC2RasterBase'])
                    # options['t1PC2RasterBase'] = fileBase+str(tile)+ext

                    # fileBase,ext = os.path.splitext(options['t2PC2RasterBase'])
                    # options['t2PC2RasterBase'] = fileBase+str(tile)+ext
                    
            # if options['useSourceRaster']:
                # fileBase,ext = os.path.splitext(options['sourceRasterBase'])    
                # options['sourceRasterBase'] = fileBase+str(tile)+ext
                
        if options['startBand'] is None: options['startBand'] = 0 
        if options['startStripe'] is None: options['startStripe'] = 0
        if options['endBand'] is None: options['endBand'] = 0
        if options['endStripe'] is None: options['endStripe'] = 0
        
        if options['radius'] <= options['blockSize']/2:
            print 'Error. Radius must be larger than half the block size.'
            raw_input('\nPress Enter to continue.') 
            exit(0)

        if float(options['blockSize'])/2 == int(options['blockSize']/2):
            print float(options['blockSize'])/2
            print int(options['blockSize']/2)
            print 'Error. Block size must be an odd number.'
            raw_input('\nPress Enter to continue.') 
            exit(0)
            
        resisRasText, fileExtension = os.path.splitext(options['resisRasterBase'])
        # if not options['useSourceRaster']:
        if len(resisRasText)>25:
                resisRasText = resisRasText[0:24]+'x'                    
        # elif len(resisRasText)>15:
            # resisRasText = resisRasText[0:14]+'x'                    

        if options['squareResistances']:
            squareText = '_sq'
        else:    
            squareText = ''
        radiusText = '_r'+str(options['radius'])
        if options['useSourceRaster']:
            fileBase,ext = os.path.splitext(options['sourceRasterBase'])
            srcText = 'src_'+ fileBase + '_'
            if len(srcText)>25:
                srcText = srcText[0:24]+'x' 
        else:
            srcText = ''

        if options['limitCalcsByCWD']:
            limText = '_lim' + str(options['cwdLimit'])
        else:
            limText = ''
        if options['negTargets']:
            negTargText = 'negTarg'
        else:
            negTargText = ''

        if options['fadeIn']:
            fadeInText = '_f'+str(int(options['fadeInDist']))
        else:
            fadeInText = ''

        if options['startBand'] > 0:
            startBandText = '_SB'+str(options['startBand'])
        else:
            startBandText = ''
        if options['endBand'] > 0:
            endBandText = '_EB'+str(options['endBand'])
        else:
            endBandText = ''
        if options['startStripe'] > 0:
            startStripeText = '_SS'+str(options['startStripe'])
        else:
            startStripeText = ''
        if options['endStripe'] > 0:
            endStripeText = '_ES'+str(options['endStripe'])
        else:
            endStripeText = ''
            
        distFunctionText = ''
        if options['useDistanceFunction'] and (options['maxDist'] or options['minDist'] or options['distEq'] is not None):
            if options['minDist']:
                distFunctionText = '_minDist'+str(options['minDist']) 
            if options['maxDist']:
                distFunctionText = distFunctionText + '_maxDist'+str(options['maxDist']) 
            if options['distEq'] is not None:
                distFunctionText = distFunctionText + '_distFn'
        if options['weightTargOnly']:
            weightText = '_targOnly'
        else:
            weightText = ''
        if options['centerGround']:
            centerText = 'cg'
        else:
            centerText = ''
        
        if options['calcMaxCur']:
            maxText = '_max'
        else:
            maxText = ''
        if options['noWeight']:
            noWeightText = 'noWeight'
        else:
            noWeightText = ''

        if options['useClimate']:
            if not options['matchClimatePCs']:

                fileBase,ext = os.path.splitext(options['climateRasterBase'])
                if len(fileBase)>15:
                    fileBase = fileBase[0:14]+'x'
                climText = 'clim_tc_'+str(options['tDiff']).replace('.','')+'_'+fileBase
                if options['absClim']:
                    climText = climText+'_abs'
            else:
                climText = '_clim_PCs_Cut'+str(options['PCWindow']).replace('.','')
        else:
            climText = ''
        if options['calcNull']:
            nullText = 'NULL_'
        else:
            nullText = ''

        if not options['useSourceRaster']:
            if options['sourceInverseR']:
                cutoffText = '_srcInvR'
            else:    
                cutoffText = '_rc'+str(options['rCutoff'])
        else:
            cutoffText = ''
        options['outputFileText'] = nullText+resisRasText + squareText + limText + '_' + srcText+climText+radiusText+'b'+str(options['blockSize'])+cutoffText +distFunctionText+weightText+noWeightText+centerText+maxText+startBandText+endBandText+startStripeText+endStripeText+negTargText+fadeInText
        options['prevFlowFile'] = None
        options['prevCumCurrentFile'] = None
        options['prevVdiffFile'] = None

        print 'resisRasterBase = ',options['resisRasterBase']
        print 'Radius = ',options['radius']
        print 'BlockSize = ',options['blockSize']
        if not options['useSourceRaster']:
            if options['sourceInverseR']:
                print 'Using inverse of resistances for sources'
            else:
                print'rCutoff = ',options['rCutoff']
        else:
            print 'sourceRasterBase = ',options['sourceRasterBase']
        print'startBand',options['startBand']
        print'endBand',options['endBand']
        print'startStripe',options['startStripe']
        print'endStripe',options['endStripe']
        if options['calcNull']:
            print'calcNull',options['calcNull']
        if options['fadeIn']:
            print 'FadeIn Distance',options['fadeInDist']
        if options['useDistanceFunction']:
            print 'minDist',options['minDist'] 
            print 'maxDist',options['maxDist'] 

        options['outputDir'] = os.path.join(options['projectDir'],options['outputDirBase'])
        options['scratchDir'] = os.path.join(options['projectDir'],options['scratchDirBase'])
        print 'project Dir',options['projectDir']
        print 'output Dir',options['outputDirBase']
        delete_dir(options['scratchDir'])
        if not path.exists(options['scratchDir']):
            os.mkdir(options['scratchDir'])
        # else:
            # for i in range(1,50):
                # scratchDir = options['scratchDir'] + str(i)
                # if not path.exists(scratchDir):
                    # options['scratchDir'] = scratchDir
                    # os.mkdir(options['scratchDir'])   
                    # break
        if not path.exists(options['outputDir']):
            os.mkdir(options['outputDir'])
        print '\n'
        return options
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

def copy_this_file(options):
    try:
        # Save a copy of this file in output directory
        # destFile = os.path.join(options['outputDir'],'omniscape_'+options['outputFileText']+'.py')
        ft = tuple(time.localtime())
        timeNow = time.ctime()
        fileName = ('%s_%s_%s_%s%s_%s' % (ft[0], ft[1], ft[2], ft[3], ft[4], 'omniscape_'+options['outputFileText']+'.py'))# os.path.basename(sys.argv[0])))
        filePath = os.path.join(options['outputDir'],fileName)
        shutil.copyfile(sys.argv[0],filePath) 
        # shutil.copyfile(sys.argv[0],destFile) 
    except:    
        print_python_error()

def get_max_vdiff(voltMap,circleResisArray,csOptions):
    try:
        voltMap_l, voltMap_r = get_horiz_neighbors(voltMap)
        voltMap_u, voltMap_d = get_vert_neighbors(voltMap)
        # print 'voltMap_l'
        # print voltMap_l
        if options['adjustVoltages']: 
            # fixme: not completed. Aim is to 
            # see how much improvement is possible if r reduced to 1.
            # vdiff * (r-1)/r... if r = 1, no improvement possible. if r = 2, half of voltage could be reduced. etc.
            # need r's, will include r from 3 pixels. but try to calc how much v would drop if that one pixel were restored.
            # HOW MUCH WOULD V DROP ACROSS PIXEL IF R>1
            
            # But know current, so v = ir to get r?
            # 10 A going thru cell.
            # volt = 50
            # R = V/I
            # assumed R = 5
            
            
            
            # figure out where max vdiff is
            # figure out resistance across that
            
            resis_l, resis_r = get_horiz_neighbors(circleResisArray)
            resis_u, resis_d = get_vert_neighbors(circleResisArray)
            
            
            print 'adjustVoltages code not finished, nuances with resistance calcs'
            blarg
            
        vdiff_N = get_vdiff(voltMap_d,voltMap_u)
        vdiff_E = get_vdiff(voltMap_l,voltMap_r)
        del voltMap_u, voltMap_d, voltMap_l, voltMap_r 
        max_vdiff = npy.fmax(vdiff_N,vdiff_E)
        del vdiff_N,vdiff_E
        if csOptions['connect_four_neighbors_only'] == False:
            voltMap_ul, voltMap_dr = get_diag1_neighbors(voltMap)
            voltMap_ur, voltMap_dl = get_diag2_neighbors(voltMap)
            vdiff_NE = get_vdiff(voltMap_dl,voltMap_ur)
            vdiff_SE = get_vdiff(voltMap_ul,voltMap_dr)
            max_vdiff = npy.fmax(max_vdiff,vdiff_NE)
            max_vdiff = npy.fmax(max_vdiff,vdiff_SE)
            del voltMap_ul, voltMap_dr, voltMap_ur, voltMap_dl, vdiff_NE,vdiff_SE
        return max_vdiff
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()
    
def get_vdiff(voltMap,voltMap_direction):
    try:
        #FIXME: vdiff code not working right, debug using 6x6
        vdiff_direction = npy.absolute(voltMap - voltMap_direction)
        vdiff_direction = npy.where(voltMap == npy.nan, 0, npy.where(voltMap_direction == npy.nan,0,vdiff_direction))
        return vdiff_direction
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

def get_horiz_neighbors(map):
    try:
        m = map.shape[0]
        n = map.shape[1]
        
        # zeromap = npy.zeros(map.shape,dtype = 'float64') 
        zeromap = npy.empty(map.shape,dtype = 'float64') 

        zeromap[:] = npy.NAN
        map_l = zeromap.copy()
        map_r = zeromap.copy()
        
        map_l[:,1:n] = map[:, 0:(n-1)]
        map_r[:,0:(n-1)] = map[:, 1:n]
        return map_l, map_r
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

def get_vert_neighbors(map):
    try:
        m = map.shape[0]
        n = map.shape[1]
        zeromap = npy.zeros(map.shape,dtype = 'float64') 
        map_u = zeromap.copy()
        map_d = zeromap.copy()

        map_u[1:m, :] = map[0:(m-1), :]
        map_d[0:(m-1) , :] = map[1:m , :]
        
        return map_u, map_d
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

def get_diag1_neighbors(map):
    try:
        m = map.shape[0]
        n = map.shape[1]   
        zeromap = npy.zeros(map.shape,dtype = 'float64') 
        map_ul = zeromap.copy()
        map_ul[1:m,1:n] = map[0:m-1, 0:n-1]
        map_dr = zeromap.copy()
        map_dr[0:m-1, 0:n-1  ] = map[1:m , 1:n ]
        return map_ul, map_dr
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()


def get_diag2_neighbors(map):
    try:
        m = map.shape[0]
        n = map.shape[1]
        zeromap = npy.zeros(map.shape,dtype = 'float64') 
        map_ur = zeromap.copy()
        map_ur[1:m,0:n-1] = map[0:m-1, 1:n  ]
        map_dl = zeromap.copy()
        map_dl[0:m-1, 1:n  ] = map[1:m  , 0:n-1]
        return map_ur, map_dl
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

# @profile
def band(inRaster,header,centerRow, options): 
    try:
        bandRows = 1 + min(options['radius'],centerRow) + min(header['nrows'] - (centerRow+1), options['radius'])
        yMin = max(header['yllcorner'],header['yllcorner'] + ((header['nrows'] - centerRow - options['radius'] - 1) * header['cellsize']))
        LLC = arcpy.Point(header['xllcorner'],yMin)
        bandArray = arcpy.RasterToNumPyArray(inRaster,LLC,"#",bandRows,-9999)  
        if bandArray.shape[1] > header['ncols']: #FIXME: temp fix because random raster yields bands that are 1 too long for some reason      
            print'Band array is too long'
            print 'ncols',header['ncols']
            print 'shape',bandArray.shape[1]
            print 'last col min,max', npy.min(bandArray[:,header['ncols']]), npy.max(bandArray[:,header['ncols']])
            print 'Col 0 min,max', npy.min(bandArray[:,0]), npy.max(bandArray[:,0])
            # bandArray = bandArray[:,1:header['ncols']] #CHANGED BHM 12/15/15
            bandArray = bandArray[:,1:bandArray.shape[1]] #CHANGED BHM 12/15/15  
        return bandArray
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

@profile
def addData(cumCurrentArray, currentArray, centerRow, centerCol, options):
    try:
        #xprint 'adding data for center row, col, rad'
        #xprint centerRow,centerCol,options['radius']
        minRow = max(0,centerRow-options['radius'])
        maxRow = min(minRow+cumCurrentArray.shape[0]-1, centerRow+options['radius'])
        minCol = max(0,centerCol-options['radius'])
        # fixme: check if next line needs mincol added like in maxrow line
        maxCol = min(cumCurrentArray.shape[1]-1, centerCol+options['radius'])

        fullCurrentArray = npy.zeros(cumCurrentArray.shape,dtype = 'float64')
        fullCurrentArray[minRow:maxRow+1,minCol:maxCol+1] = currentArray
        if options['calcMaxCur']:
            cumCurrentArray = npy.where(cumCurrentArray < fullCurrentArray, fullCurrentArray, cumCurrentArray)
        else:
            cumCurrentArray += fullCurrentArray
        return cumCurrentArray
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

@profile
def addData_arcpy(cumCurrentRaster, currentRaster):
    try:
        descData = arcpy.Describe(cumCurrentRaster)
        arcpy.env.extent = descData.Extent    
        arcpy.env.extent = "MAXOF"
        #fixme: put statements below in single nested con
        newCumCurrentRaster = arcpy.sa.Con(arcpy.sa.IsNull(currentRaster),cumCurrentRaster,arcpy.sa.Plus(currentRaster,cumCurrentRaster))
        newCumCurrentRaster2 = arcpy.sa.Con(arcpy.sa.IsNull(newCumCurrentRaster),currentRaster,newCumCurrentRaster)
        return newCumCurrentRaster2
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

        
def trim(array, centerRow, centerCol, nrows, ncols, options):
        radius = options['radius']
        startRow=startCol=0
        if centerRow<radius:
            startRow = radius-centerRow
            # print 'startRow',startRow
        if centerCol<radius:
            startCol = radius-centerCol
            # print 'startCol',startCol
        endRow = array.shape[0]-1
        endCol = array.shape[1]-1
        if centerRow + radius > nrows-1:
            nTrimRows = centerRow + radius - (nrows-1)
            # print 'nTrimRows',nTrimRows
            endRow = endRow - nTrimRows
            # print 'ENDROW=',endRow
        if centerCol + radius > ncols-1:
            nTrimCols = centerCol + radius - (ncols-1)
            # print 'nTrimCols',nTrimCols
            endCol = endCol - nTrimCols
            # print 'ENDCOL=',endCol
            
        # startCol = max(centerCol - options['radius'],0)
        # endCol = min(centerCol + options['radius'],array.shape[1]-1)
        arraySmall = array[startRow:endRow+1,startCol:endCol+1]

        return arraySmall
    



@profile
def circ(array, rowArray, colArray, centerRow, centerCol, radius):
    try:
        startRow = max(centerRow - radius,0)
        endRow = min(centerRow + radius,array.shape[0]-1)
        startCol = max(centerCol - radius,0)
        endCol = min(centerCol + radius,array.shape[1]-1)
        colArraySmall = (colArray[startRow:endRow+1,startCol:endCol+1]).astype('float64')
        rowArraySmall = (rowArray[startRow:endRow+1,startCol:endCol+1]).astype('float64')
        arraySmall = array[startRow:endRow+1,startCol:endCol+1]
        
        distArray = npy.sqrt(npy.multiply((rowArraySmall - centerRow), (rowArraySmall- centerRow)) + npy.multiply((colArraySmall-centerCol), (colArraySmall-centerCol))) #causes memory error

        del rowArraySmall,colArraySmall
        arrayMasked = npy.where(distArray <= radius, arraySmall, -9999) #fixme do nans here
        gc.collect()

        return arrayMasked
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

# @profile        
def center_block(array, blockSize, centerRow, centerCol):
    try:
        # returns array of same shape as array (i.e. entire radius), but everything outside of center blkock is zero
        startRow = centerRow - ((blockSize-1)/2)
        endRow = centerRow + ((blockSize-1)/2)
        startCol = centerCol - ((blockSize-1)/2)
        endCol = centerCol + ((blockSize-1)/2)
        blockArray = npy.zeros(array.shape, dtype = 'float64')# - 9999  # replace -9999 with nan?  
        blockArray[startRow:endRow+1,startCol:endCol+1] = array[startRow:endRow+1,startCol:endCol+1]
        return blockArray            
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

def get_subset_header(array, fullHeader, options, centerRow, centerCol):
    try:
        subsetHeader = {}
        llrow = min(centerRow + options['radius'], fullHeader['nrows']-1)
        llcol = max(centerCol - options['radius'], 0)
       
        diffy = fullHeader['nrows'] - 1 - llrow
        diffx = llcol
        yllcorner = fullHeader['yllcorner'] + fullHeader['cellsize'] * diffy
        xllcorner = fullHeader['xllcorner'] + fullHeader['cellsize'] * diffx

        subsetHeader['ncols'] = array.shape[1]
        subsetHeader['nrows'] = array.shape[0]
        subsetHeader['xllcorner'] = xllcorner
        subsetHeader['yllcorner'] = yllcorner
        subsetHeader['cellsize'] = fullHeader['cellsize']
        subsetHeader['nodata'] = fullHeader['nodata']
        return subsetHeader
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

def get_cs_path():
    """Returns path to Circuitscape installation """
    try:
        envList = ["ProgramW6432", "ProgramFiles", "ProgramFiles(x86)"]
        for x in range (0,len(envList)):
            try:
                pfPath = os.environ[envList[x]]
                csPath = os.path.join(pfPath,'Circuitscape\\cs_run.exe')
                if os.path.exists(csPath): return csPath
            except: pass
        return 'D:\\Program Files\\Circuitscape\\cs_run.exe'
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()
    
def elapsed_time(start_time):
    """Returns elapsed time given a start time"""
    try:
        now = time.clock()
        elapsed = now - start_time
        secs = float((int(100*elapsed)))/100
        mins = int(elapsed / 60)
        hours = int(mins / 60)
        mins = mins - hours * 60
        secs = secs - mins * 60 - hours * 3600
        if mins == 0:
            print('That took ' + str(secs) + ' seconds.\n')
        elif hours == 0:
            print('That took ' + str(mins) + ' minutes and ' +
                              str(secs) + ' seconds.\n')
        else:
            print('That took ' + str(hours) + ' hours ' +
                              str(mins) + ' minutes and ' + str(secs) +
                              ' seconds.\n')
        return now
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

@profile        
def ascii_grid_reader(filename, nodata,data_type):
    """Reads rasters saved as ASCII grids or numpy arrays into Circuitscape."""
    try:   
        if nodata == False:
            pmap = npy.loadtxt(filename, skiprows = 5, dtype = data_type)
        else:
            pmap = npy.loadtxt(filename, skiprows = 6, dtype = data_type)
            pmap = npy.where(pmap == nodata, -9999, pmap)

        return pmap
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()


def get_header(filename):
    try:
        header = {}
        descData = arcpy.Describe(filename)
        cellsize = descData.meanCellHeight
        extent = descData.Extent
        xllcorner = extent.XMin
        yllcorner = extent.YMin
        yulcorner = extent.YMax
        xlrcorner = extent.XMax
        nrows = int((yulcorner-yllcorner)/cellsize)
        ncols = int((xlrcorner-xllcorner)/cellsize)
        header = {}
        header['ncols'] = ncols
        header['nrows'] = nrows
        header['xllcorner'] = xllcorner
        header['yllcorner'] = yllcorner
        header['cellsize'] = cellsize
        header['nodata'] = -9999 #if (nodata == False) else nodata 
        return header 
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

@profile    
def ascii_grid_writer(file_name, data, header, options):
    """Writes rasters to ASCII grid or numpy formats."""     
    try:
        f = gzip.open(file_name+'.gz', 'w') if options['compress'] else open(file_name, 'w')    
        f.write('ncols         ' + str(header['ncols']) + '\n')
        f.write('nrows         ' + str(header['nrows']) + '\n')
        f.write('xllcorner     ' + str(header['xllcorner']) + '\n')
        f.write('yllcorner     ' + str(header['yllcorner']) + '\n')
        f.write('cellsize      ' + str(header['cellsize']) + '\n')
        f.write('NODATA_value  ' + str(header['nodata']) + '\n')
         
        delimiter = ''
        fmt = ['%.10g ']*header['ncols']
        fmt = delimiter.join(fmt)
        fmt += '\n'
        for row in data:
            f.write(fmt % tuple(row))     
        f.close()
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

@profile    
def export_ras_to_npy(raster,npyFile):
    try:    
        descData = arcpy.Describe(raster)
        cellSize = descData.meanCellHeight
        extent = descData.Extent
        spatialReference = descData.spatialReference
        
        pnt = arcpy.Point(extent.XMin,extent.YMin)
        outData = arcpy.RasterToNumPyArray(raster,"#","#","#",-9999)
        #outData = npy.where(outData == noDataVal,-9999,outData)
        if npy.array_equiv(outData, outData.astype('int32')):
            outData = outData.astype('int32')
        npy.save(npyFile, outData)
        write_header(raster,outData,npyFile)
                
        numElements = (outData.shape[0] * outData.shape[1])
        #rows,cols = npy.where(outData != -9999)
        numNodes = (npy.where(outData != -9999, 1, 0)).sum() 
        #numZeros = (npy.where(outData != -9999, 1, 0)).sum() 
        #del rows
        
        del outData
        return numElements, numNodes
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

@profile
def import_npy_to_ras(npyFile,baseRaster,outRasterPath):
    try:
        npyArray = npy.load(npyFile, mmap_mode = None)
        npyArray = npyArray.astype('float32')
        descData = arcpy.Describe(baseRaster)
        cellSize = descData.meanCellHeight
        extent = descData.Extent
        spatialReference = descData.spatialReference
        
        pnt = arcpy.Point(extent.XMin,extent.YMin)
        newRaster = arcpy.NumPyArrayToRaster(npyArray,pnt,
                                             cellSize,cellSize,-9999)
        if outRasterPath is not None:
            newRaster.save(outRasterPath)
            return
        else:
            return newRaster
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()
        

def write_header(raster,numpyArray,numpyFile):
    try:    
        ncols = numpyArray.shape[1]
        nrows = numpyArray.shape[0]
        descData = arcpy.Describe(raster)
        cellSize = descData.meanCellHeight
        extent = descData.Extent
        xllcorner = extent.XMin
        yllcorner = extent.YMin
        nodata = -9999
        fileBase, fileExtension = path.splitext(numpyFile)
        headerFile = fileBase + '.hdr'
        f = False
        f = open(headerFile, 'w')

        f.write('ncols         ' + str(ncols) + '\n')
        f.write('nrows         ' + str(nrows) + '\n')
        f.write('xllcorner     ' + str(xllcorner) + '\n')
        f.write('yllcorner     ' + str(yllcorner) + '\n')
        f.write('cellsize      ' + str(cellSize) + '\n')
        f.write('NODATA_value  ' + str(nodata) + '\n')

        f.close()
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

def print_failure(numResistanceNodes, memFlag, sleepTime):
    try:
        print('\nCircuitscape failed. See error information above.')
        if memFlag == True:
            totMem, availMem = lu.get_mem()                    
            print('Note: Circuitscape can only solve 2-3 million nodes')
            print('per gigabyte of available RAM. Your resistance raster had ')
            print(str(int(numResistanceNodes)) + ' nodes.\n')  
            print('Total physical RAM on your machine is ~' 
                   + str(totMem) 
                   + ' GB. \nAvailable memory is ~'
                   + str(availMem) + ' GB. \n')
        print('Trying again in ' + str(sleepTime) + ' seconds.')
        lu.snooze(sleepTime)                    
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()


def setCircuitscapeOptions(resisAsciiFile,sourceAsciiFile,groundAsciiFile,outputFile):
    """Sets default options for calling Circuitscape.

    """
    try:
        csOptions = {}
        csOptions['data_type'] = 'raster'
        csOptions['version'] = 'unknown'
        csOptions['low_memory_mode'] = False
        csOptions['scenario'] = 'advanced'
        csOptions['habitat_file'] = '(Browse for a habitat map file)'
        csOptions['habitat_map_is_resistances'] = True
        csOptions['point_file'] = ('(Browse for file with '
                              'locations of focal points or areas)')
        csOptions['point_file_contains_polygons'] = True
        csOptions['connect_four_neighbors_only'] = False
        csOptions['connect_using_avg_resistances'] = True
        csOptions['use_polygons'] = False
        csOptions['polygon_file'] = '(Browse for a short-circuit region file)'
        csOptions['source_file'] = '(Browse for a current source file)'
        csOptions['ground_file'] = '(Browse for a ground point file)'
        csOptions['ground_file_is_resistances'] = True
        csOptions['use_unit_currents'] = False
        csOptions['use_direct_grounds'] = False
        csOptions['remove_src_or_gnd'] = 'rmvsrc'
        csOptions['output_file'] = '(Choose a base name for output files)'
        csOptions['write_cur_maps'] = True
        csOptions['write_cum_cur_map_only'] = True
        csOptions['log_transform_maps'] = False
        csOptions['write_volt_maps'] = False
        csOptions['solver'] = 'cg+amg'
        csOptions['compress_grids'] = False
        csOptions['print_timings'] = False
        csOptions['use_mask'] = False
        csOptions['mask_file'] = 'None'
        csOptions['use_included_pairs'] = False
        csOptions['included_pairs_file'] = 'None'
        csOptions['use_variable_source_strengths'] = False
        csOptions['variable_source_file'] = 'None'
        csOptions['write_max_cur_maps'] = False
        csOptions['set_focal_node_currents_to_zero'] = True
        csOptions['print_timings'] = False
        csOptions['output_file'] = outputFile
           
        csOptions['remove_src_or_gnd'] = 'rmvsrc'#'keepall'
        csOptions['habitat_file'] = resisAsciiFile
        csOptions['source_file'] = sourceAsciiFile
        csOptions['ground_file'] = groundAsciiFile

        return csOptions
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

        
def writeCircuitscapeConfigFile(configFile, csOptions):
    """Creates a configuration file for calling Circuitscape.

    """
    try:
        config = ConfigParser.ConfigParser()

        sections = {}
        section = 'Version'
        sections['version'] = section

        section = 'Connection scheme for raster habitat data'
        sections['connect_four_neighbors_only'] = section
        sections['connect_using_avg_resistances'] = section

        section = 'Short circuit regions (aka polygons)'
        sections['use_polygons'] = section
        sections['polygon_file'] = section

        section = 'Options for advanced mode'
        sections['source_file'] = section
        sections['ground_file'] = section
        sections['ground_file_is_resistances'] = section
        sections['use_unit_currents'] = section
        sections['use_direct_grounds'] = section
        sections['remove_src_or_gnd'] = section

        section = 'Calculation options'
        sections['solver'] = section
        sections['print_timings'] = section
        sections['low_memory_mode'] = section

        section = 'Output options'
        sections['output_file'] = section
        sections['write_cur_maps'] = section
        sections['write_cum_cur_map_only'] = section
        sections['log_transform_maps'] = section
        sections['write_volt_maps'] = section
        sections['compress_grids'] = section
        sections['write_max_cur_maps'] = section
        sections['set_focal_node_currents_to_zero'] = section

        section = 'Mask file'
        sections['use_mask'] = section
        sections['mask_file'] = section

        section = 'Options for pairwise and one-to-all and all-to-one modes'
        sections['use_included_pairs'] = section
        sections['included_pairs_file'] = section
        sections['point_file'] = section
        sections['point_file_contains_polygons'] = section

        section = 'Options for one-to-all and all-to-one modes'
        sections['use_variable_source_strengths'] = section
        sections['variable_source_file'] = section

        section = 'Habitat raster or graph'
        sections['habitat_file'] = section
        sections['habitat_map_is_resistances'] = section

        section = "Circuitscape mode"
        sections['scenario'] = section
        sections['data_type'] = section

        if csOptions['ground_file_is_resistances'] == 'not entered':
            csOptions['ground_file_is_resistances'] = False
        if csOptions['point_file_contains_polygons'] == 'not entered':
            csOptions['point_file_contains_polygons'] = False

        for csOption in sections:
            try:
                config.add_section(sections[csOption])
            except:
                pass
        for csOption in sections:
            config.set(sections[csOption], csOption, csOptions[csOption])

        f = open(configFile, 'w')
        config.write(f)
        f.close()
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()


def call_circuitscape(CSPATH, outConfigFile):
    try:    
        memFlag = False
        failFlag = False
        print('Calling Circuitscape:')
        proc = subprocess.Popen([CSPATH, outConfigFile],
                               stdout = subprocess.PIPE, stderr = subprocess.STDOUT, 
                               shell = True)
        while proc.poll() is None:
            output = proc.stdout.readline()
            
            if 'No valid sources' in output:
                print 'no valid sources detected, skipping'
                blarg
                return

            if 'Traceback' in output:
                print("\nCircuitscape failed.")
                failFlag = True
                if 'memory' in output:
                    memFlag = True
            if ('Job took' not in output and 'Reading' not in output and 'Processing' not in output and 'laplacian' not in output and 
                    'Resistance/conductance' not in output and 'module' not in output and 'node_map' not in output 
                    and (('--' in output) or ('sec' in output) or (failFlag == True))):
                print(output.replace("\r\n",""))                
        
        # Catch any output lost if process closes too quickly
        output = proc.communicate()[0]
        for line in output.split('\r\n'):
            if 'Traceback' in line:
                print("\nCircuitscape failed.")
                if 'memory' in line:
                    memFlag = True
            if ('Processing' not in line and 'laplacian' not in line and 
                    'node_map' not in line and (('--' in line) or 
                    ('sec' in line) or (failFlag == True))):
               print("      " + str(line))#.replace("\r\n","")))              
        return memFlag
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()


def create_log_file(options):
    ft = tuple(time.localtime())
    timeNow = time.ctime()
    fileName = 'omni_LOG_'+options['outputFileText']+'.txt'
    filePath = os.path.join(options['outputDir'],fileName)
    try:
        logFile = open(filePath,'a')
    except:
        logFile = open(filePath,'w')
    if inParameters is not None:
        logFile.write('*'*70 + '\n')
        logFile.write('Omniscape log file: %s \n\n' % (toolName))
        logFile.write('Start time:\t%s \n' % (timeNow))
        logFile.write('Parameters:\n', options)
    logFile.close()
    options['logFilePath'] = filePath
    return options

def write_log(string):
    try:
        logFile = open(options['logFilePath'],'a')
    except:
        logFile = open(options['logFilePath'],'w')
    try:
        #Sometimes int objects returned for arc failures so need str below
        logFile.write(str(string) + '\n')
    except IOError:
        pass
    finally:
        logFile.close()


def close_log_file():
    timeNow = time.ctime()
    try:
        write_log('\nStop time:\t\t%s \n\n' % (timeNow))
        # cfg.logFile.close() # DMK - Log closed after each write
    except:
        pass

def print_python_error():
    """Handle python errors and provide details to user"""
    import traceback
    tb = sys.exc_info()[2]  # get the traceback object
    # tbinfo contains the error's line number and the code
    tbinfo = traceback.format_tb(tb)[0]
    line = tbinfo.split(", ")[1]

    err = traceback.format_exc().splitlines()[-1]
    msg = ("Python error on **" + line + ":")
    print msg
    print err
    if options['tileNum'] > 0:
        raw_input('\nPress Enter to continue.') 
    exit(1)


def print_geoproc_error():
    """Handle geoprocessor errors and provide details to user"""
    import traceback
    tb = sys.exc_info()[2]  # get the traceback object
    # tbinfo contains the error's line number and the code
    tbinfo = traceback.format_tb(tb)[0]
    line = tbinfo.split(", ")[1]
    msg = ("ArcGIS error on **" + line + ":")
    print(msg)
    msg=arcpy.GetMessages(2)
    print(arcpy.GetMessages(2))
    if options['tileNum'] > 0:
        raw_input('\nPress Enter to continue.')    
    exit(1)

        
if __name__ == '__main__':
    options = omniscape(options)

    # Need this outside of main loop to avoid file lock problems:
    clean_up(options)
    if tileNum>0:
        raw_input('Done. Press Enter to continue.')
    print 'Done'    
