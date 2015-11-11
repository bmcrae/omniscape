# only save one omniscape.py file
# subtract out sources from FA grid? could eliminate whiskers.

# don't do int before band saves for flow, unless speed at issue
# save int and intln flow files

# Flow mapping
# break out values into 10 classes
# grow each class by value. Nah.
# Have separate line layers for each class? Easier to have one on top of the other, transparancy etc.


# go with fade... setting whole block to average of edges has its own artifacts.
# Based on 25km comparisons, full block size fade distance works much better than half. Matches BS1 data well.


# What if using distfn- does this mean targets in block that are not in center cell will have less weight?
# either ignore or set center block source modifiers to 1
# OK- targets are already 1, modify only works on sources.

# what if sources and targs in diff components? does that mess up multiplier?
    # need to run components separately?
    # look at current flowing to ground to determine nsources? 
    # then divide current by this number, adjst nsources accordingly to get multiplier.

# what if center ground is in nodata?
# what if targets are not all in same component?

# write and read numpy grids?
# adddata for FA?
# higher random values for resisfa?


tile=-1 # For running overlapping tiles created by ArcGIS, to be stitched together later
#---------------------------------------------------------------------
# INPUTS #
#---------------------------------------------------------------------
options={}

# MOVING WINDOW AND TARGET BLOCK SIZES ####
options['radius'] = 56# In PIXELS. Search radius, with sources activated within the radius and outside of the center (target) block.
options['blockSize'] = 1 # Odd number. Targets will be square blocks of pixels with this number of pixels on a side.

options['projectDir'] = r'C:\Dropbox\Working\AdaptWest\PC_mockup'
options['resisRasterBase'] = 'HM_pow10_PNWmockup.tif'
options['outputDirBase'] = 'shortoutFolder2'

options['useSourceRaster']=False # Use a layer specifying source/target pixels to connect. 0= no source, anything positive will indicate strength of a source at that pixel. If not using a separate source raster, will use r Cutoff and consider anything with resistance lower than this to be a source/target
options['sourceRasterBase'] = '6x6Sources.asc'# Name of source raster, if using one
options['sourceInverseR']=False # Source strengths are inverse of input resistances (before squaring)

# RESISTANCES AND CUTOFF VALUES ####
options['rCutoff'] = 2000# If NOT using a source raster, everything <= this value in the ORIGINAL resistance raster will be a source (before squaring if squareResistances is True)
options['squareResistances']=True # Will square resistance raster before doing any calculations. sourceInverseR and rCutoff applied before squaring.

# CLIMATE ####
options['useClimate']=True
options['matchClimatePCs']=True  # climate method- principle components if True, present-day temperature differences if False
# parameters to match temperature differences
options['tDiff']=4 # Match pixels that differ in temperature by this value
options['tWindow']=1 #How close does target need to be to the tDiff value
options['climateRasterBase'] ='tmean_adks.asc'#'6x6t1pc2.asc'#'6x6clim.asc'#'TMEAN_NE_clip.tif'
options['absClim']=False # connect if ABS VAL of climate differs by tcutoff. Meant to help with fade.
# parameters to match current and future PCs
options['t1PC1RasterBase'] ='NORM_6190_PC1_PNWmockup.tif'
options['t1PC2RasterBase'] ='NORM_6190_PC2_PNWmockup.tif'
options['t2PC1RasterBase'] ='MIROC5_2080s_RCP85_PC1_PNWmockup.tif'
options['t2PC2RasterBase'] ='MIROC5_2080s_RCP85_PC2_PNWmockup.tif'
options['PCWindow'] = 0.9 # Euclidean distance between PCs to determine match

# DISTANCE FUNCTION#### 
options['useDistanceFunction'] = False # If true can set minimum and maximum distances or function between source and target pixels. Sources outside this range won't be activated.
options['distEq'] = "(options['radius']-dist)/options['radius']" # Distance equation for source strengths. Use NONE to ignore

options['minDist'] = None # In Pixels. Sources closer than this will have no current injected. Use None to ignore.
options['maxDist'] = None # Sources farther than this will have no current injected. Use None to ignore.
# Note: current can actually be higher in some pixels when min distances are used. Will be common with targetOnly mode because same amount
# of current is injected farther from target. But can also happen when opposing (canceling) currents occur without min dist.

# FLOW ACCUMULATION CALCULATIONS ####
options['calcFA']=False
options['addRandomResistances']=True #add random values to FA resis raster

# OPTION TO LIMIT ANALYSIS EXTENT ####  
# Bands are horizontal, with width equal to blockSize. There are approx nrows/blockSize bands in a raster. Stripes are vertical, with width equal to blocksize.
# These options allow you to only process a subset of bands and stripes. 
options['startBand'] = 0 # First band to process. Use 0 to ignore.
options['endBand'] =  0# Stops before processing this band. Use 0 to ignore.
options['startStripe'] = 0 #First stripe to process. Use 0 to ignore.
options['endStripe'] = 0 # Stops before processing this stripe. 0 to ignore.

# VOLTAGES #### Note: code not complete yet!
options['calcVoltages']=False
options['adjustVoltages']=False

# TARGETS AND WEIGHTING ####
options['weightTargOnly'] = True # Total current flow is ~ ntargs. May make sense if # dispersers limited, or number of dispersers a pixel can accept is limited. If false and noweight is false, current flow will be ~ntargs*nsources
# Recommend not changing following
options['noWeight']=False # Recommend False. 
options['centerGround']=True# Recommend True.
options['negTargets']= False # Recommend False. negative sources at targets- blocks can work better with centerground, fade out, and no neg targs.
options['subtractSources'] = False # Recommend False. 

# FADE CURRENTS #FIXME: check out divide by zero in fade calc
options['fadeIn']=True# Decreases current from each solve linearly with distance to center. Helps with artifacts from blocks. Best with center ground.
from math import *
options['fadeInDist']=(options['blockSize']-1)/2#/4.0#sqrt(2)*(options['blockSize']/2.0)#options['radius']/2 #
options['fadeConstant']=0.25 # Constant added to numerator and denominator in fade equation. Avoids zero current at center ground. Higher values mean more current at ground 

# Quantilize maps ####
options['quantilizeMaps'] = True # Additional current map will be written with values binned by quantile (e.g., 0-100 with 100 quantiles). Makes display easier (just use min-max stretch).
options['numQuantiles'] = 100

# SPECIAL FUNCTIONS
options['calcNull']=False # Calculates a 'null' result with all resistances = 1.
options['saveSourceAndTargetCounts'] = True # #Save source and target counts when doing climate analyses
options['doSourceAndTargetCountsOnly'] = True # if saving counts, DON'T calculate current
options['printTimings'] =  True
#---------------------------------------------------------------------
# END INPUTS #
#---------------------------------------------------------------------
from functools import wraps
PROF_DATA = {}
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
arcpy.env.overwriteOutput=True

def omniscape(options):
    theStart = datetime.datetime.now()  
    options = set_options_and_dirs(options)   
    copy_this_file(options)
    arcpy.env.scratchWorkspace = options['scratchDir']
    arcpy.env.Workspace = options['scratchDir']
    os.environ["TEMP"] = options['scratchDir']
    os.environ["TMP"] = options['scratchDir']    

    # Set raster paths and export to ascii if needed. FIXME: don't really need ascii, just convenient for header code for now
    resisRaster = path.join(options['projectDir'],options['resisRasterBase'])
    descData=arcpy.Describe(resisRaster)
    spatialReference=descData.spatialReference

    if options['useClimate']:
        if options['matchClimatePCs']:
            t1PC1Raster = path.join(options['projectDir'],options['t1PC1RasterBase'])
            t1PC2Raster = path.join(options['projectDir'],options['t1PC2RasterBase'])
            t2PC1Raster = path.join(options['projectDir'],options['t2PC1RasterBase'])
            t2PC2Raster = path.join(options['projectDir'],options['t2PC2RasterBase']) 
        else:    
            climateRaster = path.join(options['projectDir'],options['climateRasterBase'])
    
# FIXME TEMP            
    # resisRaster=arcpy.sa.Con(arcpy.sa.IsNull(arcpy.Raster(resisRaster)),1,1)
    # resisRaster.save('c:\\temp\\ones.tif')
    # resisRaster = 'c:\\temp\\ones.tif'
    
    if options['useSourceRaster']:
        sourceRaster = path.join(options['projectDir'],options['sourceRasterBase'])
    elif options['sourceInverseR']:
        sourceRaster = 1 / arcpy.Raster(resisRaster) 
    else:
        sourceRaster = arcpy.sa.Con((arcpy.Raster(resisRaster) <= options['rCutoff']),1,0)
   
    resisRaster = square_resistances(resisRaster,options)
    cumCurrentRaster = arcpy.sa.Con((arcpy.Raster(resisRaster) > 0),0,0)
    
    if options['calcFA']:    
        resisRasterFA = add_random_resistances(resisRaster,options)
  
    header = get_header(resisRaster)
      
    descData=arcpy.Describe(resisRaster)
    arcpy.env.extent=descData.Extent    
    
    if options['calcVoltages']:
        cumVdiffRaster  = arcpy.sa.Con((arcpy.Raster(resisRaster) > 0),0,0)   
    else: cumVdiffArray = cumVdiffRaster = None       
    if options['calcFA']:
        cumFlowRaster = arcpy.sa.Con((arcpy.Raster(resisRaster) > 0),0,0)
    else: cumFlowRaster = None
    if options['saveSourceAndTargetCounts'] and options['useClimate']: 
        cumTargetRaster = arcpy.sa.Con((arcpy.Raster(resisRaster) > 0),0,0)
        cumSourceRaster = arcpy.sa.Con((arcpy.Raster(resisRaster) > 0),0,0)
    else:
        cumTargetRaster = None
        cumSourceRaster = None
    
    iter = 0
    bandNum = 0
    start_time0 = time.clock()
    for centerRow in range((options['blockSize']-1)/2,header['nrows'],options['blockSize']):
        solveInBand=False
        bandNum += 1 
        if options['startBand'] > 0 and bandNum < options['startBand']: continue
        if options['endBand'] > 0 and bandNum >= options['endBand']: break
        if options['endBand'] > 0:
            print 'Starting band #',bandNum,' out of ',options['endBand'] - 1,' centered on row '+str(centerRow)
        else:
            print 'Starting band #',bandNum,' out of ',int(header['nrows']/options['blockSize']),' centered on row '+str(centerRow)

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
            
        cumCurrentArray = npy.zeros(bandArray.shape, dtype = 'float64') 
        if options['saveSourceAndTargetCounts'] and options['useClimate']:
            cumTargetArray = cumCurrentArray.copy() 
            cumSourceArray = cumCurrentArray.copy()
            
        if options['calcVoltages']: cumVdiffArray = cumCurrentArray.copy()
       
        subsetCenterRow = min(centerRow,options['radius'])
        
        # Check for all nodata in center band of resistance raster
        if options['blockSize']>1:
            bandCenterArray=bandArray[(subsetCenterRow-(options['blockSize']-1)/2):(subsetCenterRow+(options['blockSize']-1)/2),:]
            # if options['calcFA']:
                # bandCenterArrayFA=bandArrayFA[(subsetCenterRow-(options['blockSize']-1)/2):(subsetCenterRow+(options['blockSize']-1)/2),:]

        else:
            bandCenterArray=bandArray[subsetCenterRow,:]
            # if options['calcFA']:
                # bandCenterArrayFA=bandArrayFA[subsetCenterRow,:]            
        if npy.max(bandCenterArray) == -9999:       
            del bandCenterArray
            continue
        del bandCenterArray
        
        if options['blockSize']>1:
            sourceCenterArray=sourceBandArray[(subsetCenterRow-(options['blockSize']-1)/2):(subsetCenterRow+(options['blockSize']-1)/2+1),:]
        else:
            sourceCenterArray=sourceBandArray[subsetCenterRow,:]
            
        if npy.max(sourceCenterArray) <=0:            
            print 'no targets; continuing' #fixme: need to differently handle sources and targets
            del sourceCenterArray
            continue
        # del sourceCenterArray fixme: saving for now to quickly look for valid target areas 
        
            
        sourceCenterArraySum0 = npy.sum(npy.where(sourceCenterArray > 0, sourceCenterArray, 0), axis=0)
        
        grid = npy.indices((bandArray.shape))
        rowArray = grid[0]
        colArray = grid[1]
        del grid
        stripeNum=0
        for centerCol in range((options['blockSize']-1)/2, header['ncols'],options['blockSize']):
            iter += 1
            stripeNum += 1
            if options['startStripe'] > 0 and stripeNum < options['startStripe']: continue
            if options['endStripe'] > 0 and stripeNum >= options['endStripe']: break

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
            #xprint 'SOURCE',sourceCenterArraySum0
            circleResisArray = circ(bandArray, rowArray, colArray, subsetCenterRow, centerCol, options)
            if options['calcFA']:            
                circleResisArrayFA = circ(bandArrayFA, rowArray, colArray, subsetCenterRow, centerCol, options)
            
            sourceArray= circ(sourceBandArray, rowArray, colArray, subsetCenterRow, centerCol, options)
            sourceArray = npy.where(sourceArray < 0, 0, sourceArray)    #fixme- just set nodata sources to zero earlier    

            targetArray = center_block(sourceArray, options, subsetCenterRow, subsetCenterCol)                             
            # targetArray = npy.where(targetArray <0, 0,targetArray) 
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
            # del rowArray, colArray
            if sourceSum == 0 or targetSum==0:
                continue 
                    
            circleHeader = get_subset_header(sourceArray, header, options, centerRow, centerCol)
            yMin= circleHeader['yllcorner'] #fixme- check. coudl be something like: max(circleHeader['yllcorner'],circleHeader['yllcorner'] + ((circleHeader['nrows'] - centerRow - options['radius'] - 1) * circleHeader['cellsize']))
            if options['saveSourceAndTargetCounts'] and options['useClimate']:
                cumSourceArray = addData(cumSourceArray, sourceArray, subsetCenterRow, centerCol, options)
                cumTargetArray = addData(cumTargetArray, targetArray, subsetCenterRow, centerCol, options)
                if options['doSourceAndTargetCountsOnly']:
                    continue
            #then normalizing. 1 amp injected, 1 amp taken out
            sourceArray = sourceArray/(sourceSum+0.0)
            
            if options['negTargets']:
                targetArray = -targetArray/(targetSum+0.0) 
                sourceArray += targetArray
            print 'sourcesum, targsum',sourceSum,targetSum
            del targetArray
            
# EXPERIMENT2 #nope, still causes solver failures 
            # targetArray = targetArray*sourceSum
            # sourceArray = sourceArray*targetSum

# For FA, save as resource grid
            sourceAsciiFile = path.join(options['scratchDir'], 'source_r'+str(centerRow) + 'c' +str(centerCol)+'iter'+str(iter)+'.asc') 
            ascii_grid_writer(sourceAsciiFile, sourceArray, circleHeader, options)
            
# for FA, save as point
            groundAsciiFile = path.join(options['scratchDir'], 'ground_r'+str(centerRow) + 'c' +str(centerCol)+'iter'+str(iter)+'.asc')            
            ascii_grid_writer(groundAsciiFile, groundArray, circleHeader, options)
            # save_numpy(groundArray) #FIXME TEMP: shows that could save some time, not a huge amount.
            del groundArray 
            
            outputFN = 'target_r'+str(centerRow) + 'c' +str(centerCol)+'.out'
            outputFile = path.join(options['scratchDir'], outputFN)

# For FA, save as resistance grid
            resisAsciiFile = path.join(options['scratchDir'], 'resis_r'+str(centerRow) + 'c' +str(centerCol)+'iter'+str(iter)+'.asc') 

            csOptions = setCircuitscapeOptions(resisAsciiFile,sourceAsciiFile,groundAsciiFile,outputFile)

            if options['calcVoltages']:
                csOptions['write_volt_maps']=True
            
            if options['calcNull']:
                circleResisArray= npy.where(circleResisArray>0,1,circleResisArray)
                if options['calcFA']:
                    circleResisArrayFA= npy.where(circleResisArrayFA>0,1,circleResisArrayFA)
            ascii_grid_writer(resisAsciiFile, circleResisArray, circleHeader, options)

            print '\nDone with prep'
            start_time0=elapsed_time(start_time0)           

            if options['endBand'] > 0:
                print 'Ready to solve block row and col', centerRow, ',', centerCol, 'in stripe #',str(stripeNum),' and band#',bandNum,' out of ',options['endBand']-1
            else:
                print 'Ready to solve block row and col', centerRow, ',', centerCol, 'in stripe #',str(stripeNum),' and band#',bandNum,' out of ',int(header['nrows']/options['blockSize'])
            print 'Elapsed time so far: {0}'.format(datetime.datetime.now()-theStart)  
            
            if options['calcFA']:
                rFA2 = calc_fa(options,groundAsciiFile,circleHeader,yMin,sourceArray,circleResisArrayFA,iter)
            solveInBand = True
            curMapPath, outConfigFile = calc_current(options, csOptions, centerRow, centerCol) 
            if not path.exists(curMapPath):
                print "Can't find circuitscape output"
                exit(0)
                continue
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
            print sourceSum, targetSum, multiplier
            
            currentArray = ascii_grid_reader(curMapPath, header['nodata'], 'float64')
            
            if options['centerGround']:
                maxCur = currentArray[subsetCenterRow, subsetCenterCol]
            else:
                maxCur = npy.max(currentArray)
            if maxCur <= 0:
                print 'NO CURRENT, continuing'
                continue

            if options['subtractSources']:
                currentArray = multiplier * (currentArray - abs(sourceArray))
            else:
                currentArray = multiplier * currentArray
            #xprint 'GROUND CUR',ascii_grid_reader(curMapPath, header['nodata'], 'float64')[subsetCenterRow, subsetCenterCol]
            del sourceArray
            if options['calcFA']:
                rFA3 = arcpy.sa.Times(rFA2, int(multiplier)) #fixme: may be faster to convert rfa or rfa2 into array and do numpy calcs for entire band, then add band in
                centralityFile= path.join(options['scratchDir'],"FA_iter" + str(iter) + '.tif')
                # rFA3.save (centralityFile)    # fixme: not needed really 
                cumFlowRaster = addData_arcpy(rFA3, cumFlowRaster)
                
#fixme- may not be compatible with subsources, weights, polygons, etc:
#fixme: dont' mult by dists, but dist/options['radius']?
            
            if options['fadeIn']:
#temp    
                # subsetHeader = get_subset_header(currentArray, header, options, centerRow, centerCol)
                # yMin= max(subsetHeader['yllcorner'],subsetHeader['yllcorner'] + ((subsetHeader['nrows'] - subsetCenterRow - options['radius'] - 1) * subsetHeader['cellsize']))
                # LLC = arcpy.Point(subsetHeader['xllcorner'],yMin)

                # currentRaster = arcpy.NumPyArrayToRaster(currentArray,LLC, subsetHeader['cellsize'],subsetHeader['cellsize'],-9999)
                # currentRaster.save('c:\\temp\\curRas_iter'+str(iter) + '.tif')
#temp

                currentArray = fade_currents(currentArray, subsetCenterRow, subsetCenterCol, options)  
                
#FIXME- drop fadearray above
                
                # currentRaster = arcpy.NumPyArrayToRaster(currentArray,LLC, subsetHeader['cellsize'],subsetHeader['cellsize'],-9999)
                # currentRaster.save('c:\\temp\\curRasFADE_iter'+str(iter) + '.tif')
                # fadeRaster = arcpy.NumPyArrayToRaster(fadeArray,LLC, subsetHeader['cellsize'],subsetHeader['cellsize'],-9999)
                # fadeRaster.save('c:\\temp\\fadeRas_iter'+str(iter) + '.tif')

#put in fn once striping solved below
            if options['calcVoltages']:
                voltMapPath = path.join(options['scratchDir'], 'target_r'+str(centerRow) + 'c' +str(centerCol) + '_voltmap.asc')
                voltMap = (ascii_grid_reader(voltMapPath, header['nodata'], 'float64'))
                voltMap = npy.where(circleResisArray<0,npy.nan,voltMap)
                #xprint 'voltmap'
                #xprint voltMap
                
                max_vdiff = get_max_vdiff(voltMap,circleResisArray,csOptions)
                max_vdiff=npy.where(npy.isnan(max_vdiff),0,max_vdiff)
                #xprint 'vidiff'
                #xprint max_vdiff
                # blarg
# fixme temp fix for striping at edges
                max_vdiff[0,::]=0
                max_vdiff[:,0]=0
                max_vdiff[:,circleHeader['ncols']-1]=0
                max_vdiff[circleHeader['nrows']-1,:]=0
                # vDiffPath = path.join(options['scratchDir'], 'target_r'+str(centerRow) + 'c' +str(centerCol) + '_VDIFF.asc')
                # ascii_grid_writer(vDiffPath, max_vdiff, circleHeader, options)

                if options['fadeIn']: #Fixme: do we want to fade voltages?
                    max_vdiff = npy.where(max_vdiff>0,npy.multiply(max_vdiff,fadeArray),max_vdiff)
                    del fadeArray

                if options['weightTargOnly']:
                    max_vdiff = targetSum * max_vdiff
                elif options['noWeight']:
                    pass
                else:    
                    max_vdiff = sourceSum * targetSum * max_vdiff
                cumVdiffArray = addData(cumVdiffArray, max_vdiff, subsetCenterRow, centerCol, options)
           
                delete_data(voltMapPath)
            del circleResisArray    
            delete_data(curMapPath)
            delete_data(outConfigFile)
            delete_data(groundAsciiFile)
            delete_data(resisAsciiFile)
            delete_data(sourceAsciiFile)
            
            start_time1 = time.clock()
            cumCurrentArray = addData(cumCurrentArray, currentArray, subsetCenterRow, centerCol, options)
            del currentArray

# FIXME: do all raster stuff for each band, not each solve
            # cumCurrentRaster = addData_arcpy(cumCurrentRaster, currentRaster)
# FIXME: move file writing to band iteration


        del rowArray,colArray
        del sourceCenterArray,bandArray,sourceBandArray
        if options['calcFA']:
            del bandArrayFA
            
        if options['saveSourceAndTargetCounts']:
            yMin= max(header['yllcorner'],header['yllcorner'] + ((header['nrows'] - centerRow - options['radius'] - 1) * header['cellsize']))
            LLC = arcpy.Point(header['xllcorner'],yMin)
            bandSourceRaster = arcpy.NumPyArrayToRaster(cumSourceArray,LLC, header['cellsize'],header['cellsize'],-9999)
            cumSourceRaster = addData_arcpy(cumSourceRaster, bandSourceRaster)       
        
            del cumSourceArray, bandSourceRaster
            bandTargetRaster = arcpy.NumPyArrayToRaster(cumTargetArray,LLC, header['cellsize'],header['cellsize'],-9999)
            cumTargetRaster = addData_arcpy(cumTargetRaster, bandTargetRaster)       
            del cumTargetArray, bandTargetRaster
        
        if solveInBand: #
            # bandRows=min(options['radius']*2+1,options['radius']+centerRow+1)
            yMin= max(header['yllcorner'],header['yllcorner'] + ((header['nrows'] - centerRow - options['radius'] - 1) * header['cellsize']))
            LLC = arcpy.Point(header['xllcorner'],yMin)

            bandCurrentRaster = arcpy.NumPyArrayToRaster(cumCurrentArray,LLC, header['cellsize'],header['cellsize'],-9999)
                       
            cumCurrentRaster = addData_arcpy(cumCurrentRaster, bandCurrentRaster)
            del bandCurrentRaster,cumCurrentArray
            if options['calcVoltages']:
                bandVdiffRaster = arcpy.NumPyArrayToRaster(cumVdiffArray,LLC,header['cellsize'],header['cellsize'],-9999)
                cumVdiffRaster = addData_arcpy(cumVdiffRaster, bandVdiffRaster)
                del cumVdiffArray, bandVdiffRaster
                
            options = write_temp_maps(options,bandNum,cumCurrentRaster,cumVdiffRaster)
            
            
            #fixme: move to write_temp_maps, add final one too.
            if options['calcFA']:
                cumFlowFile = os.path.join(options['outputDir'], 'BAND'+str(bandNum)+'flow_' + options['outputFileText']+'.tif')
                cumFlowRasterInt = arcpy.sa.Int(cumFlowRaster)
                cumFlowRasterInt.save(cumFlowFile)
                delete_data(options['prevFlowFile'])
                options['prevFlowFile'] = cumFlowFile
            
        print 'Done with band #',bandNum,'.Elapsed time so far: {0}'.format(datetime.datetime.now()-theStart)  

    print 'Done with solves.'
    
    
    print_prof_data()

    # write_final_maps(options,cumCurrentRaster,cumVdiffRaster,cumFlowRaster,cumSourceRaster,cumTargetRaster) 

    write_final_maps(options,cumCurrentRaster,cumVdiffRaster,cumFlowRaster,cumSourceRaster,cumTargetRaster) 

    
    
    
    
    
#fixme temp    
    clean_up(options)
    #xprint locals()

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

# To evaluate speed differences. definitely faster than saving and reading asciis. may speed up cs too?    
# @profile
# def save_numpy(array):
    # npyFile=r'c:\temp\temp.npy'
    # npy.save(npyFile, array)     
    # npyArray = npy.load(npyFile, mmap_mode=None)

@profile
def calc_current(options, csOptions, centerRow, centerCol): 
    configFN = 'config_target_b'+str(centerRow) + 'c' +str(centerCol)+'.ini'
    outConfigFile = path.join(options['scratchDir'], configFN)
    writeCircuitscapeConfigFile(outConfigFile, csOptions)  
    solveInBand=True
    CSPATH = get_cs_path()                            
    start_time1 = time.clock()
    call_circuitscape(CSPATH, outConfigFile)
    start_time1 = elapsed_time(start_time1)
    start_time0 = time.clock()

    delete_data(configFN)
    configFN2 = 'target_r'+str(centerRow) + 'c' +str(centerCol)+'.ini' #fixme?
    delete_data(path.join(options['scratchDir'], configFN2))

    curMapPath = path.join(options['scratchDir'], 'target_r'+str(centerRow) + 'c' +str(centerCol) + '_curmap.asc')
    return curMapPath, outConfigFile

@profile
def calc_fa(options,groundAsciiFile,circleHeader,yMin,sourceArray,circleResisArrayFA, iter):
    # fixme: convoluted way to get a shapefile via raster conversion. simplify
    # FA point
    outPoint = path.join(options['scratchDir'], 'point_iter'+str(iter)+'.shp')
    if arcpy.Exists(outPoint):
        arcpy.Delete_management(outPoint)
    field = "VALUE"
    arcpy.RasterToPoint_conversion(groundAsciiFile, outPoint, field)    
    # FA source strength
    # yMin= circleHeader['yllcorner'] #fixme- check. coudl be something like: max(circleHeader['yllcorner'],circleHeader['yllcorner'] + ((circleHeader['nrows'] - centerRow - options['radius'] - 1) * circleHeader['cellsize']))
    LLC = arcpy.Point(circleHeader['xllcorner'],yMin)
    sourceRasFA = arcpy.NumPyArrayToRaster(sourceArray,LLC, circleHeader['cellsize'],circleHeader['cellsize'],-9999)    
    # sourceRasFA.save(path.join(options['scratchDir'], 'sourceRasFA_iter'+str(iter)+'.tif')) #fixme: not needed
    resisRasFA = arcpy.NumPyArrayToRaster(circleResisArrayFA,LLC, circleHeader['cellsize'],circleHeader['cellsize'],-9999)    
    del circleResisArrayFA
    # resisRasFA.save('c:\\temp\\xresisRasFA'+'iter'+str(iter)+'.tif') #fixme not needed
    rOut = 'FA'
    
    print '\nStarting flow accumulation'
    start_timeFA = time.clock()
    rCB = arcpy.sa.CostBackLink(outPoint, resisRasFA, '#')#, rOut + "d" + str(FID) ) #BHM can limit cwd here 999999999
    arcpy.Delete_management(outPoint)
    # rCB.save('c:\\temp\\xxrCB'+'iter'+str(iter)+'.tif') #fixme not needed
    rFD = arcpy.sa.Con( rCB > 0, arcpy.sa.Power(2, (rCB - 1))) 
    # rFD.save('c:\\temp\\xxrFD'+'iter'+str(iter)+'.tif') #fixme not needed
    
    rFA = arcpy.sa.Plus(arcpy.sa.FlowAccumulation( rFD, sourceRasFA ), sourceRasFA) #Note! added in source strengths 8/15/15 
    
    # rFA.save('c:\\temp\\xxrFA'+'iter'+str(iter)+'.tif') #fixme not needed
    rFA2=arcpy.sa.Con(arcpy.sa.IsNull(rFA),0,rFA)
    start_timeFA = elapsed_time(start_timeFA)
    return rFA2
    

@profile
def fade_currents(currentArray, subsetCenterRow, subsetCenterCol, options):
    grid = npy.indices((currentArray.shape))
    subsetRowArray = grid[0]
    subsetColArray = grid[1]
    del grid
    subsetDistArray = npy.sqrt(npy.multiply(subsetRowArray - subsetCenterRow, subsetRowArray- subsetCenterRow) + npy.multiply(subsetColArray-subsetCenterCol, subsetColArray-subsetCenterCol))           
    fadeArray=npy.where(subsetDistArray<options['fadeInDist'],(subsetDistArray+options['fadeConstant'])/(options['fadeInDist']+options['fadeConstant']),1)
    currentArray=npy.where(currentArray>0,npy.multiply(currentArray,fadeArray),currentArray)
    return currentArray 
    
    
def square_resistances(resisRaster,options):
        if options['squareResistances']: #fixme:put in fn
            fileBase, fileExtension = path.splitext(options['resisRasterBase'])
            resisRasterSQ = path.join(options['scratchDir'],'sq_'+fileBase+'.tif')
            outResistanceRaster = arcpy.sa.Times(resisRaster, resisRaster) 
            outResistanceRaster.save(resisRasterSQ)           
            resisRaster = resisRasterSQ    
        return resisRaster

@profile        
def match_climate_pcs(sourceArray, targetArray, t1PC1BandArray, t1PC2BandArray, t2PC1BandArray, t2PC2BandArray, rowArray, colArray, subsetCenterRow, centerCol, options):                
    print 'Starting climate using PCs'
    start_timeClimate = time.clock()
    t1PC1Array= circ(t1PC1BandArray, rowArray, colArray, subsetCenterRow, centerCol, options)              
    t1PC2Array= circ(t1PC2BandArray, rowArray, colArray, subsetCenterRow, centerCol, options)              
    t2PC1Array= circ(t2PC1BandArray, rowArray, colArray, subsetCenterRow, centerCol, options)              
    t2PC2Array= circ(t2PC2BandArray, rowArray, colArray, subsetCenterRow, centerCol, options)              
    sourceArrayNew=npy.zeros(sourceArray.shape,dtype='float64')
    # targetArrayNew=npy.zeros(sourceArray.shape,dtype='float64')
    
    # indices of valid targets
    tRows, tCols=npy.where(targetArray)
    targetTally=0

    for i in range(0,len(tRows)):
        t2PC1Target=t2PC1Array[tRows[i],tCols[i]]     #Targets are future. Fixme: would need to calculate t1pcs for targets if doing absclim analogy
        t2PC2Target=t2PC2Array[tRows[i],tCols[i]] 
        if t2PC1Target==-9999 or t2PC2Target == -9999: #fixme: do nans
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
        PCCutoffArray=npy.where((PCDistArray>= -options['PCWindow']) & (PCDistArray<=options['PCWindow']),1,0) 
        if len(sourceArray) < 10:
            print 'PCCutoffArray'
            print PCCutoffArray
        sourceArrayTarget_i=npy.multiply(sourceArray,PCCutoffArray)# sources for this target cell
        if npy.max(sourceArrayTarget_i)>0:
            targetTally+=1 # There's a valid source for this target. Increment target count for scaling current in options['weightTargOnly'] setting
        
        # For climate, may have more or fewer targets. Need to scale sources by ntargs, and targets by nsources
        sourceArrayNew=sourceArrayNew+sourceArrayTarget_i
        
#?/ not sure how to handle                                        
#This is where sourcesum and targetsum become the same                    
# putting in current equal to nsource*ntarg, taking same out. then normalizing to unit current. then mult by nsource*ntarg.
        #taking same out. 
        targetArray[tRows[i],tCols[i]] =sourceArrayTarget_i.sum()                   
    del tRows, tCols
           
    # take sum of sourceArrayTarget, set target strength to that...
    if len(sourceArray) < 10:
        print 'source new'
        print sourceArrayNew
        print 'target new'
        print targetArray
        
    sourceArray=sourceArrayNew
    del sourceArrayNew
    sourceSum = sourceArray.sum()
    targetSum=targetArray.sum() 
    #xprint 'targetTally',targetTally
    # Note that we now have identical source sums and targetsums
    start_timeClimate = elapsed_time(start_timeClimate) #Fixme: can climate be sped up? Taking up to 1 second with blocksize 25, radius 100           

    return sourceArray, targetArray, sourceSum, targetSum, targetTally

@profile
def match_temperature_diffs(sourceArray, targetArray, climateBandArray, rowArray, colArray, subsetCenterRow, centerCol, options):                
    # Fixme: need to add a quick search for overall potential for matches within bands
    print 'Starting climate using temperature differential'
    start_timeClimate = time.clock()
    climateArray= circ(climateBandArray, rowArray, colArray, subsetCenterRow, centerCol, options)              
    sourceArrayNew=npy.zeros(sourceArray.shape,dtype='float64')
    
    # indices of valid targets
    tRows, tCols=npy.where(targetArray)
    if len(sourceArray) < 10:
        print 'clim'
        print climateArray
    targetTally=0
    for i in range(0,len(tRows)):
        tTarget=climateArray[tRows[i],tCols[i]] 
        if tTarget==-9999:
            targetArray[tRows[i],tCols[i]] = 0
            continue
        tDiffArray=npy.where(climateArray == -9999, 0, climateArray-tTarget) #target is cooler. # fixme- could save a bit of time by removing where- maybe could just mask outinvalid- or uses nans, or remove -999 values later
       
        if options['absClim']:
            tCutoffArray=npy.where((abs(tDiffArray)>=options['tDiff']-options['tWindow']) & (abs(tDiffArray)<=options['tDiff']+options['tWindow']),1,0)
        else:
            tCutoffArray=npy.where((tDiffArray>=options['tDiff']-options['tWindow']) & (tDiffArray<=options['tDiff']+options['tWindow']),1,0) 
        if len(sourceArray) < 10:
            print 'tCutoffArray'
            print tCutoffArray
        sourceArrayTarget_i=npy.multiply(sourceArray,tCutoffArray)# sources for this target cell
        if npy.max(sourceArrayTarget_i)>0:
            targetTally+=1 # There's a valid source for this target. Increment target count for scaling current in options['weightTargOnly'] setting
        
        # For climate, may have more or fewer targets. Need to scale sources by ntargs, and targets by nsources
        sourceArrayNew=sourceArrayNew+sourceArrayTarget_i
        
#?/ not sure how to handle                                        
#This is where sourcesum and targetsum become the same                    
# putting in current equal to nsource*ntarg, taking same out. then normalizing to unit current. then mult by nsource*ntarg.
        #taking same out. 
        targetArray[tRows[i],tCols[i]] =sourceArrayTarget_i.sum()                   
    del tRows, tCols
        
    # take sum of sourceArrayTarget, set target strength to that...
    if len(sourceArray) < 10:
        print 'source new'
        print sourceArrayNew
        print 'target new'
        print targetArray
    sourceArray=sourceArrayNew
    del sourceArrayNew
    sourceSum = sourceArray.sum()
    targetSum=targetArray.sum()
    print 'targetTally',targetTally
    # Note that we now have identical source sums and targetsums
    start_timeClimate = elapsed_time(start_timeClimate) #Fixme: can climate be sped up? Taking up to 1 second with blocksize 25, radius 100           
    return sourceArray, targetArray, sourceSum, targetSum, targetTally
    
@profile
def modify_source_strengths_by_distance(sourceArray, subsetCenterRow, subsetCenterCol, options):
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
            distanceFunctionArray = npy.ones(sourceArray.shape,dtype='float64') #FIXME! Need distance function. Right now just doing min and max. 

        #Array to modify source array- includes distance function, max and min distances
        if not options['maxDist']:
            distanceModifierArray=npy.where(subsetDistArray>=options['minDist'],distanceFunctionArray,0)
        elif not options['minDist']:
            distanceModifierArray=npy.where(subsetDistArray<=options['maxDist'],distanceFunctionArray,0)
        else:
            distanceModifierArray=npy.where((subsetDistArray<=options['maxDist']) & (subsetDistArray>=options['minDist']),distanceFunctionArray,0)
   
        sourceArrayNew=npy.where(sourceArray>0,npy.multiply(sourceArray,distanceModifierArray),sourceArray)
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
    
    
def add_random_resistances(resisRaster, options):
    if options['addRandomResistances']: 
        fileBase, fileExtension = path.splitext(options['resisRasterBase'])
        resisRasterFA = path.join(options['scratchDir'],'rand_'+fileBase+'.tif')
        descData=arcpy.Describe(resisRaster)
        extent=descData.Extent
        cellSize=descData.MeanCellHeight
        seedValue = 1

        outRandomRaster = arcpy.sa.CreateRandomRaster(seedValue, cellSize, extent) 
        outResisRaster = arcpy.sa.Plus(outRandomRaster, resisRaster) 
        outResisRaster.save(resisRasterFA)
        # resisRaster=resisRasterFA #Fixme: can this just be in memory isntead?
    return resisRasterFA
    
def quantilize(raster):
    try:
        if raster is not None:       
            inArray = arcpy.RasterToNumPyArray(raster,"#","#","#",-9999)
            interval = 100.0/options['numQuantiles']
            breakSequence = seq(interval, 100-interval, interval)
            quantileArray = npy.zeros_like(inArray, dtype='int32')        
            quantileArray = npy.where(inArray > 0,options['numQuantiles'],quantileArray)
                       
            inVector = inArray.flatten()
            ind=npy.where(inVector > 0) 
            inVector = inVector[ind] #Eliminate nodata and zeros

            # quantilize                      
            if len(inVector)==0:
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
            
            descData=arcpy.Describe(raster)
            cellSize=descData.meanCellHeight
            extent=descData.Extent
            spatialReference=descData.spatialReference
            
            pnt=arcpy.Point(extent.XMin,extent.YMin)
            quantileRaster = arcpy.NumPyArrayToRaster(quantileArray,pnt,
                                                 cellSize,cellSize,-9999)
            arcpy.DefineProjection_management(quantileRaster,spatialReference)

            del quantileArray
            return quantileRaster
        else: return None
    except: 
        print 'Failed to quantilize'
        return None
      

        
def seq(start, stop, step=1):
    n = int(round((stop - start)/float(step)))
    if n > 1:
        return([start + step*i for i in range(n+1)])

def delete_data(file):
    try:
        if os.path.isfile(file):
            os.remove(file)
            gc.collect()
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
    print 'Creating final output maps'  
    quantCurrentRaster = quantVdiffRaster = quantFlowRaster = None
    if options['quantilizeMaps']:
        quantCurrentRaster  = quantilize(cumCurrentRaster)
        # quantVdiffRaster = quantilize(cumVdiffRaster) 
        # quantFlowRaster = quantilize(cumFlowRaster)

    sumFile = path.join(options['outputDir'], 'cur_' + options['outputFileText']+'.tif')
    print 'Writing:\n',sumFile
    outRasterString = sumFile
    cumCurrentRaster.save(sumFile)
    del cumCurrentRaster
    delete_data(options['prevCumCurrentFile'])
    if quantCurrentRaster is not None:
        sumFile = path.join(options['outputDir'], 'cur_PCT_' + options['outputFileText']+'.tif')
        outRasterString = outRasterString + '; ' + sumFile
        print 'Writing:\n',sumFile
        quantCurrentRaster.save(sumFile)
        del quantCurrentRaster
    
    if options['calcFA']:       
        sumFile = path.join(options['outputDir'], 'flow_' + options['outputFileText']+'.tif')
        outRasterString = outRasterString + '; ' + sumFile
        cumFlowRaster2 = arcpy.sa.Int(cumFlowRaster)
        cumFlowRaster3 = arcpy.sa.Con(cumFlowRaster2 > 0,cumFlowRaster2)#(arcpy.Raster(cumFlowRaster))
        print 'Writing:\n',sumFile
        cumFlowRaster3.save(sumFile)
        
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
    
    
@profile    
def write_temp_maps(options,bandNum,cumCurrentRaster,cumVdiffRaster): 
    print 'Writing temporary grids...'
    cumCurrentFile = os.path.join(options['outputDir'], 'BAND'+str(bandNum)+'cur_' + options['outputFileText']+'.tif')
    try:
        cumCurrentRaster.save(cumCurrentFile)
    except:
        try:
            print 'Error writing'
            time.sleep(5)                
            cumCurrentFile = os.path.join(options['outputDir'],'BAND'+str(bandNum)+options['resisRasterBase']+'_sumtemp2.tif')
            cumCurrentRaster.save(cumCurrentFile)

        except:
            print 'Second error writing- may need to reinstate creation of new scratch directory if fails again'
            time.sleep(5)
            cumCurrentFile = os.path.join(options['outputDir'],'BAND'+str(bandNum)+options['resisRasterBase']+'_sumtemp3.tif')
            cumCurrentRaster.save(cumCurrentFile)                       
            
    delete_data(options['prevCumCurrentFile'])
    options['prevCumCurrentFile'] = cumCurrentFile

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
        delete_data(options['prevVdiffFile'])
        options['prevVdiffFile'] = cumVdiffFile
    return options
    
            
def clean_up(options):
    print 'Cleaning up...'
    # try:
        # delete_dir(options['scratchDir']) #can't figure out why but this causes a python crash. Works if in main function and not here.
    # except:
        # pass
    for filename in glob.glob(os.path.join(options['outputDir'],'band*.*')) :
        try:    
            os.remove( filename )
        except:
            pass
        
        
def set_options_and_dirs(options):
    """derives output base filename by combining input options. Also handles filenames for tiling    """    
    if options['calcNull']:
        options['outputDirBase'] = options['outputDirBase']+'_NULL'
    options['projectDirBase'] = 'scratch'+options['outputDirBase']
    options['compress'] = False
    
    if tile >=0: 
        options['outputDirBase']=options['outputDirBase']+str(tile)
        options['projectDirBase']=options['projectDirBase']+str(tile)
        fileBase,ext=os.path.splitext(options['resisRasterBase'])
        options['resisRasterBase']=fileBase+str(tile)+ext
        if options['useClimate']:
            if not options['matchClimatePCs']:
                fileBase,ext=os.path.splitext(options['climateRasterBase'])
                options['climateRasterBase']=fileBase+str(tile)+ext
            else:
                fileBase,ext=os.path.splitext(options['t1PC1RasterBase'])
                options['t1PC1RasterBase']=fileBase+str(tile)+ext
        
                fileBase,ext=os.path.splitext(options['t2PC1RasterBase'])
                options['t2PC1RasterBase']=fileBase+str(tile)+ext

                fileBase,ext=os.path.splitext(options['t1PC2RasterBase'])
                options['t1PC2RasterBase']=fileBase+str(tile)+ext

                fileBase,ext=os.path.splitext(options['t2PC2RasterBase'])
                options['t2PC2RasterBase']=fileBase+str(tile)+ext
                
        if options['useSourceRaster']:
            fileBase,ext=os.path.splitext(options['sourceRasterBase'])    
            options['sourceRasterBase']=fileBase+str(tile)+ext
            
    if options['startBand'] is None: options['startBand'] = 0 
    if options['startStripe'] is None: options['startStripe'] = 0
    if options['endBand'] is None: options['endBand'] = 0
    if options['endStripe'] is None: options['endStripe'] = 0
    
    if options['radius'] <= options['blockSize']/2:
        print 'Error. Radius must be larger than half the block size.'
        exit(0)

    if float(options['blockSize'])/2 == int(options['blockSize']/2):
        print float(options['blockSize'])/2
        print int(options['blockSize']/2)
        print 'Error. Block size must be an odd number.'
        exit(0)
        
    resisRasText, fileExtension = os.path.splitext(options['resisRasterBase'])
    if not options['useSourceRaster']:
        if len(resisRasText)>25:
            resisRasText=resisRasText[0:24]+'x'                    
    elif len(resisRasText)>15:
        resisRasText=resisRasText[0:14]+'x'                    

    if options['squareResistances']:
        squareText='_sq'
    else:    
        squareText=''
    radiusText = '_r'+str(options['radius'])
    if options['useSourceRaster']:
        fileBase,ext=os.path.splitext(options['sourceRasterBase'])
        srcText = 'srcRas_'+ fileBase + '_'
    else:
        srcText = ''
    if options['negTargets']:
        negTargText = 'negTarg'
    else:
        negTargText = ''

    if options['fadeIn']:
        fadeInText='_fad'+str(int(options['fadeInDist']))
    else:
        fadeInText=''

    if options['startBand'] > 0:
        startBandText = '_SB'+str(options['startBand'])
    else:
        startBandText=''
    if options['endBand'] > 0:
        endBandText = '_EB'+str(options['endBand'])
    else:
        endBandText=''
    if options['startStripe'] > 0:
        startStripeText = '_SS'+str(options['startStripe'])
    else:
        startStripeText=''
    if options['endStripe'] > 0:
        endStripeText = '_ES'+str(options['endStripe'])
    else:
        endStripeText=''
        
    distFunctionText=''
    if options['useDistanceFunction'] and (options['maxDist'] or options['minDist'] or options['distEq'] is not None):
        if options['minDist']:
            distFunctionText= '_minDist'+str(options['minDist']) 
        if options['maxDist']:
            distFunctionText= distFunctionText + '_maxDist'+str(options['maxDist']) 
        if options['distEq'] is not None:
            distFunctionText= distFunctionText + '_distFn'
    if options['weightTargOnly']:
        weightText='_targOnly'
    else:
        weightText=''
    if options['centerGround']:
        centerText='cg'
    else:
        centerText=''
    

    if options['noWeight']:
        noWeightText='noWeight'
    else:
        noWeightText=''
    if options['subtractSources']:
        subtrText = 'subSrc'            
    else:
        subtrText = ''

    if options['useClimate']:
        if not options['matchClimatePCs']:

            fileBase,ext=os.path.splitext(options['climateRasterBase'])
            if len(fileBase)>15:
                fileBase=fileBase[0:14]+'x'
            climText='clim_tc_'+str(options['tDiff']).replace('.','')+'_'+fileBase
            if options['absClim']:
                climText=climText+'_abs'
        else:
            climText='_clim_PCs_Cut'+str(options['PCWindow']).replace('.','')
    else:
        climText=''
    if options['calcNull']:
        nullText='NULL_'
    else:
        nullText=''

    if not options['useSourceRaster']:
        if options['sourceInverseR']:
            cutoffText = '_srcInvR'
        else:    
            cutoffText ='_rc'+str(options['rCutoff'])
    else:
        cutoffText = ''
    options['outputFileText'] = nullText+resisRasText + squareText + '_'+srcText+climText+radiusText+'bl'+str(options['blockSize'])+cutoffText +distFunctionText+weightText+noWeightText+subtrText+centerText+startBandText+endBandText+startStripeText+endStripeText+negTargText+fadeInText
    options['prevFlowFile']=None
    options['prevCumCurrentFile']=None
    options['prevVdiffFile'] = None

    print 'resisRasterBase=',options['resisRasterBase']
    print 'Radius=',options['radius']
    print 'BlockSize=',options['blockSize']
    if not options['useSourceRaster']:
        if options['sourceInverseR']:
            print 'Using inverse of resistances for sources'
        else:
            print'rCutoff=',options['rCutoff']
    else:
        print 'sourceRasterBase=',options['sourceRasterBase']
    print'startBand',options['startBand']
    print'endBand',options['endBand']
    print'startStripe',options['startStripe']
    print'endStripe',options['endStripe']
    if options['calcNull']:
        print'calcNull',options['calcNull']
    if options['fadeIn']:
        print 'fadeIn Distance',options['fadeInDist']
    if options['useDistanceFunction']:
        print 'minDist',options['minDist'] 
        print 'maxDist',options['maxDist'] 

    options['outputDir']=os.path.join(options['projectDir'],options['outputDirBase'])
    options['scratchDir']=os.path.join(options['projectDir'],options['projectDirBase'])
    print 'project Dir',options['projectDir']
    print 'output Dir',options['outputDirBase']
    if not path.exists(options['projectDir']):
        os.mkdir(options['projectDir'])
    if not path.exists(options['scratchDir']):
        os.mkdir(options['scratchDir'])
    if not path.exists(options['outputDir']):
        os.mkdir(options['outputDir'])
    print '\n'
    return options

def copy_this_file(options):
    # Save a copy of this file in output directory
    destFile=os.path.join(options['outputDir'],'omniscape_'+options['outputFileText']+'.py')
    ft = tuple(time.localtime())
    timeNow = time.ctime()
    fileName = ('%s_%s_%s_%s%s_%s' % (ft[0], ft[1], ft[2], ft[3], ft[4], os.path.basename(sys.argv[0])))
    filePath = os.path.join(options['outputDir'],fileName)
    shutil.copyfile(sys.argv[0],filePath) 
    shutil.copyfile(sys.argv[0],destFile) 

            
def get_max_vdiff(voltMap,circleResisArray,csOptions):
    voltMap_l, voltMap_r = get_horiz_neighbors(voltMap)
    voltMap_u, voltMap_d = get_vert_neighbors(voltMap)
    print 'voltMap_l'
    print voltMap_l
    if options['adjustVoltages']: #fixme: not completed. Aim is to 
        # see how much improvement is possible if r reduced to 1.
        # vdiff * (r-1)/r... if r=1, no improvement possible. if r=2, half of voltage could be reduced. etc.
        # need r's, will include r from 3 pixels. but try to calc how much v would drop if that one pixel were restored.
        # HOW MUCH WOULD V DROP ACROSS PIXEL IF R>1
        
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
    
def get_vdiff(voltMap,voltMap_direction):
#FIXME: vdiff code not working right, debug using 6x6
    vdiff_direction = npy.absolute(voltMap - voltMap_direction)
    vdiff_direction = npy.where(voltMap == npy.nan, 0, npy.where(voltMap_direction == npy.nan,0,vdiff_direction))
    return vdiff_direction

def get_horiz_neighbors(map):
    m = map.shape[0]
    n = map.shape[1]
    
    # zeromap = npy.zeros(map.shape,dtype = 'float64') 
    zeromap = npy.empty(map.shape,dtype = 'float64') 

    zeromap[:] = npy.NAN
    map_l  = zeromap.copy()
    map_r  = zeromap.copy()
    
    map_l[:,1:n] = map[:, 0:(n-1)]
    map_r[:,0:(n-1)] = map[:, 1:n]
    return map_l, map_r

def get_vert_neighbors(map):
    m = map.shape[0]
    n = map.shape[1]
    zeromap = npy.zeros(map.shape,dtype = 'float64') 
    map_u  = zeromap.copy()
    map_d  = zeromap.copy()

    map_u[1:m, :] = map[0:(m-1), :]
    map_d[0:(m-1) , :] = map[1:m , :]
    
    return map_u, map_d

def get_diag1_neighbors(map):
    m = map.shape[0]
    n = map.shape[1]   
    zeromap = npy.zeros(map.shape,dtype = 'float64') 
    map_ul  = zeromap.copy()
    map_ul[1:m,1:n] = map[0:m-1, 0:n-1]
    map_dr  = zeromap.copy()
    map_dr[0:m-1, 0:n-1  ] = map[1:m , 1:n ]
    return map_ul, map_dr


def get_diag2_neighbors(map):
    m = map.shape[0]
    n = map.shape[1]
    zeromap = npy.zeros(map.shape,dtype = 'float64') 
    map_ur  = zeromap.copy()
    map_ur[1:m,0:n-1] = map[0:m-1, 1:n  ]
    map_dl  = zeromap.copy()
    map_dl[0:m-1, 1:n  ] = map[1:m  , 0:n-1]
    return map_ur, map_dl

# @profile
def band(inRaster,header,centerRow, options):
#FIXME put bounds on bandrows- minrow=0, maxrow=nrows...
    # bandRows=min(options['radius']*2+1,options['radius']+centerRow+1,header['nrows'])
   
    bandRows=1 + min(options['radius'],centerRow) + min(header['nrows'] - (centerRow+1), options['radius'])
    yMin= max(header['yllcorner'],header['yllcorner'] + ((header['nrows'] - centerRow - options['radius'] - 1) * header['cellsize']))
    LLC = arcpy.Point(header['xllcorner'],yMin)
    # arcpy.Raster(inRaster).save('c:\\temp\\bandInRas')
    bandArray = arcpy.RasterToNumPyArray(inRaster,LLC,"#",bandRows,-9999)
    # newRaster = arcpy.NumPyArrayToRaster(bandArray,LLC,
                                         # header['cellsize'],header['cellsize'],-9999)
    # newRaster.save('c:\\temp\\bandarray')
    return bandArray

# @profile
def addData(cumCurrentArray, currentArray, centerRow, centerCol, options):
    #xprint 'adding data for center row, col, rad'
    #xprint centerRow,centerCol,options['radius']
    minRow = max(0,centerRow-options['radius'])
    maxRow = min(minRow+cumCurrentArray.shape[0]-1, centerRow+options['radius'])
    minCol = max(0,centerCol-options['radius'])
    # fixme: check if next line needs mincol added like in maxrow line
    maxCol = min(cumCurrentArray.shape[1]-1, centerCol+options['radius'])

    fullCurrentArray = npy.zeros(cumCurrentArray.shape,dtype='float64')
    fullCurrentArray[minRow:maxRow+1,minCol:maxCol+1]=currentArray
    cumCurrentArray += fullCurrentArray
    return cumCurrentArray

@profile
def addData_arcpy(cumCurrentRaster, currentRaster):
    descData=arcpy.Describe(cumCurrentRaster)
    arcpy.env.extent=descData.Extent    
    newCumCurrentRaster = arcpy.sa.Con(arcpy.sa.IsNull(currentRaster),cumCurrentRaster,arcpy.sa.Plus(currentRaster,cumCurrentRaster))
    return newCumCurrentRaster

@profile
def circ(array, rowArray, colArray, centerRow, centerCol, options):
        
        #xprint 'rowArray',rowArray.shape
        #xprint rowArray
        #xprint 'centerRow',centerRow
        #xprint 'colArray',colArray.shape
        #xprint colArray
        #xprint 'centerCol',centerCol
        # st0=time.clock()
        # distArray0 = npy.multiply((rowArray - centerRow), (rowArray- centerRow)) + npy.multiply((colArray-centerCol), (colArray-centerCol))
       
        # distArray = npy.sqrt(distArray0)
        # st0=elapsed_time(st0)
        #xprint 'x2'

        startRow = max(centerRow - options['radius'],0)
        endRow = min(centerRow + options['radius'],array.shape[0]-1)
        startCol = max(centerCol - options['radius'],0)
        endCol = min(centerCol + options['radius'],array.shape[1]-1)

        colArraySmall = (colArray[startRow:endRow+1,startCol:endCol+1]).astype('float64')
        rowArraySmall = (rowArray[startRow:endRow+1,startCol:endCol+1]).astype('float64')
        arraySmall=array[startRow:endRow+1,startCol:endCol+1]
        
        #xprint 'arraySmall',arraySmall
        #xprint'colsmall',colArraySmall
        #xprint'rowsmall',rowArraySmall

        # fastest to do this piecemeal
        # gc.collect()
        #xprint 'cr',centerRow
        #xprint 'cc',centerCol
        #xprint'r0',rowArraySmall
        #xprint'c0',colArraySmall

        # # Code below could be used if processing large radii creates memory problems
        # rowArraySmall-=centerRow
        # colArraySmall-=centerCol
        # rowArraySmall*=rowArraySmall
        # colArraySmall*=colArraySmall
        # colArraySmall+=rowArraySmall       
        # npy.sqrt(colArraySmall, colArraySmall) # take sqrt in place
        #xprint'sqrt'
        # now = time.clock()
        distArray = npy.sqrt(npy.multiply((rowArraySmall - centerRow), (rowArraySmall- centerRow)) + npy.multiply((colArraySmall-centerCol), (colArraySmall-centerCol))) #causes memory error
        # now=elapsed_time(now)
        del rowArraySmall,colArraySmall
        #xprint 'dist',distArray
        
        
        arrayMasked = npy.where(distArray <=options['radius'], arraySmall, -9999) #fixme do nans here

        

        gc.collect()
        # st0=elapsed_time(st0)
        # distArray = npy.sqrt(npy.multiply((rowArray - centerRow), (rowArray- centerRow)) + npy.multiply((colArray-centerCol), (colArray-centerCol))) #causes memory error
        #xprint arrayMasked
        return arrayMasked

# @profile        
def center_block(array, options, centerRow, centerCol):
        # returns array of same shape as array (i.e. entire radius), but everything outside of center blkock is zero
        startRow = centerRow - ((options['blockSize']-1)/2)
        endRow = centerRow + ((options['blockSize']-1)/2)
        startCol = centerCol - ((options['blockSize']-1)/2)
        endCol = centerCol + ((options['blockSize']-1)/2)
        blockArray = npy.zeros(array.shape, dtype = 'float64')# - 9999  # replace -9999 with nan?  
        blockArray[startRow:endRow+1,startCol:endCol+1] = array[startRow:endRow+1,startCol:endCol+1]
        return blockArray            

@profile
def center_mask(array, options, centerRow, centerCol): #fixme not used
        print 'masking',centerRow,centerCol
        startRow = centerRow - ((options['blockSize']-1)/2)
        endRow = centerRow + ((options['blockSize']-1)/2)
        startCol = centerCol - ((options['blockSize']-1)/2)
        endCol = centerCol + ((options['blockSize']-1)/2)
        maskedArray = array
        maskedArray[startRow:endRow+1,startCol:endCol+1] = 0
        return maskedArray            

@profile
def get_center_block(array, options, centerRow, centerCol): #fixme not used
        print 'getting center block',centerRow,centerCol
        startRow = centerRow - ((options['blockSize']-1)/2)
        endRow = centerRow + ((options['blockSize']-1)/2)
        startCol = centerCol - ((options['blockSize']-1)/2)
        endCol = centerCol + ((options['blockSize']-1)/2)
        blockArray = npy.where(array > 0, -9999, -9999)
        blockArray[startRow:endRow+1,startCol:endCol+1] = 1
        return blockArray            


def get_subset_header(array, fullHeader, options, centerRow, centerCol):
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

def get_cs_path():
    """Returns path to Circuitscape installation """
    envList = ["ProgramW6432", "ProgramFiles", "ProgramFiles(x86)"]
    for x in range (0,len(envList)):
        try:
            pfPath = os.environ[envList[x]]
            csPath = os.path.join(pfPath,'Circuitscape\\cs_run.exe')
            if os.path.exists(csPath): return csPath
        except: pass
    return 'D:\\Program Files\\Circuitscape\\cs_run.exe'
    
def elapsed_time(start_time):
        """Returns elapsed time given a start time"""
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

@profile        
def ascii_grid_reader(filename, nodata,data_type):
        """Reads rasters saved as ASCII grids or numpy arrays into Circuitscape."""
       
        if nodata == False:
            pmap = npy.loadtxt(filename, skiprows=5, dtype=data_type)
        else:
            pmap = npy.loadtxt(filename, skiprows=6, dtype=data_type)
            pmap = npy.where(pmap==nodata, -9999, pmap)

        return pmap


def get_header(filename):
    header = {}
    descData=arcpy.Describe(filename)
    cellsize=descData.meanCellHeight
    extent=descData.Extent
    xllcorner = extent.XMin
    yllcorner = extent.YMin
    yulcorner = extent.YMax
    xlrcorner = extent.XMax
    nrows=int((yulcorner-yllcorner)/cellsize)
    ncols=int((xlrcorner-xllcorner)/cellsize)
    header = {}
    header['ncols'] = ncols
    header['nrows'] = nrows
    header['xllcorner'] = xllcorner
    header['yllcorner'] = yllcorner
    header['cellsize'] = cellsize
    header['nodata'] = -9999 #if (nodata == False) else nodata 
    return header 

@profile    
def ascii_grid_writer(file_name, data, header, options):
    """Writes rasters to ASCII grid or numpy formats."""     

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

@profile    
def export_ras_to_npy(raster,npyFile):
    
    descData=arcpy.Describe(raster)
    cellSize=descData.meanCellHeight
    extent=descData.Extent
    spatialReference=descData.spatialReference
    
    pnt=arcpy.Point(extent.XMin,extent.YMin)
    outData = arcpy.RasterToNumPyArray(raster,"#","#","#",-9999)
    #outData = npy.where(outData==noDataVal,-9999,outData)
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

@profile
def import_npy_to_ras(npyFile,baseRaster,outRasterPath):
    npyArray = npy.load(npyFile, mmap_mode=None)
    npyArray=npyArray.astype('float32')
    descData=arcpy.Describe(baseRaster)
    cellSize=descData.meanCellHeight
    extent=descData.Extent
    spatialReference=descData.spatialReference
    
    pnt=arcpy.Point(extent.XMin,extent.YMin)
    newRaster = arcpy.NumPyArrayToRaster(npyArray,pnt,
                                         cellSize,cellSize,-9999)
    if outRasterPath is not None:
        newRaster.save(outRasterPath)
        return
    else:
        return newRaster
        

def write_header(raster,numpyArray,numpyFile):
    
    ncols=numpyArray.shape[1]
    nrows=numpyArray.shape[0]
    descData=arcpy.Describe(raster)
    cellSize=descData.meanCellHeight
    extent=descData.Extent
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

def print_failure(numResistanceNodes, memFlag, sleepTime):
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


def setCircuitscapeOptions(resisAsciiFile,sourceAsciiFile,groundAsciiFile,outputFile):
    """Sets default options for calling Circuitscape.

    """
    csOptions = {}

    csOptions['data_type']='raster'
    csOptions['version']='unknown'
    csOptions['low_memory_mode']=False
    csOptions['scenario']='advanced'
    csOptions['habitat_file']='(Browse for a habitat map file)'
    csOptions['habitat_map_is_resistances']=True
    csOptions['point_file']=('(Browse for file with '
                          'locations of focal points or areas)')
    csOptions['point_file_contains_polygons']=True
    csOptions['connect_four_neighbors_only']=False
    csOptions['connect_using_avg_resistances']=True
    csOptions['use_polygons']=False
    csOptions['polygon_file']='(Browse for a short-circuit region file)'
    csOptions['source_file']='(Browse for a current source file)'
    csOptions['ground_file']='(Browse for a ground point file)'
    csOptions['ground_file_is_resistances']=True
    csOptions['use_unit_currents']=False
    csOptions['use_direct_grounds']=False
    csOptions['remove_src_or_gnd']='rmvsrc'
    csOptions['output_file']='(Choose a base name for output files)'
    csOptions['write_cur_maps']=True
    csOptions['write_cum_cur_map_only']=True
    csOptions['log_transform_maps']=False
    csOptions['write_volt_maps']=False
    csOptions['solver']='cg+amg'
    csOptions['compress_grids']=False
    csOptions['print_timings']=False
    csOptions['use_mask']=False
    csOptions['mask_file']='None'
    csOptions['use_included_pairs']=False
    csOptions['included_pairs_file']='None'
    csOptions['use_variable_source_strengths']=False
    csOptions['variable_source_file']='None'
    csOptions['write_max_cur_maps']=False
    csOptions['set_focal_node_currents_to_zero']=True
    csOptions['print_timings']=False
    csOptions['output_file'] = outputFile
       
    csOptions['remove_src_or_gnd']='rmvsrc'#'keepall'
    csOptions['habitat_file'] = resisAsciiFile
    csOptions['source_file'] = sourceAsciiFile
    csOptions['ground_file']= groundAsciiFile

    return csOptions

def writeCircuitscapeConfigFile(configFile, csOptions):
    """Creates a configuration file for calling Circuitscape.

    """
    config = ConfigParser.ConfigParser()

    sections={}
    section='Version'
    sections['version']=section

    section='Connection scheme for raster habitat data'
    sections['connect_four_neighbors_only']=section
    sections['connect_using_avg_resistances']=section

    section='Short circuit regions (aka polygons)'
    sections['use_polygons']=section
    sections['polygon_file']=section

    section='Options for advanced mode'
    sections['source_file']=section
    sections['ground_file']=section
    sections['ground_file_is_resistances']=section
    sections['use_unit_currents']=section
    sections['use_direct_grounds']=section
    sections['remove_src_or_gnd']=section

    section='Calculation options'
    sections['solver']=section
    sections['print_timings']=section
    sections['low_memory_mode']=section

    section='Output options'
    sections['output_file']=section
    sections['write_cur_maps']=section
    sections['write_cum_cur_map_only']=section
    sections['log_transform_maps']=section
    sections['write_volt_maps']=section
    sections['compress_grids']=section
    sections['write_max_cur_maps']=section
    sections['set_focal_node_currents_to_zero']=section

    section='Mask file'
    sections['use_mask']=section
    sections['mask_file']=section

    section='Options for pairwise and one-to-all and all-to-one modes'
    sections['use_included_pairs']=section
    sections['included_pairs_file']=section
    sections['point_file']=section
    sections['point_file_contains_polygons']=section

    section='Options for one-to-all and all-to-one modes'
    sections['use_variable_source_strengths']=section
    sections['variable_source_file']=section

    section='Habitat raster or graph'
    sections['habitat_file']=section
    sections['habitat_map_is_resistances']=section

    section="Circuitscape mode"
    sections['scenario']=section
    sections['data_type']=section

    if csOptions['ground_file_is_resistances']=='not entered':
        csOptions['ground_file_is_resistances'] = False
    if csOptions['point_file_contains_polygons']=='not entered':
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


def call_circuitscape(CSPATH, outConfigFile):
    
    memFlag = False
    failFlag = False
    print('Calling Circuitscape:')
    proc = subprocess.Popen([CSPATH, outConfigFile],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, 
                           shell=True)
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
    output=proc.communicate()[0]
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

def gprint(string):
    gp.addmessage(string)
    try:
        if cfg.LOGMESSAGES:
            write_log(string)
    except:
        pass

def create_log_file(options):
    ft = tuple(time.localtime())
    timeNow = time.ctime()
    fileName = 'omni_LOG_'+options['outputFileText']+'.txt'
    filePath = os.path.join(options['outputDir'],fileName)
    try:
        logFile=open(filePath,'a')
    except:
        logFile=open(filePath,'w')
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
        logFile=open(options['logFilePath'],'a')
    except:
        logFile=open(options['logFilePath'],'w')
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

if __name__ == '__main__':
    omniscape(options)

print 'Done'    

# beamtiles
# for tile, take 3x3 subset
# beam up, down, etc

# OR wedges.... could test by running 4x, with code that snuffs out sources outsid of each wedge
    # Doesn't seem to help, outer circles show up more


# # VOLTAGE###
# dist matrix, get everything within options['radius']
# fill options['radius'] with max diff
# BUT WOULD NEED TO DO MOVING WINDOW ACROSS ENTIRE MATRIX LIKE FOC STATS
# http://stackoverflow.com/questions/10996769/pixel-neighbors-in-2d-array-image-using-python
# https://geonet.esri.com/thread/16548
# ndimage.filters:
# http://gis.stackexchange.com/questions/34306/focal-statistics-with-variable-options['radius']
# circle edges showing up... could use a kernel that fades to outside.
# Somewhat but not terribly worse than without fading.
#norm voltages doesn't seem to be good- more circle edges showing up
    #could do center ground currents  without voltages...
# fading (normalization) matches VERY well with standard results using different color ramps (quintile, equal int, geom)



#speed


# Focal range to ID potential analysis points
# create center band
# next valid target
# lookup
#DONT USE ascii inputs, slows bands etc way down
# Add del statements at end of bands, at continue statements

# create 1/0 sourceraster at start, make sourcecenterarray early on
# if a centerarray has a 1, find 1st instance using sum and argmax

    # else:
        # sourceRaster = arcpy.sa.Con((arcpy.Raster(resisRaster) > options['rCutoff']),0,1)

        
# blockMap-map of targets with at least one valid target
# iterate through blockmap
# bandnum, stripenum
# focalmax
# generalize

# hsplit bandarray

# just take band, get stripes with valid targets

# failing on cospatial:Traceback (most recent call last):
  # File "D:\GIS_DATA\NACR\McRae\circuitscapeBraidedThruway\circleFlowBand.py", line 1600, in <module>
    # circleFlowBand()
  # File "D:\GIS_DATA\NACR\McRae\circuitscapeBraidedThruway\circleFlowBand.py", line 580, in circleFlowBand
    # grid = npy.indices((bandArray.shape))
  # File "D:\Python27\ArcGIS10.2\lib\site-packages\numpy\core\numeric.py", line 1602, in indices
    # res = empty((N,)+dimensions, dtype=dtype)
# MemoryError

# Tiling-manual for now
# def lookup_stripe
# jumpstripe

# speed up code using below?
# import numpy as np
# N = 5
# aa = np.arange(-N,N).reshape(2,5)
#xprint aa
# aaa = np.where(aa>N/2,1,0)
#xprint aaa
# sum = np.sum(aaa, axis=0)
# ind = np.argmax(sum==1)
#xprint 'looking for '+str(N/2)
#xprint ind
# blarg



# # SOLVE TIMES
    # #r=300:15
    #r=150: 2-5
# source strength ppn to tdiff?
# Next- redo 10km with 1C change, then 50km
# climate code too slow when iterating through targets

#NEED TO FIX SQUARES. FADE STARTING AT OUTSIDE EDGE?
    # This works great. Could be problematic for short distance connections and climate (since that only goes one direction and doesn't add two full opposite cur maps)
    # For directional climate, run twice, once with targets at center and once with sources at center?
    # For non-dir climate, just run once activating anything that differs by cutoff?
    # fadeindist of options['radius']/2 also gets rid of blocking. But nice symmetry to options['radius']- add 2 wedges, get a line.
    
# CLIMATE MULTISCALE- do long distance movements with big donut hole at coarse cell sizes. Get latitudinal movements, etc.
    # Do short distance movemetns to get elevation etc.
    # run at several radii, weight radii the same, have current cancel each other out

# Speed up search for non-nodata stripes by finding 1st instance in band?

# Next: climate function- if 0 < deltaT < 5, i = deltaT. If 5 < deltaT < 10, i=10-deltaT
# for climate, scaling sources by # targets and v-v. Just multiply current map by targetsum, i.e. use target only?
# budget climate dispersers based on distance from source?

# set high R to areas with climate gradients (would affect contouring movements as well)?
# slope- could square it to make steep gradients more costly?
    # would this result in a pth following gentle gradient more than zero for a long ways then one pixel steep?
    # no- same cost unless gave proportionally higher cost to higher gradients
    #APPEARS TO DO WELL! TRY AGAINST CLIMATE LM ON MAC
        # NOT BAD. TRY WITH ORIGINAL CLIMATE LAYER INSTEAD OF RESAMPLED.
    # How to scale- would need to multiply and add 1?
    # 0--1
    # 0.1--1.1 or 101 
    #TRY
    
    # NOT HARD to add climate to CS. Use slope as placeholder, estimate 2 days to modify _construct_g_graph in compute_base
    
    
# # MULTISCALE- coarser cell and block size for lng distances????
# ADD COARSENING feature?

# flipSrcTarg- set targets as sources and v-v

# get sources, targs
# if options['useClimate']
    # sourceArrayNew=zeros
    # # srow, scol= where(sourceArray>0)
    # trows, tcols=where(targetArray>0)
    # for each target point
        # tDiffArray=temp-tTarget
        # sourceArrayNew=sourceArrayNew+sourceArray*f(tDiff)
            # i ~ tdiff, quality
            # for now: i=1 if tdiff > 1, otherwise 0
        # each source = mean (f(tdiffs))

# TIME TO DO CLIMATE?
# only way to do directional/faunal flow.
# would need to track vectors, not just mags!
#zero connectivity in plains tho!!!

 
# fixme: fadein is causing negative numbers in voltages and currents with large radii
    
# Idea- current injected ~ numtargs.
# But then mask out target block and put ground in center.
# vastly increases run time:
            # rows,cols=npy.where(targetArray > 0)
            # groundArray=npy.zeros(sourceArray.shape,dtype='int32')
            # groundArray[rows[0],cols[0]] = 100000

# fade in doing trick. test results against much finer run, like 15


#vdiff 
# consider sum too- highlight single pixels causing barriers in mult dir's
# consider looking at maximum difference edge-to-edge of window... same idea. pixel in middle of big diff highlighted
    # i.e., range of vup, vdn etc.


# try 1/r

# re-try w2w, but with permeable water?

# Long term
# current:
    # -incorporate cost distance or resist dist
    # -distance
    # -temp diff
# connected components (not needed if center ground?)


#  try having the entire landscape supply and receive current (except water). could bring out connective areas with no weighting  by what's connected.
  
# donut could help a bit by reducing build up and emphasizing long distance conn (lowering current in everglades)
    # ANS NOT REALLY

# 1)set block cur to ns*ntarg before adding
# 2)Can w2w correct? looking to connect each point to each edge. ideally edge pixels inject as much current as flowing in

#looking like large blocks (even 11 pix) are problematic when have intact patches.
# try w2w with SE data- how well does it approximate this
    # great in some areas, but terrible for eg everglades. small exit area in southern tip of fl. at minimum would need to overlap tiles

# could do 50km tiles, overlapping...

# Short circuit block? with high resis?.. ANS doesnt' seem to help, get squares
# attribute entire block with average block current?

# NULL model to get rid of edges???? ANS still not great... see circles more, weird distn.

# fixme- create file base at beginning, use throughout
# fixme- stamp all nodata inputs as nodata in final
# 

# TODO
# modularize to run bands one at a time on server?
# fixme- need connected components. for comp in components, if has srces and targs, scale appropriately


# Question: why have halo? is it for stitching? if so how would overlap be handled?
# answer- makes things much easier. don't need to adjust center row, col when smashed 
# time trials- stitching together small grids with arcpy vs array
# numpy is MUCH faster.


# import math 

# # options['resisRasterBase'] = 'resis100km.asc'
# # options['resisRasterBase'] = 'ecolresis_100m.asc'

# # options['resisRasterBase'] = 'ecol_resis3_nodata2_10kbarriers2.asc'
# # options['resisRasterBase'] = 'EcologyResistances100k.asc'
# # options['sourceRasterBase'] = 'EcologySourcePoints.asc'
# # options['resisRasterBase'] = 'ecolres_halo.asc'
# # options['resisRasterBase'] = 'nlcdres_ne_NODATA.ASC'
# # options['resisRasterBase'] = 'resistance_surface_se.asc'
# # options['resisRasterBase'] = 'resistance_surface_se_nodata_270ma.asc'
# # options['resisRasterBase'] = 'tip_peninsula_resis_270m.ASC'
# # options['resisRasterBase'] = 'klamath2_r_md2.asc'
# # options['resisRasterBase'] = 'EcologyResistances.asc'
# options['resisRasterBase'] = '5x5R.asc'
# # options['resisRasterBase'] = 'ecol_resis3_nodata2_10kbarriers2.asc'
# options['resisRasterBase'] = 'ecol_resis3_nodata2.asc'

# options['useSourceRaster']=False
# options['sourceRasterBase'] = '5x5_src.asc'
# options['sourceRasterBase'] = 'proof_points39.asc'

# options['useClimate']=False
# options[tCutoff']=3
# options['climateRasterBase'] ='5x5Temperature.asc'
# options['climateRasterBase'] = 'ecologytemperature.asc'

# options['projectDir'] = 'C:\\Dropbox\\Working\\Circuitscape_BraidedThruway\\CF_PROOF'
# options['outputDirBase'] = 'CF_5x5'
# options['outputDirBase'] = 'CF_climEcol'

# options['projectDirBase'] = 'scratch'

# options['projectDir'] = 'C:\\dropbox\\Working\\Circuitscape_BraidedThruway\\CF_NE_NLCD2'
# options['outputDirBase'] = 'NY'
# # options['resisRasterBase'] = 'ny_nlcd_500m.asc'
# options['resisRasterBase'] = 'ny_nlcd_166m_water5_nodata.asc'
# # options['resisRasterBase'] = 'ny_nlcd_166m_water5_nodata_adks2.asc'



# options['projectDir'] = 'C:\\dropbox\\Working\\Circuitscape_BraidedThruway\\CF_NE_CLIM_Coarse'
# options['outputDirBase'] = 'adks450m'
# # options['resisRasterBase'] = 'ny_nlcd_500m.asc'
# options['resisRasterBase'] = 'resis_sq_adks_450m.ASC'
# options['resisRasterBase'] = 'resis_sq_adks_450m.tif'

# options['projectDir'] = 'C:\\temp'

# options['outputDirBase'] = 'neTest'
# # options['resisRasterBase'] = 'ny_nlcd_500m.asc'
# options['resisRasterBase'] = 'resistance.tif'
# options['useClimate']=True
# options[tCutoff']=0.1
# options['climateRasterBase'] ='tmean.tif'


# options['projectDir'] = 'D:\\GIS_DATA\\NACR\\McRae\\circuitscapeBraidedThruway\\CF_NE'
# options['outputDirBase'] = 'tile_test7'
# # options['resisRasterBase'] = 'ny_nlcd_500m.asc'
# options['resisRasterBase'] = 'resis_chesapeak2.tif'
# options['useClimate']=True
# options[tCutoff']=0.5
# options['climateRasterBase'] ='tmean_ne_chesapeak.tif'
# tile=1
# xxxx


# options['projectDir'] = 'C:\\dropbox\\Working\\Circuitscape_BraidedThruway\\CF_PROOF'
# # options['resisRasterBase'] = 'EcologyResistances.asc'
# options['resisRasterBase'] = '6x9R.asc'
# options['rCutoff']=1
# options['useClimate']=False
# options['resisRasterBase'] = 'ecol_resis3_nodata2_10kbarriers2.asc'
# options['resisRasterBase'] = 'ecol_resis3_nodata2.asc'


# options['projectDir'] = 'C:\\Dropbox\\Working\\Circuitscape_BraidedThruway\\CF_stgr'
# options['outputDir'] = 'C:\\Dropbox\\Working\\Circuitscape_BraidedThruway\\CF_stgr\\cf_stgr_srcs'
# options['scratchDir'] = 'C:\\Dropbox\\Working\\Circuitscape_BraidedThruway\\CF_stgr\\scratch'

# options['resisRasterBase'] = 'typh_r.ASC'
# options['sourceRasterBase'] = 'typh_sources.ASC'
# options['projectDir'] = 'c:\\temp\\CF_BAND'
# options['outputDir'] = 'c:\\temp\\CF_BAND\\csbandecolrVr100k_VDacross'
# options['scratchDir'] = 'c:\\temp\\CF_BAND\\scratchbandecolr'


# PROBLEM: if a source is disconnected, then source and targ i's don't match!
# Consider just tying targets to ground. otherwise need to calculate pairwise r's 
# OR do connected components

# --or subset for each solve? would mean replacing circ code with arcpy. use bounding circle code?
# NEXT STEPS
# massive tiles
# This results in memory error for 10k x 10k raster:
# pmap = npy.loadtxt(filename, skiprows=6, usecols=range(2,3), dtype=data_type)
# Instead, could try load using memory mapping?
#       No, doesn't work. Can't create npy file in first place.
# Add iterative
# Add indep source/targ layer
# f(distance)
# f(temp)
# f(qual)


# fix negative current failures. Need for when there is barrier in block!? Also getting HOLES in target patches   
    #fixed,but would be best if there was ONE ground, wouldn't likely screw up flow pattern
    # NO- massively increases proc time
# TRY- ONLY do long distance movements- e.g. 20-50km?
# i = f(dist), f(quality)
#try larger blocks, see how results hold up (WELL)
#scan for blocks with targets first
   # where 
    # max in block

    # circle resis array (SUBSET of options['radius']*2 on a side, ND outside circle)
    # focal block- small or size of circle?
        # >>size of circle for now, keep simple, inspect in arcgis

# 1) pinchpoints
# 2) how connected is a pixel? record how many sources there are in window
# 3) where is a pixel relative to the network?
# 4) where is a pixel relative to clusters?
# 5) does a pixel help connect other pixels to each other? to edge of map?
# 6) how much flow through an area

#IDEAS
# have 10x10 focal window with targets. More targets in circle = more source amps

# block 10
# skip last rows and cols for now, can probably scale things to use them with smaller block later.

# compare with wall2wall. Having sources throughout grid and targets at edge of square (x4 dirs) or edge of circle

# clip res and other grids to be 2r high and wide
# add ability to read numpy source and ground files to CS
# subtract flow into target pixel (or track SEPARATELY)
# subtract flow out of source pixs (same as above?)
# subtract null
# i injected = f(distance). Could tail off, could only be on outer ring
# i injected = f(quality)
# targets = f(resil)
# i injected = f(temp difference)


            #temp targetArray = -targetArray*sourceSum
            #temp- causing failures? sourceArray = sourceArray*tar
            #temp sourceArray = sourceArray+targetArray

# If tiling slows down w2w, and focal window speeds this up, may be comparable.
# suorce amperage = # targets
# count up sources, each target has -N amperage
# step 11
# set sources
# mask out targets from sources

"""omniscape.py- just a set of functions using omniOptions
    iterates through blocks
    returns name of current (&voltage) map(s) 
omniscape_block.py- just what's required to run a single block? omniscape iterates through blocks and calls this.
run_omniscape.py- inputs. populates omniOptions
omniscape_range_shift and such- custom code to call using temp data or whatever

omniscape_climate
PROBLEM: right now it is fast because we only solve a block once. All valid sources with matching targets in block are invoked simultaneously
tdiff is many-to-one since target just needs to be >tdiff. CHANGE THIS? SHOW 1deg or 2deg or 4deg flow? Then could use bins.
?SOLUTION: either discrete bins OR climate module(s) called to modify source and target arrays JUST BEFORE solve.

BIG PROBLEM??
    Pairing- amount of current depends on number of targets or source-target pairs. But what about differences in climate pairings?
    ---one source, many matching destinations, one destination many matching sources, 

PRE-SOLVE CLIMATE MODULE
modify sources, targets, based on criteria in options:
mode='gradient': present day temp
mode='velocity': present and future climate var's. Assume 2 for now.
for velocity:
   block climate variables (up to 4 rasters)
   calculate euc climate distances between each center target cell and all other cells in block. use these to modify target & source strengths


PRE-RUN CLIMATE MODULE
take climate inputs
divide into bins, can be as fine as you like. 
    save present_bins, future_bins, then could solve blocks only once. or
for climateBin
    map present bin location
    map future matches
    --or-- map future bin location, present matches?
    export all sources (present) and targets (future). call omniscape.py

    
from omniscape import *

# Omniscape base code
    # inputs:
        # sources (strengths) 
        # targets (strengths)
        # resistance (can use cutoff to create binary src or target layers)
        # block size
        # targonly- amount of current is based on sum of target strengths in block. source strengths proportional to source strength layer (and distance function if appl)
        # options['radius']
        # calc_current
        # calc_flow_accum
        # distance function for source strengths
            # from math import *
            # x = 2
            # string = "sin(x)*x**2"
            # test = eval(string)
            #xprint test
        
        # fade??
# Climate calling code
    # inputs:
        # temp, or pc1, pc2 at T1, T2
        # for each temp or pc1/pc2 combo:
            # identify sources and targets
            # call omniscape_base 
# Partial out fade? Would get 3 outputs- raw, fade, raw-fade. Or just fade and raw-fade.  
# Add memory error catch, try again.
saveFade- sums up fades and saves. Will take more memory.



# sourceRaster
# useSourceRasterCutoff #binary 1/0 if above(?) this value or have negative values to denote below this value
    # srcRasCutoffVal=

# INDEP SOURCES AND TARGS    
# If not options['useSourceRaster'] or useTargetRaster:
    # sourceraster, targraster = everything < sRCutoff, tRCutoff in resisRaster
    # if options['useSourceRaster'] and not useTargetRaster:
        # sourceRaster,  = sourceRaster (strength = raster values)
        # targraster= everything < tRCutoff in resisRaster
    # if useTargetRaster and not sourceRaster:
        # targs = targetRaster (strength = raster values?)
        # sources = everything < sRCutoff in resisRaster
    # total flow could be 
        # 1) Sum(targetstrengths)*sum(sourcestrengths)
        # 2) Sum(targetstrengths) for targetOnly
            
    # targOnly: flow = ntargs if there is a source, else 0
    
    # sourceOnly? flow = nsources if there is a target, else 0
        # Doesn't really make sense? ntargs affected by block size, nsources not really

# Total flow ntargs*nsources, ntargs, sum(sourcestrenghts), sum(targetstrengths) 
# Indiv sources = source strength, tdiff, geogdist, effR       
"""        
        
    # ALSO NEED CLIMATE TARGONLY
   # larger areas = more dispersers absorbed
   # approximate results with bandsize=1
   # total current = ntargs that differ with at least one source
   # for each target, tally if there's any sources
    # after this line: tCutoffArray=npy.where(tDiffArray>=options[tCutoff'],1,0)

#4/16/15:
# removed polygonblocks
# remoed maskblocks
# removed fillblocks
# remove noweight?- doesn't make sense, since would produce diff results with diff block sizes.

# 8/10/15:
# removed fadeout
# 8/11
# removed straddleblocks
# removed donut
#removed fadevolt- may reinstate
#removed bufferdist- may reinstate
#removed logTrans- would give diff values for different block sizes I believe

# #todo:

# source and target counts? get 

# Fix striping in volt 
# scipy sparse csr_matrix can be 10x+ faster with where, but expensive to create and don't help with math operations. experiment with nans and avoiding where()
# use nans for nodata values
           
# Dave takes int(ln(flow)) and does a thin. also runs 2nd lcp through FA layer to clean up.
# separate source and target layers?

# Build in simple rewrite in case of write error. Could just do new scratch directories and a 'try_x_' in output files.
#automatically do int of Fa results and convert to polylines
#re-scale FA results 1-100?
#place ground node at center of mass of target block?

#Note: because of mis-alignment of block centers and source pixels, can get different flow lines when a pixel is a source and whena pixel is a target. Can happen even with blocksize=3, will increase with larger blocks.
# Could help to have ground at center of mass of target block

# Things to vary for initial experiments:
# Resistance cutoff for source cells (rCutoff)
# block size (blockSize)
# search radius for cells to connect (radius)
# PC distance cutoff (PCWindow)
# current ~ n targets (weightTargOnly = True) or current ~ n targets * n sources (weightTargOnly = False) 
# resistance raster scoring (I'm using R = (1+HM)^10)

# NOTE: this will fail if your project directory is in dropbox and is syncing.

#experiment that didn't go anywhere- fade with a nudge is better.
# @profile
# def apportion_currents(circleResisArray,currentArray,subsetCenterRow,subsetCenterCol,options,multiplier):
    # #fixme: test
    # centerResistanceArray = center_block(circleResisArray, options, subsetCenterRow, subsetCenterCol) 
    # centerConductanceArray = npy.where(centerResistanceArray>0,1/centerResistanceArray, 0)
    # centerCurrentArray = center_block(currentArray, options, subsetCenterRow, subsetCenterCol) 
    
    # sumCenterConductances = centerConductanceArray.sum() #need NAN here

    # print 'sumcond',sumCenterConductances
    # shareConductanceArray = centerConductanceArray/sumCenterConductances

    # # sumCenterCurrents = current at ground, or sum of currents 1 pixel away
    # # sumCenterCurrents= multiplier
    
    # #sumCenterCurrents= sum of ALL currents in center
    # sumCenterCurrents=centerCurrentArray.sum()  
    
    # print 'sumcur', sumCenterCurrents
    # centerCurrentArrayNew = shareConductanceArray * sumCenterCurrents
    # currentArray = npy.where(centerCurrentArrayNew>0,centerCurrentArrayNew, currentArray)
    # print 'centerResis'
    # print centerResistanceArray
    # print 'centerConduc'
    # print centerConductanceArray
    # print 'centerCurrentArrayNew'
    # print centerCurrentArrayNew
    # return currentArray
# Instead of fading, partition current in block according to resistance?
# doesn't work, neither does taking average current flowing in to edge cells of block and setting block to that value

#large block sizes can lead to edge effects I think. edges of study area lack local connections. Shows up with distance function results (200km rad)
# maybe not- #targonly reduces edge effects. above run was not targonly.
