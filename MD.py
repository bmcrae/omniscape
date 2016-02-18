3200000	km2 (about 2/5 area of conus)
3.2E+12	m2
3555555556	30m pixels
3555.555556	
	
3.5 billion pixels	
*250 thousand calcs	
	
8.88889E+14	calcs
888.8888889	
888 trillion calcs	

load bands for all N rasters
search band1
when data pixel found, do circ for all rasters

def std_euc(array1, array2, array3, array4, subsetCenterRow, subsetCenterCol):
    varC1 = options['varC1'] etc.
    C1 diffs, c2 diffs, etc...
    450m ~707 pixels= 250k comparisons
            pretty quick- see manual SEucDist code at end.
    
    
    






tileNum = 10 # 0 to ignore. ONLY WORKS ON SERVERS FOR ME. For running tiled bands to be stitched together later.
numTiles = 8 # number of horizontal tiles to break into if tileNum > 0.
#---------------------------------------------------------------------
# INPUTS #
#---------------------------------------------------------------------
options = {}

# MOVING WINDOW AND TARGET BLOCK SIZES ####
options['radius'] = 278# In PIXELS. Search radius, with sources activated within the radius and outside of the center (target) block.

options['projectDir'] = r'C:\DATADRIVE\DUKE_PNW_DATA\PNW_OmniScape' #r'D:\GIS_Data\bmcrae\Duke_PNW_Omniscape\d8_omniscape'# this is where all the input data are, and where out directory will be created.
options['c1RasterBase'] = 'R_d8_clpF_180m.tif'
options['c2RasterBase'] = 'NORM_6190_PC1_PNWmockup.tif'
options['c3RasterBase'] = 'NORM_6190_PC2_PNWmockup.tif'
options['c4RasterBase'] = 'MIROC5_2080s_RCP85_PC1_PNWmockup.tif'
options['c5RasterBase'] = 'MIROC5_2080s_RCP85_PC2_PNWmockup.tif'

options['c1Var'] = 10
options['c2Var'] = 10
options['c3Var'] = 10
options['c4Var'] = 10
options['c5Var'] = 10

options['outputDirBase'] = 'd8_50km'
options['printTimings'] = True

options['cleanUpBandFiles'] = True # set to false for debugging, will keep scratch files and individual band current maps

# OPTION TO LIMIT ANALYSIS EXTENT ####  
# Bands are horizontal, with width equal to blockSize. There are approx nrows/blockSize bands in a raster. Stripes are vertical, with width equal to blocksize.
# These options allow you to only process a subset of bands and stripes. 
options['startBand'] = 0# First band to process. Use 0 to ignore.
options['endBand'] = 0 # Stops after processing this band. Use 0 to ignore.
options['startStripe'] = 0 # First stripe to process. Use 0 to ignore.
options['endStripe'] = 0 # Stops after processing this stripe. 0 to ignore.


    
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

def topoclimate(options):
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
        

        # Set raster paths and export to ascii if needed. FIXME: don't really need ascii, just convenient for header code for now
        c1Raster = path.join(options['projectDir'],options['c1RasterBase'])
        descData = arcpy.Describe(c1Raster)
        spatialReference = descData.spatialReference
        options['cellSize'] = descData.MeanCellHeight
        
        c2Raster = c3Raster = c4Raster = c5Raster = None
        if options['c2RasterBase'] is not None:
            c2Raster = path.join(options['projectDir'],options['c2RasterBase'])
        if options['c3RasterBase'] is not None:
            c3Raster = path.join(options['projectDir'],options['c3RasterBase'])
        if options['c4RasterBase'] is not None:
            c4Raster = path.join(options['projectDir'],options['c4RasterBase'])
        if options['c5RasterBase'] is not None:
            c5Raster = path.join(options['projectDir'],options['c5RasterBase']) 
             
        header = get_header(c1Raster)          
        descData = arcpy.Describe(c1Raster)
        arcpy.env.extent = descData.Extent    
        
        
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

            c1BandArray = band(c1Raster,header,centerRow, options)

            if npy.max(c1BandArray) == -9999: continue
            c2BandArray = band(c2Raster, header, centerRow, options)
            c3BandArray = band(c3Raster, header, centerRow, options)
            c4BandArray = band(c4Raster, header, centerRow, options)
            c5BandArray = band(c5Raster, header, centerRow, options)
                    
    #fixme: getting extra nodata row on top of climate band?
                #xprint 'climband'
                #xprint climateBandArray
            bandDistArray = npy.zeros(c1BandArray.shape, dtype = 'float64') 
            
            subsetCenterRow = min(centerRow,options['radius'])
            
            # Check for all nodata in center band of C1 raster
            if options['blockSize']>1:
                bandCenterArray = c1BandArray[(subsetCenterRow-(options['blockSize']-1)/2):(subsetCenterRow+(options['blockSize']-1)/2),:]
            else:
                bandCenterArray = c1BandArray[subsetCenterRow,:]
             
            if npy.max(bandCenterArray) == -9999:       
                del bandCenterArray
                continue
            del bandCenterArray
            
            
            grid = npy.indices((c1BandArray.shape))
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
                c1CircleArray = circ(c1BandArray, rowArray, colArray, subsetCenterRow, centerCol, options)

                        
                circleHeader = get_subset_header(sourceArray, header, options, centerRow, centerCol)
                yMin = circleHeader['yllcorner'] #fixme- check. coudl be something like: max(circleHeader['yllcorner'],circleHeader['yllcorner'] + ((circleHeader['nrows'] - centerRow - options['radius'] - 1) * circleHeader['cellsize']))
                
                
                
            
                print '\nDone with prep'
                start_time0 = elapsed_time(start_time0)           

                # print 'Solving stripe #' + str(stripeNum) + '/' + approxEndStripe,' and band #' + str(bandNum) + '/' + approxEndBand+'\n'

                # maxSolvesDone=(bandNum-options['startBand']+1)*numstripes + (stripeNum-options['startStripe']+1)
                pctDone = report_pct_done(iter, maxNumSolvesToDo, -1)

                print 'Elapsed time so far: {0}'.format(datetime.datetime.now()-theStart)  
                
                solveInBand = True
                        
                if options['calcCurrent']:

                    curMapPath, outConfigFile = calc_current(options, csOptions, centerRow, centerCol) 

--snip--
                    if options['centerGround']:
                        maxCur = currentArray[subsetCenterRow, subsetCenterCol]
                    else:
                        maxCur = npy.max(currentArray)
                    if maxCur <= 0:
                        print 'NO CURRENT, continuing'
                        continue
                    
        
        
        
                del circlec1Array               
                start_time1 = time.clock()

    # FIXME: do all raster stuff for each band, not each solve
                # cumCurrentRaster = addData_arcpy(cumCurrentRaster, currentRaster)
    # FIXME: move file writing to band iteration

        
                
            
            
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
        
                  
                options = write_temp_maps(options,bandNum,cumCurrentRaster,cumVdiffRaster,cumFlowRaster)
                                         
            # print 'Done with band #',bandNum,',
            print 'Elapsed time so far: {0}'.format(datetime.datetime.now()-theStart)  
        print 'Done with solves.'  
        print_prof_data()
        write_final_maps(options,cumCurrentRaster,cumVdiffRaster,cumFlowRaster,cumSourceRaster,cumTargetRaster) 
        #xprint locals()
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
def get_raster_max(raster):
    try:
        maxObject = arcpy.GetRasterProperties_management(raster, "MAXIMUM") 
        max = float(str(maxObject.getOutput(0)))
        return max
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

            fileBase, fileExtension = path.splitext(options['c1RasterBase'])
            c1RasterSQ = path.join(options['scratchDir'],'sq_'+fileBase+'.tif')
            outResistanceRaster = arcpy.sa.Times(inRaster, inRaster) 
            outResistanceRaster.save(c1RasterSQ)           
            inRaster = c1RasterSQ    
        return inRaster
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

        
def sqr_root_inputs(c1Raster, sourceRaster,options):
    try:
        descData = arcpy.Describe(c1Raster)
        arcpy.env.extent = descData.Extent    

        fileBase, fileExtension = path.splitext(options['c1RasterBase'])
        c1RasterSQRT = path.join(options['scratchDir'],'sqrt_'+fileBase+'.tif')
        outResistanceRaster = arcpy.sa.SquareRoot(c1Raster) 
        outResistanceRaster.save(c1RasterSQRT)           
        c1Raster = c1RasterSQRT    
            
        
        fileBase, fileExtension = path.splitext(options['sourceRasterBase'])
        sourceRasterSQRT = path.join(options['scratchDir'],'sqrt_'+fileBase+'.tif')
        outSourceRaster = arcpy.sa.SquareRoot(sourceRaster) 
        outSourceRaster.save(sourceRasterSQRT)           
        sourceRaster = sourceRasterSQRT            
        return c1Raster, sourceRaster
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()
        
            
@profile        
def match_climate_pcs(sourceArray, targetArray, t1PC1BandArray, t1PC2BandArray, t2PC1BandArray, t2PC2BandArray, rowArray, colArray, subsetCenterRow, centerCol, options):                
    try:
        print 'Starting climate using PCs'
        start_timeClimate = time.clock()
        t1PC1Array = circ(t1PC1BandArray, rowArray, colArray, subsetCenterRow, centerCol, options)              
        t1PC2Array = circ(t1PC2BandArray, rowArray, colArray, subsetCenterRow, centerCol, options)              
        t2PC1Array = circ(t2PC1BandArray, rowArray, colArray, subsetCenterRow, centerCol, options)              
        t2PC2Array = circ(t2PC2BandArray, rowArray, colArray, subsetCenterRow, centerCol, options)              
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
                    cumVdiffFile = os.path.join(options['outputDir'],'BAND'+str(bandNum)+options['c1RasterBase']+'_vdiffSumtemp2.tif')
                    cumVdiffRaster.save(cumVdiffFile)
                except:
                    print 'Second error writing'
                    time.sleep(5)
                    cumVdiffFile = os.path.join(options['outputDir'],'BAND'+str(bandNum)+options['c1RasterBase']+'_vdiffSumtemp3.tif')
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
        options['blockSize'] = 1 
        if options['tileNum'] > 0:
            if options['tileNum'] > options['numTiles']:
                print 'Error: tile number (tileNum) cannot be higher than number of tiles (numTiles).'
                raw_input('\nPress Enter to continue.') 
                exit(0)    

            c1Raster = path.join(options['projectDir'],options['c1RasterBase'])
            header = get_header(c1Raster)
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
            # fileBase,ext = os.path.splitext(options['c1RasterBase'])
            # options['c1RasterBase'] = fileBase+str(tile)+ext
            # if options['useClimate']:
                # if not options['matchClimatePCs']:
                    # fileBase,ext = os.path.splitext(options['climateRasterBase'])
                    # options['climateRasterBase'] = fileBase+str(tile)+ext
                # else:
                    # fileBase,ext = os.path.splitext(options['c2RasterBase'])
                    # options['c2RasterBase'] = fileBase+str(tile)+ext
            
                    # fileBase,ext = os.path.splitext(options['c4RasterBase'])
                    # options['c4RasterBase'] = fileBase+str(tile)+ext

                    # fileBase,ext = os.path.splitext(options['c3RasterBase'])
                    # options['c3RasterBase'] = fileBase+str(tile)+ext

                    # fileBase,ext = os.path.splitext(options['c5RasterBase'])
                    # options['c5RasterBase'] = fileBase+str(tile)+ext
                    
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
            
        resisRasText, fileExtension = os.path.splitext(options['c1RasterBase'])
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

        print 'c1RasterBase = ',options['c1RasterBase']
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

def get_max_vdiff(voltMap,circlec1Array,csOptions):
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
            
            resis_l, resis_r = get_horiz_neighbors(circlec1Array)
            resis_u, resis_d = get_vert_neighbors(circlec1Array)
            
            
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
        if inRaster is None: return None
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

@profile
def circ(array, rowArray, colArray, centerRow, centerCol, options):
    try:
        startRow = max(centerRow - options['radius'],0)
        endRow = min(centerRow + options['radius'],array.shape[0]-1)
        startCol = max(centerCol - options['radius'],0)
        endCol = min(centerCol + options['radius'],array.shape[1]-1)
        colArraySmall = (colArray[startRow:endRow+1,startCol:endCol+1]).astype('float64')
        rowArraySmall = (rowArray[startRow:endRow+1,startCol:endCol+1]).astype('float64')
        arraySmall = array[startRow:endRow+1,startCol:endCol+1]
        
        distArray = npy.sqrt(npy.multiply((rowArraySmall - centerRow), (rowArraySmall- centerRow)) + npy.multiply((colArraySmall-centerCol), (colArraySmall-centerCol))) #causes memory error

        del rowArraySmall,colArraySmall
        # arrayMasked = npy.where(distArray <= options['radius'], arraySmall, -9999) #fixme do nans here
        rows,cols = npy.where((distArray <= options['radius'])&(arraySmall>-9999))
        vector=(arraySmall[rows,cols]).flatten()
        gc.collect()

        return vector
        
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

# @profile        
def center_block(array, options, centerRow, centerCol):
    try:
        # returns array of same shape as array (i.e. entire radius), but everything outside of center blkock is zero
        startRow = centerRow - ((options['blockSize']-1)/2)
        endRow = centerRow + ((options['blockSize']-1)/2)
        startCol = centerCol - ((options['blockSize']-1)/2)
        endCol = centerCol + ((options['blockSize']-1)/2)
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
    options = topoclimate(options)

    # Need this outside of main loop to avoid file lock problems:
    clean_up(options)
    if tileNum>0:
        raw_input('Done. Press Enter to continue.')
    print 'Done'    

    
    
# # ---run with arcgis win32 python... that's the only one I have scipy setup for

# import numpy as np
# import scipy.spatial

# print scipy.__version__


# # IS THIS IT???
# # http://stackoverflow.com/questions/21003272/difference-between-all-1d-points-in-array-with-python-diff
# # a = np.random.random(700)
# # print 'a0 a1 a2',a[0],a[1],a[2]
# # print a[0]-a[1]
# # print a[0]-a[2]
# # diff=np.subtract.outer(a,a)[np.tril_indices(a.shape[0],k=-1)]
# # print 'diff',diff
# # print diff.size

# Z1 = np.random.random(700)
# Z2 = np.random.random(700)
# Var1=np.var(Z1)
# Var2=np.var(Z2)

# diff1 = np.subtract.outer(Z1,Z1)[np.tril_indices(Z1.shape[0],k=-1)]
# diff2 = np.subtract.outer(Z2,Z2)[np.tril_indices(Z2.shape[0],k=-1)]
# diff1Sq = np.multiply(diff1,diff1)
# diff2Sq = np.multiply(diff2,diff2)
# diff1SqDivVar = np.divide(diff1Sq,Var1)
# diff2SqDivVar = np.divide(diff2Sq,Var2)

# SEucDist = np.sqrt(diff1SqDivVar + diff2SqDivVar)
 
# print 'Manual calc'
# print SEucDist
# print 'Manual mean',np.mean(SEucDist)


# Var = [Var1,Var2]
# Z=np.zeros((700,2),dtype='float64')
# Z[:,0]=Z1
# Z[:,1]=Z2

# Y1 = scipy.spatial.distance.pdist(Z, 'seuclidean', V=Var)#=None)
# Y2 = scipy.spatial.distance.pdist(Z, 'seuclidean', V=None)#=None)
# print 'Y1',Y1
# print 'Y2',Y2
# print 'average Y1, Y2',np.mean(Y1),np.mean(Y2)
# print Y1.size

# # a = np.array(range(5, 10))
# # b = np.array(range(1, 6))

# # res = a[np.newaxis, :] - b[:, np.newaxis]
# # print 'a',a
# # print 'b',b
# # print(res)

# blarg

# # use numpy?
# # http://stackoverflow.com/questions/26076576/subtract-all-pairs-of-values-from-two-arrays


# # for Var1, get all (x1-y1)
# # Diffxy[x,y] = Var1[x]-Var1[y]


# # http://stackoverflow.com/questions/17936587/in-numpy-find-euclidean-distance-between-each-pair-from-two-arrays

# v1=(np.array([[1,2,3,4],[0,0,0,0],[0,0,0,0],[0,0,0,0]])).transpose()
# v2=(np.array([1,3,4,6])).transpose()
# v3=np.array([[1,3,4,6],[0,0,0,0],[0,0,0,0],[0,0,0,0]])

# test=v1-v3
# print v1
# print v3

# print 'test'
# print test
# print np.multiply(test,test)
# blar
# Var = np.array([1.0,3.0])
# # Var=Var.transpose()

# print v1-v2
# print np.multiply((v1-v2),(v1-v2))
# distances = (v1-v2)^2
# distances = distances.sum(axis=-1)
# print 'dist',distances



# X=np.zeros((2,4),dtype='int32')
# X[0,:]=v1
# X[1,:]=v2
# X = X.transpose()
# print X
# Y = scipy.spatial.distance.pdist(X, 'seuclidean', V=Var)#=None)

# print Y    