#fixme: need to create a mask raster (con(c1 not null, c2 not null etc) and then create temporary c1, c2 rasters that have been masked

#fixme: add back in prev files
#fixme: lots of > 0 checks, but need NEGATIVE values.


tileNum = 0 # 0 to ignore. ONLY WORKS ON SERVERS FOR ME. For running tiled bands to be stitched together later.
numTiles = 8 # number of horizontal tiles to break into if tileNum > 0.
#---------------------------------------------------------------------
# INPUTS #
#---------------------------------------------------------------------
options = {}

# MOVING WINDOW AND TARGET BLOCK SIZES ####
options['radius'] = 450/90# In PIXELS. Search radius, with sources activated within the radius and outside of the center (target) block.
options['blockSize'] = 1

options['projectDir'] = r'C:\MD' #r'D:\GIS_Data\bmcrae\Duke_PNW_Omniscape\d8_omniscape'# this is where all the input data are, and where out directory will be created.
options['c1RasterBase'] = 'CTI_MR_90ma.tif'#
options['c2RasterBase'] = 'HLI_MR_90ma1.tif'#'c2.asc'#NORM_6190_PC1_PNWmockup.tif'
options['c3RasterBase'] = None#'c3.asc'#'NORM_6190_PC2_PNWmockup.tif'
options['c4RasterBase'] = None#'c4.asc'#'MIROC5_2080s_RCP85_PC1_PNWmockup.tif'
options['c5RasterBase'] = None#'resistances.asc'#'MIROC5_2080s_RCP85_PC2_PNWmockup.tif'


options['outputDirBase'] = 'HLICTI_MR90m'
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
from scipy.spatial.distance import euclidean, pdist, wminkowski, squareform
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
        c2Raster = c3Raster = c4Raster = c5Raster = None
        options['c2Var'] = options['c3Var'] =options['c4Var'] =options['c5Var'] = None
        c1Raster = path.join(options['projectDir'],options['c1RasterBase'])
        arcpy.CalculateStatistics_management(c1Raster, "1", "1", "#")
        options['c1Var'] = get_raster_var(c1Raster)

        c2BandArray = c3BandArray = c4BandArray = c5BandArray = None
        if options['c2RasterBase'] is not None:
            c2Raster = path.join(options['projectDir'],options['c2RasterBase'])
            arcpy.CalculateStatistics_management(c2Raster, "1", "1", "#")
            options['c2Var'] = get_raster_var(c2Raster)
        if options['c3RasterBase'] is not None:
            c3Raster = path.join(options['projectDir'],options['c3RasterBase'])
            arcpy.CalculateStatistics_management(c3Raster, "1", "1", "#")
            options['c3Var'] = get_raster_var(c3Raster)
        if options['c4RasterBase'] is not None:
            c4Raster = path.join(options['projectDir'],options['c4RasterBase'])
            arcpy.CalculateStatistics_management(c4Raster, "1", "1", "#")
            options['c4Var'] = get_raster_var(c4Raster)
        if options['c5RasterBase'] is not None:
            c5Raster = path.join(options['projectDir'],options['c5RasterBase']) 
            arcpy.CalculateStatistics_management(c5Raster, "1", "1", "#")
            options['c5Var'] = get_raster_var(c5Raster)
        
        descData = arcpy.Describe(c1Raster)
        spatialReference = descData.spatialReference
        options['cellSize'] = descData.MeanCellHeight             
        header = get_header(c1Raster)          
        descData = arcpy.Describe(c1Raster)
        arcpy.env.extent = descData.Extent    
        
        cumDistRaster = arcpy.sa.Con((arcpy.Raster(c1Raster) > 0),0,0)        
      
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
        partialResultsArray = None
#new
        breakFlag = False
        count = 0
        calcInBand = False
        for centerRow in range((options['blockSize']-1)/2,header['nrows'],options['blockSize']):
#new
            count += 1
            bandNum += 1 
            if options['startBand'] > 0 and bandNum < options['startBand']: continue
            if options['endBand'] > 0 and bandNum >= options['endBand']+1: 
                breakFlag = True
                break
            print 'Starting band #',str(bandNum)+'/'+approxEndBand,' centered on row '+str(centerRow)

            c1BandArray, LLC = band(c1Raster,header,centerRow, options)

            if npy.max(c1BandArray) == -9999: continue
            c2BandArray, dummy  = band(c2Raster, header, centerRow, options)
            c3BandArray, dummy  = band(c3Raster, header, centerRow, options)
            c4BandArray, dummy  = band(c4Raster, header, centerRow, options)
            c5BandArray, dummy  = band(c5Raster, header, centerRow, options)
#new
            if partialResultsArray is None:
# FIXME: create partialresultsarray that is SIZE of band, but not ALIGNED with band. OR, create array that is size of REMAINDER.            
# LLC is same
# size is 
# except for last band?
# that's when LLC equals the LLC of the input raster
#fixme: this doesn't work when using an end band!
                nResultsRows = options['radius']+1
                # if LLC 
                # nResultsRows > 
                partialResultsArray = npy.zeros((nResultsRows,header['ncols']), dtype = 'float64')-9999 
                # partialResultsCenterRow = centerRow
                partialResultsLLC = LLC 
                
                # print 'partialResultsArray.shape[0]'
                # print partialResultsArray.shape[0]

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
                # print 'Band #',bandNum,'Stripe #',stripeNum
    #fixme- check this:
                c2CircleVector = c3CircleVector = c4CircleVector = c5CircleVector = None
                subsetCenterCol = min(centerCol,options['radius']) 
                subsetCenterRow = min(centerRow,options['radius'])
                # time consuming- add smart search. 
                # if options['blockSize']>1:
                    # if npy.max(sourceCenterArraySum0[centerCol-(options['blockSize']-1)/2:centerCol+(options['blockSize']-1)/2+1]) <= 0: 
                        # print'No values in radius'
                        # continue
                # else:
                    # if sourceCenterArray[centerCol] <= 0: #BHM ADDED INDEX 10/12/15 
                        # print'No sources in radius'
                        # continue 
                # print c1BandArray
                # time consuming- add smart search. 
                c1CircleVector= circ(c1BandArray, rowArray, colArray, subsetCenterRow, centerCol, options)
                if c2BandArray is not None:
                    c2CircleVector = circ(c2BandArray, rowArray, colArray, subsetCenterRow, centerCol, options)
                # Was getting memory error calling circ with None arrays, setting to None above instead                    
                if c3BandArray is not None:
                    c3CircleVector = circ(c3BandArray, rowArray, colArray, subsetCenterRow, centerCol, options)
                if c4BandArray is not None:
                    c4CircleVector = circ(c4BandArray, rowArray, colArray, subsetCenterRow, centerCol, options)
                if c5BandArray is not None:
                    c5CircleVector = circ(c5BandArray, rowArray, colArray, subsetCenterRow, centerCol, options)
                # circleHeader = get_subset_header(c1CircleArray, header, options, centerRow, centerCol)
                # yMin = circleHeader['yllcorner'] #fixme- check. coudl be something like: max(circleHeader['yllcorner'],circleHeader['yllcorner'] + ((circleHeader['nrows'] - centerRow - options['radius'] - 1) * circleHeader['cellsize']))
                
            
                # print '\nDone with prep'
                # start_time0 = elapsed_time(start_time0)           

                # pctDone = report_pct_done(iter, maxNumSolvesToDo, -1)
                # print 'Elapsed time so far: {0}'.format(datetime.datetime.now()-theStart)  
                # print'subs'

                
                
                if c1BandArray[subsetCenterRow, centerCol]!=-9999:
                    calcInBand = True                        
                    stdEucDist, options = std_euc(c1CircleVector, c2CircleVector, c3CircleVector, c4CircleVector, c5CircleVector, subsetCenterRow, subsetCenterCol, options)

                    if options['blockSize']>1:
                        print'error- cannot do blocks yet'
                        exit(0)
                        partialResultsArray[(subsetCenterRow-(options['blockSize']-1)/2):(subsetCenterRow+(options['blockSize']-1)/2)+1,(centerCol-(options['blockSize']-1)/2):(centerCol+(options['blockSize']-1)/2)+1] = stdEucDist
                    else:
                        # partialResultsArray[subsetCenterRow,centerCol] = stdEucDist
                        partialResultsArray[count-1,centerCol] = stdEucDist

                        # print'subcr,cc,std'
                        # print subsetCenterRow,centerCol,stdEucDist
                        # print 'bda'
                        # print partialResultsArray
                        # print 'donebda'
                        
                        
                del c1CircleVector, c2CircleVector, c3CircleVector, c4CircleVector, c5CircleVector                              
                start_time1 = time.clock()     
            # print 'Band Dist array:'
            # print partialResultsArray
            print 'Done with band#',str(bandNum)+'/'+approxEndBand
            # pctDone = report_pct_done((bandNum-options['startBand']+1)*(stripeNum-options['startStripe']+1), maxNumSolvesToDo, -1)

            pctDone = report_pct_done(iter+1, maxNumSolvesToDo, -1)
            # if calcInBand: #
                # partialResultsArray = 
#new
            if count == partialResultsArray.shape[0]:
#new            
                #if we've filled up partialResultsArray, add it to the cumulativeraster
                # yMin = max(header['yllcorner'],header['yllcorner'] + ((header['nrows'] - partialResultsCenterRow - options['radius'] - 1) * header['cellsize']))
                # LLC = arcpy.Point(header['xllcorner'],yMin)
                # print 'pra shape, llc'
                
                # print partialResultsArray, partialResultsLLC
                partialResultsRaster = arcpy.NumPyArrayToRaster(partialResultsArray, partialResultsLLC, header['cellsize'],header['cellsize'],-9999)                          

                # # SAVING BAND RASTER REMOVES OCCASIONAL HORIZONTAL STRIPING
                tempBandFile = os.path.join(options['scratchDir'], 'justBAND'+str(bandNum)+'cur_' + options['outputFileText']+'.tif')
                partialResultsRaster.save(tempBandFile)
                delete_data(tempBandFile)
                
                cumDistRaster = addData_arcpy(cumDistRaster, partialResultsRaster)
                del partialResultsArray
                partialResultsArray = None
                calcInBand=False
                count = 0
                del partialResultsRaster
    
              
                options = write_temp_maps(options,cumDistRaster,bandNum) 
                                         
            # print 'Done with band #',bandNum,',
            print 'Elapsed time so far: {0}'.format(datetime.datetime.now()-theStart)  
        print 'DONE with bands'
        # print 'partialResultsArray, partialResultsLLC'    
        # print partialResultsArray, partialResultsLLC    
        # print 'partialres data'
        # print partialResultsArray[0:count,:]
        
        # LLC = arcpy.Point(header['xllcorner'],header['yllcorner'])
        if breakFlag:#there's an endband
            partialResultsRaster = arcpy.NumPyArrayToRaster(partialResultsArray, partialResultsLLC, header['cellsize'],header['cellsize'],-9999)                          
        else: 
            partialResultsRaster = arcpy.NumPyArrayToRaster(partialResultsArray[0:count,:], partialResultsLLC, header['cellsize'],header['cellsize'],-9999)                          
        # # SAVING BAND RASTER REMOVES OCCASIONAL HORIZONTAL STRIPING
        tempBandFile = os.path.join(options['scratchDir'], 'justBAND'+str(bandNum)+'cur_' + options['outputFileText']+'.tif')
        partialResultsRaster.save(tempBandFile)
        delete_data(tempBandFile)

        cumDistRaster = addData_arcpy(cumDistRaster, partialResultsRaster)
        del partialResultsRaster
        # cumDistRaster.save(r'c:\md\test.tif')


        print 'Done with solves.'  
        print_prof_data()
        write_final_maps(options,cumDistRaster) 
        #xprint locals()
        return options
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print c1BandArray.shape
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

@profile    
def std_euc(vector1, vector2, vector3, vector4, vector5, subsetCenterRow, subsetCenterCol,options):
    try:    

        Var =[options['c1Var']]
        sampleArray = npy.array([vector1])

        if options['numInputRasters'] == 5:
            sampleArray = npy.array([vector1, vector2, vector3, vector4, vector5])        
            Var =[options['c1Var'],options['c2Var'],options['c3Var'],options['c4Var'],options['c5Var']]
        elif options['numInputRasters'] == 4:
            sampleArray = npy.array([vector1, vector2, vector3, vector4])        
            Var =[options['c1Var'],options['c2Var'],options['c3Var'],options['c4Var']]
        elif options['numInputRasters'] == 3:
            sampleArray = npy.array([vector1, vector2, vector3])        
            Var =[options['c1Var'],options['c2Var'],options['c3Var']]    
        elif options['numInputRasters'] == 2:
            sampleArray = npy.array([vector1, vector2])        
            Var =[options['c1Var'],options['c2Var']]
        else:
            sampleArray = npy.array([vector1])        
            Var =[options['c1Var']]    
        sampleArray = sampleArray.transpose() 
        # try:

            # blarg
        try:
            distances = get_distances(sampleArray,Var)
        except:
            print "NODATA AREAS MAY NOT MATCH"
            print 'vector 1, vector 2'
            print vector1
            print vector2 
            print 'shapes'
            print vector1.shape
            print vector2.shape
            print 'var, sa'
            print Var
            print sampleArray
            print 'now shape',sampleArray.shape
            print 'previous'
            print options['prevVar']
            print options['prevSampleArray']
            blarg
        options['prevVar'] = Var
        options['prevSampleArray'] = sampleArray

        meanDist = distances.mean()

        return meanDist,options
        
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

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
def get_distances(sampleArray, Var):
    try:
        return pdist(sampleArray, 'seuclidean', V=Var)
    except:    
        print 'shape,',sampleArray.shape
        print_python_error()

    
    
# @profile
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
def get_raster_var(raster):
    try:
        stdObject = arcpy.GetRasterProperties_management(raster, "STD") 
        STD = float(str(stdObject.getOutput(0)))
        return STD*STD
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
    

def write_final_maps(options,cumDistRaster):
    try:
        # print 'Creating final output maps'  
        # quantCurrentRaster = quantVdiffRaster = quantFlowRaster = None
        # if options['quantilizeMaps']:
            # quantCurrentRaster = quantilize(cumCurrentRaster)
            # # quantVdiffRaster = quantilize(cumVdiffRaster) 
            # # quantFlowRaster = quantilize(cumFlowRaster)

            sumFile = path.join(options['outputDir'], 'dist_' + options['outputFileText']+'.tif')
            print 'Writing:\n',sumFile
            outRasterString = sumFile
            cumDistRaster.save(sumFile)
            del cumDistRaster
            # if options['cleanUpBandFiles']:
                # delete_data(options['prevCumCurrentFile'])

            print '\nSaved output files (you can just copy and paste entire path into ArcMap add data dialog):'
            print outRasterString,'\n'
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()
    
    
@profile    
def write_temp_maps(options,cumDistRaster,bandNum): 
    try:
        print 'Writing temporary grids...\n'
        sumFile = path.join(options['outputDir'], 'dist_PARTIAL' + options['outputFileText']+'_'+str(bandNum)+'.tif')
        outRasterString = sumFile
        # delete_data(sumFile)
        cumDistRaster.save(sumFile)
        delete_data(options['prevTempFile'])
        options['prevTempFile']=sumFile

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

# @profile
def set_options_and_dirs(options):
    try:
        """derives output base filename by combining input options. Also handles filenames for tiling. Checks and corrects options    """    
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
            
        

        options['scratchDirBase'] = 'scratch'+options['outputDirBase'] 
        if options['tileNum'] > 0:
            options['scratchDirBase'] = options['scratchDirBase'] + str(options['tileNum'])
        
                
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
        if options['c2RasterBase'] is None: 
            options['numInputRasters']=1
        else:
            resisRasText2, fileExtension = os.path.splitext(options['c2RasterBase'])
            resisRasText = resisRasText+'_'+resisRasText2
            if options['c3RasterBase'] is None: 
                options['numInputRasters']=2
            else:
                resisRasText3, fileExtension = os.path.splitext(options['c3RasterBase'])
                resisRasText = resisRasText+'_'+resisRasText3

                if options['c4RasterBase'] is None: 
                    options['numInputRasters']=3
                else:
                    resisRasText4, fileExtension = os.path.splitext(options['c4RasterBase'])
                    resisRasText = resisRasText+'_'+resisRasText4
                    if options['c5RasterBase'] is None: 
                        options['numInputRasters']=4
                    else:
                        resisRasText5, fileExtension = os.path.splitext(options['c5RasterBase'])
                        resisRasText = resisRasText+'_'+resisRasText5
                        options['numInputRasters']=5

        # if not options['useSourceRaster']:
        # if len(resisRasText)>25:
                # resisRasText = resisRasText[0:24]+'x'                    
        # elif len(resisRasText)>15:
            # resisRasText = resisRasText[0:14]+'x'                    


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
            
        radiusText = '_r'+str(options['radius'])
        # options['outputFileText'] = nullText+resisRasText + squareText + limText + '_' + srcText+climText+radiusText+'b'+str(options['blockSize'])+cutoffText +distFunctionText+weightText+noWeightText+centerText+maxText+startBandText+endBandText+startStripeText+endStripeText+negTargText+fadeInText
        options['outputFileText'] = resisRasText +radiusText+'b'+str(options['blockSize'])+startBandText+endBandText+startStripeText+endStripeText
        # options['prevFlowFile'] = None
        # options['prevCumCurrentFile'] = None
        # options['prevVdiffFile'] = None

        print 'c1RasterBase = ',options['c1RasterBase']
        print 'Radius = ',options['radius']
        print 'BlockSize = ',options['blockSize']
        print'startBand',options['startBand']
        print'endBand',options['endBand']
        print'startStripe',options['startStripe']
        print'endStripe',options['endStripe']

        options['outputDir'] = os.path.join(options['projectDir'],options['outputDirBase'])
        options['scratchDir'] = os.path.join(options['projectDir'],options['scratchDirBase'])
        print 'project Dir',options['projectDir']
        print 'output Dir',options['outputDirBase']
        delete_dir(options['scratchDir'])
        if not path.exists(options['scratchDir']):
            os.mkdir(options['scratchDir'])
        if not path.exists(options['outputDir']):
            os.mkdir(options['outputDir'])
        print '\n'
        options['prevTempFile']=None
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


@profile
def band(inRaster,header,centerRow, options): 
    try:
        if inRaster is None: return None, None
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
        # print 'LLC', LLC
        return bandArray, LLC
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
        
        ##Nodata value: 9
        # arcpy.Mosaic_management(currentRaster,cumCurrentRaster,"FIRST","","", "", "", "", "") #takes too long
        # return cumCurrentRaster
        newCumCurrentRaster = arcpy.sa.Con(arcpy.sa.IsNull(currentRaster),cumCurrentRaster,currentRaster)
        
        # newCumCurrentRaster2 = arcpy.sa.Con(arcpy.sa.IsNull(newCumCurrentRaster),currentRaster,newCumCurrentRaster)
        return newCumCurrentRaster
    except arcpy.ExecuteError:
        print_geoproc_error()
    except:    
        print_python_error()

@profile
def circ(array, rowArray, colArray, centerRow, centerCol, options):
    try:
        if array is None: return None

        colArraySmall, rowArraySmall, arraySmall = get_array_small(array, colArray, rowArray, centerRow, centerCol, options)
        
        distArray =  get_dist_array(rowArraySmall, centerRow, colArraySmall, centerCol)

        del rowArraySmall,colArraySmall
        # arrayMasked = npy.where(distArray <= options['radius'], arraySmall, -9999) #fixme do nans here
        vector = get_vector(distArray, options, arraySmall)
            
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

# @profile
def get_vector(distArray, options, arraySmall):
    rows,cols = npy.where((distArray <= options['radius'])&(arraySmall>-9999))
    vector=((arraySmall[rows,cols]).flatten()).transpose()
    return vector


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

        
# @profile
def get_array_small(array, colArray, rowArray, centerRow, centerCol, options):
    startRow = max(centerRow - options['radius'],0)
    endRow = min(centerRow + options['radius'],array.shape[0]-1)
    startCol = max(centerCol - options['radius'],0)
    endCol = min(centerCol + options['radius'],array.shape[1]-1)        
#DITCH THESE
    colArraySmall = (colArray[startRow:endRow+1,startCol:endCol+1]).astype('float64')
    rowArraySmall = (rowArray[startRow:endRow+1,startCol:endCol+1]).astype('float64')
    arraySmall = array[startRow:endRow+1,startCol:endCol+1]
    smallGrid = npy.indices((arraySmall.shape))
    baseRow = max(centerRow - options['radius'], 0)
    baseCol = max(centerCol - options['radius'], 0)
    
    rowArraySmall2 = smallGrid[0] + baseRow
    colArraySmall2 = smallGrid[1] + baseCol

    # print npy.sum(rowArraySmall2 - rowArraySmall)
    if npy.sum(rowArraySmall2 - rowArraySmall) != 0:
        blarg
    if npy.sum(colArraySmall2 - colArraySmall) !=0:
        print rowArraySmall
        print rowArraySmall2
        print colArraySmall
        print colArraySmall2
        print npy.sum(colArraySmall2 - colArraySmall)
        blarg2
    
        
    del smallGrid

    return colArraySmall, rowArraySmall, arraySmall
    
    
# @profile
def get_dist_array(rowArraySmall, centerRow, colArraySmall, centerCol):
    return npy.sqrt(npy.multiply((rowArraySmall - centerRow), (rowArraySmall- centerRow)) + npy.multiply((colArraySmall-centerCol), (colArraySmall-centerCol))) #causes memory error

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
    # clean_up(options)
    if tileNum>0:
        raw_input('Done. Press Enter to continue.')
    print 'Done'    

    
    
# # ---run with arcgis win32 python... that's the only one I have scipy setup for

# import numpy as np
# import scipy.spatial

# print scipy.__version__


# # IS THIS IT???
# # http://stackoverflow.com/questions/21003272/difference-between-all-1d-points-in-array-with-python-diff
# # a = npy.random.random(700)
# # print 'a0 a1 a2',a[0],a[1],a[2]
# # print a[0]-a[1]
# # print a[0]-a[2]
# # diff=npy.subtract.outer(a,a)[npy.tril_indices(a.shape[0],k=-1)]
# # print 'diff',diff
# # print diff.size

# Z1 = npy.random.random(700)
# Z2 = npy.random.random(700)
# Var1=npy.var(Z1)
# Var2=npy.var(Z2)

# diff1 = npy.subtract.outer(Z1,Z1)[npy.tril_indices(Z1.shape[0],k=-1)]
# diff2 = npy.subtract.outer(Z2,Z2)[npy.tril_indices(Z2.shape[0],k=-1)]
# diff1Sq = npy.multiply(diff1,diff1)
# diff2Sq = npy.multiply(diff2,diff2)
# diff1SqDivVar = npy.divide(diff1Sq,Var1)
# diff2SqDivVar = npy.divide(diff2Sq,Var2)

# SEucDist = npy.sqrt(diff1SqDivVar + diff2SqDivVar)
 
# print 'Manual calc'
# print SEucDist
# print 'Manual mean',npy.mean(SEucDist)


# Var = [Var1,Var2]
# Z=npy.zeros((700,2),dtype='float64')
# Z[:,0]=Z1
# Z[:,1]=Z2

# Y1 = scipy.spatial.distance.pdist(Z, 'seuclidean', V=Var)#=None)
# Y2 = scipy.spatial.distance.pdist(Z, 'seuclidean', V=None)#=None)
# print 'Y1',Y1
# print 'Y2',Y2
# print 'average Y1, Y2',npy.mean(Y1),npy.mean(Y2)
# print Y1.size

# # a = npy.array(range(5, 10))
# # b = npy.array(range(1, 6))

# # res = a[npy.newaxis, :] - b[:, npy.newaxis]
# # print 'a',a
# # print 'b',b
# # print(res)

# blarg

# # use numpy?
# # http://stackoverflow.com/questions/26076576/subtract-all-pairs-of-values-from-two-arrays


# # for Var1, get all (x1-y1)
# # Diffxy[x,y] = Var1[x]-Var1[y]


# # http://stackoverflow.com/questions/17936587/in-numpy-find-euclidean-distance-between-each-pair-from-two-arrays

# loc1=(npy.array([[1,2,3,4,5],[0,0,0,0],[0,0,0,0],[0,0,0,0]])).transpose()
# loc2=(npy.array([1,3,4,6],[0,0,0,0],[0,0,0,0],[0,0,0,0]])).transpose()
# loc3=(npy.array([[1,3,4,6],[0,0,0,0],[0,0,0,0],[0,0,0,0]])).transpose()

# test=loc1-loc2
# print loc1
# print loc2
# print v1
# print v3

# print 'test'
# print test
# print npy.multiply(test,test)
# blar
# Var = npy.array([1.0,3.0])
# # Var=Var.transpose()

# print v1-v2
# print npy.multiply((v1-v2),(v1-v2))
# distances = (v1-v2)^2
# distances = distances.sum(axis=-1)
# print 'dist',distances



# X=npy.zeros((2,4),dtype='int32')
# X[0,:]=v1
# X[1,:]=v2
# X = X.transpose()
# print X
# Y = scipy.spatial.distance.pdist(X, 'seuclidean', V=Var)#=None)

# print Y    




# import numpy as np
# from scipy.spatial.distance import euclidean, pdist, wminkowski, squareform

# # loc1=(npy.array([[1,3,4,5,6],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]])).transpose()
# # loc2=(npy.array([[1,3,4,6,6],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]])).transpose()
# # loc3=(npy.array([[1,3,4,6,6],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]])).transpose()
# # dist=scipy.spatial.distance.euclidean(loc1,loc2)
# s1=npy.array([0.132,-0.156,1.576,0])
# s2=npy.array([-0.802,0.036,-1.979,0])
# s3=npy.array([0.413,-0.988,-1.268,0])
# s4=npy.array([1.72,-0.668,-0.557,0])
# s5=npy.array([-0.288,-0.860,0.154,0])
# # list_of_locs = [s1, s2, s3, s4, s5]

# # make a 4x3 matrix from list of objects
# X = npy.array([s1, s2, s3, s4, s5])
# print X
# #calculate pairwise distances, using weighted Minkowski norm
# distances = pdist(X,wminkowski,2, [1,1,1,1])
# distances2 =pdist(X, 'seuclidean', V=[1.,1.,1.,1.])
# print distances2


# print distances-distances2

# #make a square matrix from result
# distances_as_2d_matrix = squareform(distances)

# print distances
# # Then take mean of distances to get final answer?
# meanDist = distances.mean()
# print meanDist

# # with variances entered separately

# s1=npy.array([4.8,72,3.5])
# s2=npy.array([2.8,75,2.5])
# s3=npy.array([5.4,59,2.7])
# s4=npy.array([8.2,64,2.9])
# s5=npy.array([3.9,61,3.1])
# # s4=npy.array([1.72,-0.668,-0.557,0])
# # s5=npy.array([-0.288,-0.860,0.154,0])

# Var = [2.141*2.141,15.615*15.615,0.281*0.281]
# X = npy.array([s1, s2, s3, s4, s5])
# distances3 =pdist(X, 'seuclidean', V=Var)
# print 'D3'
# print distances3
# print 'Mean'
# meanDist = distances3.mean()
# print meanDist

# # exit(0)

# print distances_as_2d_matrix
# exit(0)
# compute distance matrix
# dist = npy.zeros((3,3), dtype = 'float64')
# for i in range (0,3):
# for j in range(1,3):
# dist[i,j]=
# test=loc1-loc2
# print 'loc1',loc1
# print 'loc2',loc2
# print 'dif12',test

# loc1v1=(npy.array([[1],[0],[2],[5]])).transpose()
# loc2v1=(npy.array([[1],[2],[3],[4]])).transpose()
# print loc1v1
# print loc2v1
# dist=scipy.spatial.distance.euclidean(loc1v1,loc2v1)
# print dist
# exit(0)
# a1 = npy.array([1,2,npy.nan])
# a2 = npy.array([1,2,3])
# print npy.nansum(a1-a2)
#fixme: step 1 = create a nodata mask to set all inputs to nodata when 1 is nodata.




# 3200000	km2 (about 2/5 area of conus)
# 3.2E+12	m2
# 3555555556	30m pixels
# 3555.555556	
	
# 3.5 billion pixels	
# *250 thousand calcs	
	
# 8.88889E+14	calcs
# 888.8888889	
# 888 trillion calcs	

# load bands for all N rasters
# search band1
# when data pixel found, do circ for all rasters. nodata = nan, 

# BLOCKING- keeping code in for now, anticipate only running with BS=1

# def std_euc(array1, array2, array3, array4, subsetCenterRow, subsetCenterCol):
#RETURNS A SINGLE VALUE 
    # varC1 = options['varC1'] etc.
    # (Variances can perhaps be calculated on a SAMPLING of global values- or getrasterproperties)
    # C1 diffs, c2 diffs, etc...
    # 450m ~707 pixels= 250k comparisons
            # pretty quick- see manual SEucDist code at end.
    
    
    






