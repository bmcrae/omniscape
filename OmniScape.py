# Add get_output_basename(options) which does ALL text and returns a single string for use in all output files.

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
            # print test
        
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
        
        
        
        
        
        
tile=-1
#---------------------------------------------------------------------
# BASIC INPUTS #
#---------------------------------------------------------------------
options={}
options['projectDir'] = r'C:\Dropbox\Working\Dickson11States\PewOmniAnalysis072515' 
options['resisRasterBase'] = 'hmv8w_90x2s_plus1_Pow10_810m_NHD_slp.asc'
options['outputDirBase'] = 'test'

options['useSourceRaster']=True
options['sourceRasterBase'] = 'rDissolved_PAD1_3_ICUNI_IV_NLCS_gt5kac.tif'

options['useClimate']=False
options['tCutoff']=0.25
options['climateRasterBase'] ='TMEAN_NE_clip.tif'

# SPECIAL FUNCTIONS
options['calcNull']=False

options['radius'] = 100 # in PIXELS
options['blockSize'] = 55 #Odd number
options['rCutoff'] = 1
options['squareResistances']=False

# VOLTAGES
options['calcVoltages']=False
options['adjustVoltages']=False

# Targets and weighting #
options['centerGround']=True # doesnt' really change results, pretty much identical
options['negTargets']= False # negative sources at targets- blocks can work better with centerground, fade out, and no neg targs.
options['weightTargOnly'] = False # total current flow is ~ ntargs. may make sense if # dispersers limited
options['noWeight']=False
options['subtractSources'] = False

# FADE CURRENTS 
options['fadeIn']=True # linear decrease with distance to center. best with center ground
from math import *
options['fadeInDist']=20#options['radius']/2#math.sqrt(2)*(options['blockSize']/2)

# ANALYSIS WINDOW ####
options['startBand'] = 5 # bands are horizontal, with width equal to options['blockSize']. approx nrows/options['blockSize'] bands in a raster.
options['endBand'] = 6 # stops before processing this band
options['startStripe'] = 9#stripes are vertical, with width equal to options['blockSize']
options['endStripe']=15 # 0 to ignore this

# CLIMATE ####
options['absClim']=False # connect if ABS VAL of climate differs by tcutoff. Meant to help with fade.

if options['calcNull']:
    options['outputDirBase'] = options['outputDirBase']+'_NULL'
options['projectDirBase'] = 'scratch'+options['outputDirBase']
options['compress'] = False

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
    options = getOutputFilename(options)
    
    theStart = datetime.datetime.now()       
    options = makeDirs(options)
    copyThisFile(options)

    arcpy.env.scratchWorkspace = options['scratchDir']
    arcpy.env.Workspace = options['scratchDir']
    os.environ["TEMP"] = options['scratchDir']
    os.environ["TMP"] = options['scratchDir']    

    # Set raster paths and export to ascii if needed. FIXME: don't really need ascii, just convenient for header code for now
    resisRaster = path.join(options['projectDir'],options['resisRasterBase'])
    fileBase, fileExtension = os.path.splitext(resisRaster)   

    if options['useClimate']:
        climateRaster = path.join(options['projectDir'],options['climateRasterBase'])
    if options['useSourceRaster']:
        sourceRaster = path.join(options['projectDir'],options['sourceRasterBase'])
    else:
        sourceRaster = arcpy.sa.Con((arcpy.Raster(resisRaster) <= options['rCutoff']),1,0)
    header = get_header(resisRaster)
      
    descData=arcpy.Describe(resisRaster)
    arcpy.env.extent=descData.Extent    

    cumCurrentRaster = arcpy.sa.Con((arcpy.Raster(resisRaster) > 0),0,0)
    if options['calcVoltages']:
        cumVdiffRaster  = arcpy.sa.Con((arcpy.Raster(resisRaster) > 0),0,0)   
    
    prevSumFile = None
    
    iter = 0
    bandNum = 0
    prevCumCurrentFile=None
    prevVdiffFile = None
   

    start_time0 = time.clock()
    for centerRow in range((options['blockSize']-1)/2,header['nrows'],options['blockSize']):
        solveInBand=False
        bandNum += 1 
        if options['startBand'] > 0 and bandNum < options['startBand']: continue
        if options['endBand'] > 0 and bandNum >= options['endBand']: break

        print 'Starting band #',bandNum,' out of ',int(header['nrows']/options['blockSize']),' centered on row '+str(centerRow)

        bandArray = band(resisRaster,header,centerRow, options)
        if options['squareResistances']:
            bandArray=npy.multiply(bandArray,bandArray)

        

        sourceBandArray = band(sourceRaster,header,centerRow, options)
        if options['useClimate']:
            climateBandArray = band(climateRaster,header,centerRow, options)
#fixme: getting extra nodata row on top of climate band?

        cumCurrentArray = npy.zeros(bandArray.shape, dtype = 'float64') 
        cumVdiffArray = cumCurrentArray.copy()

        if npy.max(bandArray) == -9999: continue
        
        subsetCenterRow = min(centerRow,options['radius'])
        if options['blockSize']>1:
            bandCenterArray=bandArray[(subsetCenterRow-(options['blockSize']-1)/2):(subsetCenterRow+(options['blockSize']-1)/2),:]
        else:
            bandCenterArray=bandArray[subsetCenterRow,:]
        if npy.max(bandCenterArray) == -9999:       
            del bandCenterArray
            continue

        if npy.max(sourceBandArray) <= 0:       
            print 'no sources in band; continuing'
            continue
        
        if options['blockSize']>1:
            sourceCenterArray=sourceBandArray[(subsetCenterRow-(options['blockSize']-1)/2):(subsetCenterRow+(options['blockSize']-1)/2+1),:]
        else:
            sourceCenterArray=sourceBandArray[subsetCenterRow,:]

        if npy.max(sourceCenterArray) <=0:
            print 'no targets; continuing' #fixme: need to differently handle sources and targets
            del sourceCenterArray
            continue
        # del sourceCenterArray fixme: saving for now to quickly look for valid target areas 
        del bandCenterArray
            
        sourceCenterArraySum0 = npy.sum(npy.where(sourceCenterArray > 0, sourceCenterArray, 0), axis=0)
        del sourceCenterArray
        
        grid = npy.indices((bandArray.shape))
        rowArray = grid[0]
        colArray = grid[1]
        del grid
        stripeNum=0
        for centerCol in range((options['blockSize']-1)/2, header['ncols'],options['blockSize']):
            iter = iter + 1
            stripeNum=stripeNum+1
            if options['startStripe'] > 0 and stripeNum < options['startStripe']:
                # print'cont at stripe',stripeNum
                continue
            if options['endStripe'] > 0 and stripeNum >= options['endStripe']:
                # print'break at stripe',stripeNum
                break
                
            print 'Band #',bandNum,'Stripe #',stripeNum
#fixme- check this:
            subsetCenterCol = min(centerCol,options['radius']) 
            subsetCenterRow = min(centerRow,options['radius'])
        
            if npy.max(sourceCenterArraySum0[centerCol-(options['blockSize']-1)/2:centerCol+(options['blockSize']-1)/2+1]) <= 0:
                # print 'no targets in stripe#',stripeNum
                continue
            

        
            # time consuming- add smart search
            circleResisArray = circ(bandArray, rowArray, colArray, subsetCenterRow, centerCol, options)
            if options['useSourceRaster']:
                sourceArray= circ(sourceBandArray, rowArray, colArray, subsetCenterRow, centerCol, options)
                sourceArray = npy.where(sourceArray < 0, 0, sourceArray)    
                
                
            # stamp everything outside block to nodata
#FIXME not working with large radii
            # print 'cirlce shape',circleResisArray.shape          
            # print 'cirlce min,max',npy.min(circleResisArray),npy.max(circleResisArray)
            # print 'subsetcenterrow, col',subsetCenterRow, subsetCenterCol
            blockResisArray = center_block(circleResisArray, options, subsetCenterRow, subsetCenterCol)                         
            # print 'cirlce center val,block center val',circleResisArray[subsetCenterRow,subsetCenterCol],blockResisArray[subsetCenterRow,subsetCenterCol]

            if options['radius'] <6:
                print  'circ'
                print circleResisArray.astype('int32')
                print  'block'
                print blockResisArray.astype('int32')

            if options['useSourceRaster']:
                targetArray = center_block(sourceArray, options, subsetCenterRow, subsetCenterCol)                             
                targetArray = npy.where(targetArray <0, 0,targetArray) #fixme- just set nodata sources to zero earlier    
            else:
                targetArray = npy.where(blockResisArray <= options['rCutoff'], 1,0)    
                targetArray = npy.where(blockResisArray < 0, 0, targetArray)    

            targetSum = targetArray.sum()
            if targetSum == 0:
                continue
            
            if not options['useSourceRaster']:
                sourceArray = npy.where(circleResisArray <= options['rCutoff'], 1,0)    
                sourceArray = npy.where(circleResisArray < 0, 0, sourceArray)  

            sourceArray = npy.where(targetArray > 0, 0, sourceArray)    
            sourceSum = sourceArray.sum()
            if sourceSum == 0:
                print'ss0'
                continue

            
# test- small block in center
            if options['centerGround']:
                groundArray = npy.where(targetArray > 0, -9999, -9999)
                groundArray[subsetCenterRow, subsetCenterCol] = 0 
            else:
                groundArray = npy.where(targetArray > 0, 10000000, -9999)

# end test
# Idea- current injected ~ numtargs.
# But then mask out target block and put ground in center.
# vastly increases run time:
            # rows,cols=npy.where(targetArray > 0)
            # groundArray=npy.zeros(sourceArray.shape,dtype='int32')
            # groundArray[rows[0],cols[0]] = 100000

# EXPERIMENT
            if len(sourceArray) < 10:
                print 'raw sources and targets:'
                print sourceArray
                print targetArray
                
            if options['useClimate']:
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
                        continue
                    tDiffArray=npy.where(climateArray == -9999, 0, climateArray-tTarget) #target is cooler
                    if options['absClim']:
                        tCutoffArray=npy.where(abs(tDiffArray)>=options['tCutoff'],1,0)
                    else:
                        tCutoffArray=npy.where(tDiffArray>=options['tCutoff'],1,0)
                        # for now: i=1 if tdiff > 1, otherwise 0
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
                # print 'whoa... do we really want identical source and target sums?'
               
            if sourceSum == 0 or targetSum==0:
                
                # print'ss0 or ts0'
                continue
            #then normalizing. 1 amp injected, 1 amp taken out
            targetArray = -targetArray/(targetSum+0.0) 
            sourceArray = sourceArray/(sourceSum+0.0)

            if options['negTargets']:
                sourceArray += targetArray
            print 'sourcesum, targsum',sourceSum,targetSum
# EXPERIMENT2 #nope, still causes solver failures 
            # targetArray = targetArray*sourceSum
            # sourceArray = sourceArray*targetSum

            circleHeader = get_subset_header(sourceArray, header, options, centerRow, centerCol)


            sourceAsciiFile = path.join(options['scratchDir'], 'source_r'+str(centerRow) + 'c' +str(centerCol)+'iter'+str(iter)+'.asc') 
            ascii_grid_writer(sourceAsciiFile, sourceArray, circleHeader, options)

            groundAsciiFile = path.join(options['scratchDir'], 'ground_r'+str(centerRow) + 'c' +str(centerCol)+'iter'+str(iter)+'.asc')            
            ascii_grid_writer(groundAsciiFile, groundArray, circleHeader, options)

            outputFN = 'target_r'+str(centerRow) + 'c' +str(centerCol)+'.out'
            outputFile = path.join(options['scratchDir'], outputFN)
            
            resisAsciiFile = path.join(options['scratchDir'], 'resis_r'+str(centerRow) + 'c' +str(centerCol)+'iter'+str(iter)+'.asc') 

            csOptions = setCircuitscapeOptions(resisAsciiFile,sourceAsciiFile,groundAsciiFile,outputFile)

            if options['calcVoltages']:
                csOptions['write_volt_maps']=True

            if options['calcNull']:
                circleResisArray= npy.where(circleResisArray>0,1,circleResisArray)
            
            ascii_grid_writer(resisAsciiFile, circleResisArray, circleHeader, options)

            configFN = 'config_target_b'+str(centerRow) + 'c' +str(centerCol)+'.ini'
            outConfigFile = path.join(options['scratchDir'], configFN)
            writeCircuitscapeConfigFile(outConfigFile, csOptions)  
            
            print '\nDone with prep'
            start_time0=elapsed_time(start_time0)
            print 'Ready to solve block row and col', centerRow, ',', centerCol, 'in band#',bandNum,' out of ',int(header['nrows']/options['blockSize'])
            print 'Elapsed time so far: {0}'.format(datetime.datetime.now()-theStart)  
            solveInBand=True
            CSPATH = get_cs_path()                            
            start_time1 = time.clock()
            call_circuitscape(CSPATH, outConfigFile)
            start_time1 = elapsed_time(start_time1)
            start_time0 = time.clock()

            delete_data(configFN)
            configFN2 = 'target_r'+str(centerRow) + 'c' +str(centerCol)+'.ini'
            delete_data(path.join(options['scratchDir'], configFN2))

            curMap = path.join(options['scratchDir'], 'target_r'+str(centerRow) + 'c' +str(centerCol) + '_curmap.asc')
            
            if not path.exists(curMap):
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
                multiplier = sourceSum * targetSum 
            
            if options['subtractSources']:
                currentArray = multiplier * (ascii_grid_reader(curMap, header['nodata'], 'float64') - abs(sourceArray))
            else:
                currentArray = multiplier * (ascii_grid_reader(curMap, header['nodata'], 'float64')) 


#fixme- may not be compatible with subsources, weights, polygons, etc:
#fixme: dont' mult by dists, but dist/options['radius']?

            if options['fadeIn']:
#fixme: could do these once (at least when current array is centered in band)
                grid = npy.indices((currentArray.shape))
                subsetRowArray = grid[0]
                subsetColArray = grid[1]
                del grid

                subsetDistArray = npy.sqrt(npy.multiply(subsetRowArray - subsetCenterRow, subsetRowArray- subsetCenterRow) + npy.multiply(subsetColArray-subsetCenterCol, subsetColArray-subsetCenterCol))           
                
                fadeArray=npy.where(subsetDistArray<options['fadeInDist'],subsetDistArray/options['fadeInDist'],1)
                
                currentArray=npy.where(currentArray>0,npy.multiply(currentArray,fadeArray),currentArray)
                del subsetRowArray, subsetColArray

# fixme: do fadeArray as multiplier- just a single one for fade in and out
            




            
# fixme: use numpyarraytoraster
            # temp
            # currentRaster = path.join(options['scratchDir'],'fillblockcurrenttemp'+str(iter)+'.asc')
            # ascii_grid_writer(currentRaster, currentArray, circleHeader, options)
            #TEMP- multiplying by targetsum to get around solver failures



            if options['calcVoltages']:
                voltMapPath = path.join(options['scratchDir'], 'target_r'+str(centerRow) + 'c' +str(centerCol) + '_voltmap.asc')
                voltMap = (ascii_grid_reader(voltMapPath, header['nodata'], 'float64'))
                voltMap = npy.where(circleResisArray<0,npy.nan,voltMap)
                max_vdiff = get_max_vdiff(voltMap,circleResisArray,csOptions)
                max_vdiff=npy.where(npy.isnan(max_vdiff),0,max_vdiff)
                
                                                  

# temp fix for striping at edges
                max_vdiff[0,::]=0
                max_vdiff[:,0]=0
                max_vdiff[:,circleHeader['ncols']-1]=0
                max_vdiff[circleHeader['nrows']-1,:]=0
                # vDiffPath = path.join(options['scratchDir'], 'target_r'+str(centerRow) + 'c' +str(centerCol) + '_VDIFF.asc')
                # ascii_grid_writer(vDiffPath, max_vdiff, circleHeader, options)

                if options['fadeIn']:
                    max_vdiff = npy.where(max_vdiff>0,npy.multiply(max_vdiff,fadeArray),max_vdiff)
                    del fadeArray

                if options['weightTargOnly']:
                    max_vdiff = targetSum * max_vdiff
                elif options['noWeight']:
                    pass
                else:    
                    max_vdiff = sourceSum * targetSum * max_vdiff

                cumVdiffArray = addData(cumVdiffArray, max_vdiff, subsetCenterRow, centerCol, options)
              
            delete_data(curMap)
            
            delete_data(outConfigFile)
            delete_data(groundAsciiFile)
            delete_data(resisAsciiFile)
            delete_data(sourceAsciiFile)

            start_time1 = time.clock()

            # cumCurrentArray = addData(cumCurrentArray, currentArray, centerRow, centerCol, options['radius'])
            cumCurrentArray = addData(cumCurrentArray, currentArray, subsetCenterRow, centerCol, options)


# FIXME: do all raster stuff for each band, not each solve
            # cumCurrentRaster = addData_arcpy(cumCurrentRaster, currentRaster)
# FIXME: move file writing to band iteration


#put below in fn
        
        if solveInBand: #options, preCumCurrentFile,prevVdiffFile = writeTempMaps(options,bandNum,cumCurrentArray,prevCumCurrentFile,cumVdiffArray,prevVdiffFile,header)
            print 'Done with band #',bandNum,'.'
            print 'Writing temporary grid...'
            bandRows=min(options['radius']*2+1,options['radius']+centerRow+1)
            yMin= max(header['yllcorner'],header['yllcorner'] + ((header['nrows'] - centerRow - options['radius'] - 1) * header['cellsize']))
            LLC = arcpy.Point(header['xllcorner'],yMin)

            bandCurrentRaster = arcpy.NumPyArrayToRaster(cumCurrentArray,LLC,
                                             header['cellsize'],header['cellsize'],-9999)

            cumCurrentRaster = addData_arcpy(cumCurrentRaster, bandCurrentRaster)
            cumCurrentFile = os.path.join(options['outputDir'],'BAND'+str(bandNum)+options['resisRasterBase']+'_sumtemp.tif')
            try:
                cumCurrentRaster.save(cumCurrentFile)
            except:
                try:
                    print 'Error writing'
                    time.sleep(5)                
                    scratchDir2=os.path.join(options['projectDir'],options['projectDirBase']+'a')
                    if not path.exists(scratchDir2):
                        os.mkdir(scratchDir2)
                    arcpy.env.scratchWorkspace = scratchDir2
                    arcpy.env.Workspace = scratchDir2

                    cumCurrentFile = os.path.join(options['outputDir'],'BAND'+str(bandNum)+options['resisRasterBase']+'_sumtemp2.tif')
                    cumCurrentRaster.save(cumCurrentFile)
                    arcpy.env.scratchWorkspace = options['scratchDir']
                    arcpy.env.Workspace = options['scratchDir']

                except:
                    print 'Second error writing'
                    time.sleep(5)
                    scratchDir3=os.path.join(options['projectDir'],options['projectDirBase']+'b')
                    if not path.exists(scratchDir3):
                        os.mkdir(scratchDir3)
                    arcpy.env.scratchWorkspace = scratchDir3
                    arcpy.env.Workspace = scratchDir3
                    cumCurrentFile = os.path.join(options['outputDir'],'BAND'+str(bandNum)+options['resisRasterBase']+'_sumtemp3.tif')
                    cumCurrentRaster.save(cumCurrentFile)                       
                    arcpy.env.scratchWorkspace = options['scratchDir']
                    arcpy.env.Workspace = options['scratchDir']
                    
            delete_data(prevCumCurrentFile)
            prevCumCurrentFile = cumCurrentFile

            if options['calcVoltages']:
                bandVdiffRaster = arcpy.NumPyArrayToRaster(cumVdiffArray,LLC,header['cellsize'],header['cellsize'],-9999)
                cumVdiffRaster = addData_arcpy(cumVdiffRaster, bandVdiffRaster)
                cumVdiffFile = os.path.join(options['outputDir'],'BAND'+str(bandNum)+options['resisRasterBase']+'_vdiffSumtemp.tif')
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
                    
                delete_data(prevVdiffFile)
                
                prevVdiffFile = cumVdiffFile
        
        print 'Done with band #',bandNum,'.Elapsed time so far: {0}'.format(datetime.datetime.now()-theStart)  

    print 'solved thru centerrow',centerRow


    for fl in glob.glob(path.join(options['outputDir'],'*.xml')):
        try:
            os.remove(fl)
        except:
            pass
    for fl in glob.glob(path.join(options['outputDir'],'*.tfw')):
        try:
            os.remove(fl)
        except:
            pass


    sumFile = path.join(options['outputDir'], 'cur_sum_' + options['outputFileText'])
    print 'writing:\n',sumFile
    cumCurrentRaster.save(sumFile)
    delete_data(prevCumCurrentFile)
    if options['calcVoltages']:
        sumFile = path.join(options['outputDir'], 'vDiff_sum_' + options['outputFileText'])
        print 'writing:\n',sumFile
        cumVdiffRaster.save(sumFile)
        delete_data(prevVdiffFile)

    delete_dir(options['scratchDir'])
    for filename in glob.glob(os.path.join(options['outputDir'],'band*.*')) :
        try:    
            os.remove( filename )
        except:
            pass

def getOutputFilename(options):
    """derives output base filename by combining input options. Also handles filenames for tiling    """    
    
    if tile >=0:
        options['outputDirBase']=options['outputDirBase']+str(tile)
        options['projectDirBase']=options['projectDirBase']+str(tile)
        fileBase,ext=os.path.splitext(options['resisRasterBase'])
        options['resisRasterBase']=fileBase+str(tile)+ext
        if options['useClimate']:
            fileBase,ext=os.path.splitext(options['climateRasterBase'])
            options['climateRasterBase']=fileBase+str(tile)+ext
        if options['useSourceRaster']:
            fileBase,ext=os.path.splitext(options['sourceRasterBase'])    
            options['sourceRasterBase']=fileBase+str(tile)+ext
    
    resisRasText, fileExtension = os.path.splitext(options['resisRasterBase'])
    if len(resisRasText)>10:
        resisRasText=resisRasText[0:14]+'Trnc'                    

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

    if options['weightTargOnly']:
        weightText='targOnly'
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
        fileBase,ext=os.path.splitext(options['climateRasterBase'])
        if len(fileBase)>10:
            fileBase=fileBase[0:14]+'Trnc'
        climText='clim_tc_'+str(options['tCutoff']).replace('.','')+'_'+fileBase
        if options['absClim']:
            climText=climText+'_abs'
    else:
        climText=''
    if options['calcNull']:
        nullText='NULL_'
    else:
        nullText=''

    options['outputFileText'] = nullText+resisRasText + '__'+srcText+climText+radiusText+'bl'+str(options['blockSize'])+'_cutoff'+str(options['rCutoff']) +weightText+noWeightText+subtrText+centerText+startBandText+endBandText+startStripeText+endStripeText+negTargText+fadeInText+'.tif'
    print 'resisRasterBase=',options['resisRasterBase']
    print 'Radius=',options['radius']
    print 'BlockSize=',options['blockSize']
    print'Rcutoff=',options['rCutoff']
    print'startBand',options['startBand']
    print'endBand',options['endBand']
    print'startStripe',options['startStripe']
    print'endStripe',options['endStripe']
    print'calcNull',options['calcNull']
    print 'fadeIn', options['fadeIn']

    return options

def copyThisFile(options):
    # Save a copy of this file in output directory
    destFile=os.path.join(options['outputDir'],os.path.basename(sys.argv[0]))
    ft = tuple(time.localtime())
    timeNow = time.ctime()
    fileName = ('%s_%s_%s_%s%s_%s' % (ft[0], ft[1], ft[2], ft[3], ft[4], os.path.basename(sys.argv[0])))
    filePath = os.path.join(options['outputDir'],fileName)
    shutil.copyfile(sys.argv[0],filePath) 

def makeDirs(options):
    options['outputDir']=os.path.join(options['projectDir'],options['outputDirBase'])
    options['scratchDir']=os.path.join(options['projectDir'],options['projectDirBase'])
    if not path.exists(options['projectDir']):
        os.mkdir(options['projectDir'])
    if not path.exists(options['scratchDir']):
        os.mkdir(options['scratchDir'])
    if not path.exists(options['outputDir']):
        os.mkdir(options['outputDir'])
    return options

            
def get_max_vdiff(voltMap,circleResisArray,csOptions):
    voltMap_l, voltMap_r = get_horiz_neighbors(voltMap)
    voltMap_u, voltMap_d = get_vert_neighbors(voltMap)

    if options['adjustVoltages']:
        resis_l, resis_r = get_horiz_neighbors(circleResisArray)
        resis_u, resis_d = get_vert_neighbors(circleResisArray)
        print 'adjustVoltages code not finished, nuances with resistance calcs'
        blarg
        
    vdiff_N = get_vdiff(voltMap_d,voltMap_u)
    # vdiff_S = get_vdiff(voltMap,voltMap_d)
    vdiff_E = get_vdiff(voltMap_l,voltMap_r)
    # vdiff_W = get_vdiff(voltMap,voltMap_l)
    del voltMap_u, voltMap_d, voltMap_l, voltMap_r 
    max_vdiff = npy.fmax(vdiff_N,vdiff_E)
    # max_vdiff = npy.fmax(max_vdiff,vdiff_E)
    # max_vdiff = npy.fmax(max_vdiff,vdiff_W)
    del vdiff_N,vdiff_E
    if csOptions['connect_four_neighbors_only'] == False:
        voltMap_ul, voltMap_dr = get_diag1_neighbors(voltMap)
        voltMap_ur, voltMap_dl = get_diag2_neighbors(voltMap)
        vdiff_NE = get_vdiff(voltMap_dl,voltMap_ur)
        vdiff_SE = get_vdiff(voltMap_ul,voltMap_dr)
        # vdiff_SW = get_vdiff(voltMap_ur,voltMap_dl)
        # vdiff_NW = get_vdiff(voltMap_dr,voltMap_ul)
        max_vdiff = npy.fmax(max_vdiff,vdiff_NE)
        # max_vdiff = npy.fmax(max_vdiff,vdiff_SE)
        max_vdiff = npy.fmax(max_vdiff,vdiff_SE)
        # max_vdiff = npy.fmax(max_vdiff,vdiff_NW)
        del voltMap_ul, voltMap_dr, voltMap_ur, voltMap_dl, vdiff_NE,vdiff_SE#,vdiff_SW,vdiff_NW
    return max_vdiff
    
def get_vdiff(voltMap,voltMap_direction):
    vdiff_direction = npy.absolute(voltMap - voltMap_direction)
    vdiff_direction = npy.where(voltMap == npy.nan, 0, npy.where(voltMap_direction == npy.nan,0,vdiff_direction))
    return vdiff_direction

def get_horiz_neighbors(map):
    m = map.shape[0]
    n = map.shape[1]
    zeromap = npy.zeros(map.shape,dtype = 'float64') 
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

def addData(cumCurrentArray, currentArray, centerRow, centerCol, options):
    # print 'adding data for center row, col, rad'
    # print centerRow,centerCol,options['radius']
    minRow = max(0,centerRow-options['radius'])
    maxRow = min(minRow+cumCurrentArray.shape[0]-1, centerRow+options['radius'])
    minCol = max(0,centerCol-options['radius'])
    # fixme: check if next line needs mincol added like in maxrow line
    maxCol = min(cumCurrentArray.shape[1]-1, centerCol+options['radius'])
    # print 'minRow,maxRow,minCol,maxCol'
    # print minRow,maxRow,minCol,maxCol
    # print 'currentArray.shape',currentArray.shape
    # print 'cumCurrentArray.shape',cumCurrentArray.shape
    fullCurrentArray = npy.zeros(cumCurrentArray.shape,dtype='float64')
    fullCurrentArray[minRow:maxRow+1,minCol:maxCol+1]=currentArray
    cumCurrentArray += fullCurrentArray
    return cumCurrentArray

def addData_arcpy(cumCurrentRaster, currentRaster):
    descData=arcpy.Describe(cumCurrentRaster)
    arcpy.env.extent=descData.Extent    
    newCumCurrentRaster = arcpy.sa.Con(arcpy.sa.IsNull(currentRaster),cumCurrentRaster,arcpy.sa.Plus(currentRaster,cumCurrentRaster))
    return newCumCurrentRaster


def circ(array, rowArray, colArray, centerRow, centerCol, options):
        distArray = npy.sqrt(npy.multiply(rowArray - centerRow, rowArray- centerRow) + npy.multiply(colArray-centerCol, colArray-centerCol))           
        arrayMasked = npy.where(distArray <= options['radius'], array, -9999)
        del distArray
        startRow = max(centerRow - options['radius'],0)
        endRow = min(centerRow + options['radius'],array.shape[0]-1)
        startCol = max(centerCol - options['radius'],0)
        endCol = min(centerCol + options['radius'],array.shape[1]-1)
        circleArray = arrayMasked[startRow:endRow+1,startCol:endCol+1]
        del arrayMasked
        return circleArray

def center_block(array, options, centerRow, centerCol):
        startRow = centerRow - ((options['blockSize']-1)/2)
        endRow = centerRow + ((options['blockSize']-1)/2)
        startCol = centerCol - ((options['blockSize']-1)/2)
        endCol = centerCol + ((options['blockSize']-1)/2)
        blockArray = npy.zeros(array.shape, dtype = 'float64') - 9999    
        blockArray[startRow:endRow+1,startCol:endCol+1] = array[startRow:endRow+1,startCol:endCol+1]
        return blockArray            

def center_mask(array, options, centerRow, centerCol):
        print 'masking',centerRow,centerCol
        startRow = centerRow - ((options['blockSize']-1)/2)
        endRow = centerRow + ((options['blockSize']-1)/2)
        startCol = centerCol - ((options['blockSize']-1)/2)
        endCol = centerCol + ((options['blockSize']-1)/2)
        maskedArray = array
        maskedArray[startRow:endRow+1,startCol:endCol+1] = 0
        return maskedArray            

def get_center_block(array, options, centerRow, centerCol):
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
        secs = int(elapsed)
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
    newRaster.save(outRasterPath)
    return
        

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
    print('     Calling Circuitscape:')
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
            print("      " + output.replace("\r\n",""))                
    
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

if __name__ == '__main__':
    omniscape(options)


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
# print aa
# aaa = np.where(aa>N/2,1,0)
# print aaa
# sum = np.sum(aaa, axis=0)
# ind = np.argmax(sum==1)
# print 'looking for '+str(N/2)
# print ind
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