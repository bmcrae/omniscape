# Omniscape base code
    # inputs:
        # sources (strengths) 
        # targets (strengths)
        # resistance
        # block size
        # targonly
        # radius
        # distance function for source strengths
        # fade?
# Climate calling code
    # inputs:
        # temp, or pc1, pc2 at T1, T2
        # for each temp or pc1/pc2 combo:
            # identify sources and targets
            # call omniscape_base 

# Partial out fade? Would get 3 outputs- raw, fade, raw-fade. Or just fade and raw-fade.

#4/16/15:
# removed polygonblocks
# remoed maskblocks
# removed fillblocks
# remove noweight?- doesn't make sense, since would produce diff results with diff block sizes.

# sourceRaster
# useSourceRasterCutoff #binary 1/0 if above(?) this value or have negative values to denote below this value
    # srcRasCutoffVal=

# INDEP SOURCES AND TARGS    
# If not useSourceRaster or useTargetRaster:
    # sourceraster, targraster = everything < sRCutoff, tRCutoff in resisRaster
    # if useSourceRaster and not useTargetRaster:
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
        
        
        
        
        
        
        
tile=-1
#---------------------------------------------------------------------
# BASIC INPUTS #
#---------------------------------------------------------------------

projectDir = r'C:\Dropbox\Working\Dickson11States\PewOmniAnalysis072515' 
resisRasterBase = 'hmv8w_90x2s_plus1_Pow10_810m_NHD_slp.asc'
outputDirBase = 'Omni_pow10'

useSourceRaster=True
sourceRasterBase = 'rDissolved_PAD1_3_ICUNI_IV_NLCS_gt5kac.tif'

useClimate=False
tCutoff=0.25
climateRasterBase ='TMEAN_NE_clip.tif'

# SPECIAL FUNCTIONS
calcNull=False

radius = 10000 # in PIXELS
blockSize = 55 #Odd number
rCutoff = 1
squareResistances=False
straddleBlocks=False

# VOLTAGES
calcVoltages=False
adjustVoltages=False

# Targets and weighting #
centerGround=True # doesnt' really change results, pretty much identical
negTargets= False # negative sources at targets- blocks can work better with centerground, fade out, and no neg targs.
weightTargOnly = True # total current flow is ~ ntargs. may make sense if # dispersers limited
noWeight=False
logTrans=False

# FADE CURRENTS 
fadeIn=True # linear decrease with distance to center. best with center ground
import math
fadeInDist=20#radius/2#math.sqrt(2)*(blockSize/2)
fadeOut=False
fadeOutDist = radius-15#(2.0/3)*(radius-(math.sqrt(2)*(blockSize/2)))+(math.sqrt(2)*blockSize/2)


# FADE VOLTAGES
fadeVoltages=False #multiply voltages by distance from center independently of current fading
bufferDist=0 # number pixels at outer edge NOT to have current sources


# ANALYSIS WINDOW ####
startBand = 0 # bands are horizontal, with width equal to blockSize. approx nrows/blockSize bands in a raster.
stopBand = 0
startStripe = 0 #stripes are vertical, with width equal to blockSize
stopStripe=0 # 0 to ignore this

# CLIMATE ####
absClim=False # connect if ABS VAL of climate differs by tcutoff. Meant to help with fade.

# DONUT ####
donutRadius = -95
maskDonut=-100

if calcNull:
    outputDirBase = outputDirBase+'_NULL'
scratchDirBase = 'scratch'+outputDirBase


import os 
if tile >=0:
    outputDirBase=outputDirBase+str(tile)
    scratchDirBase=scratchDirBase+str(tile)
    fileBase,ext=os.path.splitext(resisRasterBase)
    resisRasterBase=fileBase+str(tile)+ext
    if useClimate:
        fileBase,ext=os.path.splitext(climateRasterBase)
        climateRasterBase=fileBase+str(tile)+ext
    if useSourceRaster:
        fileBase,ext=os.path.splitext(sourceRasterBase)    
        sourceRasterBase=fileBase+str(tile)+ext

voltText=''
radiusText = '_r'+str(radius)
if useSourceRaster:
    fileBase,ext=os.path.splitext(sourceRasterBase)
    srcText = 'srcRas_'+ fileBase + '_'
else:
    srcText = ''
if negTargets:
    negTargText = 'negTarg'
else:
    negTargText = ''

if bufferDist>0:
    print 'bufferDist=',bufferDist
    bufferText='buff'+str(bufferDist)
else:
    bufferText=''
if fadeIn:
    fadeInText='fadeIn'+str(int(fadeInDist))
else:
    fadeInText=''
if fadeOut:
    fadeOutText='fadeOut'+str(int(fadeOutDist))
else:
    fadeOutText=''
if fadeVoltages:
    fadeVoltText='fadeVolt'
elif fadeIn or fadeOut:
    fadeVoltText=fadeInText+fadeOutText
else:
    fadeVoltText=''
if maskDonut>0:
    maskDonutText='maskDonut'+str(maskDonut)
else:
    maskDonutText=''
if straddleBlocks:
    straddleText='strad'
else:
    straddleText=''

if weightTargOnly:
    weightText='targOnly'
else:
    weightText=''
if centerGround:
    centerText='cg'
else:
    centerText=''
subtractSources = False

if noWeight:
    noWeightText='noWeight'
else:
    noWeightText=''
compress = False
if donutRadius>0:
   donutText = 'donut'+str(donutRadius)
else:
    donutText=''
if subtractSources:
    subtrText = 'subSrc'            
else:
    subtrText = ''

if useClimate:
    fileBase,ext=os.path.splitext(climateRasterBase)
    if len(fileBase)>10:
        fileBase=fileBase[0:14]+'Trnc'
    climText='clim_tc_'+str(tCutoff).replace('.','')+'_'+fileBase
    if absClim:
        climText=climText+'_abs'

else:
    climText=''
if logTrans:
    logText='LN_'
else:
    logText=''
if calcNull:
    nullText='NULL_'
else:
    nullText=''
       
print 'resisRasterBase=',resisRasterBase
print 'Radius=',radius
print 'BlockSize=',blockSize
print'Rcutoff=',rCutoff
if donutRadius > 0:
    print 'donutRadius=',donutRadius


print'startBand',startBand
print'stopBand',stopBand
print'startStripe',startStripe
print'stopStripe',stopStripe
print'calcNull',calcNull
print 'fadeIn', fadeIn


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
def circleFlowBand():

    theStart = datetime.datetime.now()
    outputDir=os.path.join(projectDir,outputDirBase)
    scratchDir=os.path.join(projectDir,scratchDirBase)
    if not path.exists(projectDir):
        os.mkdir(projectDir)
    if not path.exists(scratchDir):
        os.mkdir(scratchDir)
    if not path.exists(outputDir):
        os.mkdir(outputDir)
    
    # Save a copy of this file in output directory
    destFile=os.path.join(outputDir,os.path.basename(sys.argv[0]))
    ft = tuple(time.localtime())
    timeNow = time.ctime()
    fileName = ('%s_%s_%s_%s%s_%s' % (ft[0], ft[1], ft[2], ft[3], ft[4], os.path.basename(sys.argv[0])))
    filePath = os.path.join(outputDir,fileName)
    shutil.copyfile(sys.argv[0],filePath) 

    arcpy.env.scratchWorkspace = scratchDir
    arcpy.env.Workspace = scratchDir
    os.environ["TEMP"] = scratchDir
    os.environ["TMP"] = scratchDir    

    # Set raster paths and export to ascii if needed. FIXME: don't really need ascii, just convenient for header code for now
    resisRaster = path.join(projectDir,resisRasterBase)
    fileBase, fileExtension = os.path.splitext(resisRaster)   
    # if fileExtension.lower() != '.asc':
        # print '\nExporting resistance raster to ASCII. FIXME: can bypass this with new header code.'
        # asciiResisRaster=path.join(projectDir,fileBase+'_export.asc')
        # arcpy.RasterToASCII_conversion(resisRaster, asciiResisRaster)
        # # resisRaster = asciiResisRaster
    # else:
        # asciiResisRaster = resisRaster
    if useClimate:
        climateRaster = path.join(projectDir,climateRasterBase)

        
        
        # fileBase, fileExtension = os.path.splitext(climateRaster)   
        # if fileExtension.lower() != 'asc':
            # asciiClimateRaster=path.join(projectDir,fileBase+'_export.asc')
            # arcpy.RasterToASCII_conversion(climateRaster, asciiClimateRaster)
            # climateRaster = asciiClimateRaster
    if useSourceRaster:
        sourceRaster = path.join(projectDir,sourceRasterBase)
        # fileBase, fileExtension = os.path.splitext(sourceRaster)   
        # if fileExtension.lower() != 'asc':
            # asciiSourceRaster=path.join(projectDir,fileBase+'_export.asc')
            # arcpy.RasterToASCII_conversion(sourceRaster, asciiSourceRaster)
            # sourceRaster = asciiSourceRaster
    else:
        sourceRaster = arcpy.sa.Con((arcpy.Raster(resisRaster) <= rCutoff),1,0)
    header = get_header(resisRaster)
      
    descData=arcpy.Describe(resisRaster)
    arcpy.env.extent=descData.Extent    

    cumCurrentRaster = arcpy.sa.Con((arcpy.Raster(resisRaster) > 0),0,0)
    if calcVoltages:
        cumVdiffRaster  = arcpy.sa.Con((arcpy.Raster(resisRaster) > 0),0,0)   
    
    prevSumFile = None
    
    
    iter = 0
    bandNum = 0
    prevCumCurrentFile=None
    prevVdiffFile = None
    print'startBand',startBand
    print'stopBand',stopBand
    print'startStripe',startStripe
    print'stopStripe',stopStripe
    
    if straddleBlocks:
        step=blockSize/2
    else:
        step=blockSize
    start_time0 = time.clock()
    for centerRow in range((blockSize-1)/2,header['nrows'],step):
        solveInBand=False
        bandNum = bandNum + 1
        skipBand = False
        if startBand > 0 and bandNum < startBand:
            continue
        if stopBand > 0 and bandNum >= stopBand:
            print 'breaking'
            break
        print 'Starting band #',bandNum,' out of ',int(header['nrows']/blockSize),' centered on row '+str(centerRow)

        bandArray = band(resisRaster,header,centerRow, radius)
        if squareResistances:
            bandArray=npy.multiply(bandArray,bandArray)

        
#TEMP   
        # bandArray=npy.where(bandArray<0,100000,bandArray)
        # if useSourceRaster:
        sourceBandArray = band(sourceRaster,header,centerRow, radius)
        # print 'centerRow',centerRow
        # print 'radius',radius
        # print 'sourceBandArray'
        # print sourceBandArray
        if useClimate:
            climateBandArray = band(climateRaster,header,centerRow, radius)
#fixme: getting extra nodata row on top of climate band?

        cumCurrentArray = npy.zeros(bandArray.shape, dtype = 'float64') 
        cumVdiffArray = cumCurrentArray.copy()
# temp
        # bandCurrentArray = npy.zeros(bandArray.shape, dtype = 'float64') 
# temp
        if npy.max(bandArray) == -9999:
            print 'nodata band; continuing'
            continue
        subsetCenterRow = min(centerRow,radius)
        if blockSize>1:
            bandCenterArray=bandArray[(subsetCenterRow-(blockSize-1)/2):(subsetCenterRow+(blockSize-1)/2),:]
        else:
            bandCenterArray=bandArray[subsetCenterRow,:]
        if npy.max(bandCenterArray) == -9999:       
            print 'nodata center band; continuing' 
            del bandCenterArray
            continue

        # if useSourceRaster:
        if npy.max(sourceBandArray) <= 0:       
            print 'no sources in band; continuing'
            continue
        
        if blockSize>1:
            sourceCenterArray=sourceBandArray[(subsetCenterRow-(blockSize-1)/2):(subsetCenterRow+(blockSize-1)/2+1),:]
        else:
            sourceCenterArray=sourceBandArray[subsetCenterRow,:]

        if npy.max(sourceCenterArray) <=0:
            print 'no targets; continuing' #fixme: need to differently handle sources and targets
            del sourceCenterArray
            continue
        # del sourceCenterArray fixme: saving for now to quickly look for valid target areas 
        del bandCenterArray
            
        # sum of sources across band
        # print 'sourceCenterArray'
        # print sourceCenterArray
        sourceCenterArraySum0 = npy.sum(npy.where(sourceCenterArray > 0, sourceCenterArray, 0), axis=0)
        # print 'sum0'
        # print sourceCenterArraySum0
        del sourceCenterArray
        
        grid = npy.indices((bandArray.shape))
        rowArray = grid[0]
        colArray = grid[1]
        del grid
        stripeNum=0
        for centerCol in range((blockSize-1)/2, header['ncols'],step):
            iter = iter + 1
            stripeNum=stripeNum+1
            if startStripe > 0 and stripeNum < startStripe:
                # print'cont at stripe',stripeNum
                continue
            if stopStripe > 0 and stripeNum >= stopStripe:
                # print'break at stripe',stripeNum
                break
                
            print 'Band #',bandNum,'Stripe #',stripeNum
#fixme- check this:
            subsetCenterCol = min(centerCol,radius) 
            subsetCenterRow = min(centerRow,radius)
        
            # print 'source0sum for stripe:'
            # print sourceCenterArraySum0[centerCol-(blockSize-1)/2:centerCol+(blockSize-1)/2+1]
            if npy.max(sourceCenterArraySum0[centerCol-(blockSize-1)/2:centerCol+(blockSize-1)/2+1]) <= 0:
                # print 'no targets in stripe#',stripeNum
                continue
            

        
            # time consuming- add smart search
            circleResisArray = circ(bandArray, rowArray, colArray, subsetCenterRow, centerCol, radius)
            if useSourceRaster:
                sourceArray= circ(sourceBandArray, rowArray, colArray, subsetCenterRow, centerCol, radius)
                sourceArray = npy.where(sourceArray < 0, 0, sourceArray)    
                
                
            # stamp everything outside block to nodata
#FIXME not working with large radii
            # print 'cirlce shape',circleResisArray.shape          
            # print 'cirlce min,max',npy.min(circleResisArray),npy.max(circleResisArray)
            # print 'subsetcenterrow, col',subsetCenterRow, subsetCenterCol
            blockResisArray = center_block(circleResisArray, blockSize, subsetCenterRow, subsetCenterCol)                         
            # print 'cirlce center val,block center val',circleResisArray[subsetCenterRow,subsetCenterCol],blockResisArray[subsetCenterRow,subsetCenterCol]

            if radius <6:
                print  'circ'
                print circleResisArray.astype('int32')
                print  'block'
                print blockResisArray.astype('int32')

            if useSourceRaster:
                targetArray = center_block(sourceArray, blockSize, subsetCenterRow, subsetCenterCol)                             
                targetArray = npy.where(targetArray <0, 0,targetArray) #fixme- just set nodata sources to zero earlier    
            else:
                targetArray = npy.where(blockResisArray <= rCutoff, 1,0)    
                targetArray = npy.where(blockResisArray < 0, 0, targetArray)    

            targetSum = targetArray.sum()
            if targetSum == 0:
                # print 'ts0'
                continue
            
            if not useSourceRaster:
                sourceArray = npy.where(circleResisArray <= rCutoff, 1,0)    
                if bufferDist > 0:
                    grid = npy.indices((circleResisArray.shape))
                    subsetRowArray = grid[0]
                    subsetColArray = grid[1]
                    del grid
                    subsetDistArray = npy.sqrt(npy.multiply(subsetRowArray - subsetCenterRow, subsetRowArray- subsetCenterRow) + npy.multiply(subsetColArray-subsetCenterCol, subsetColArray-subsetCenterCol))           
                    sourceArray = npy.where((radius-subsetDistArray)<bufferDist,0,sourceArray)                

                sourceArray = npy.where(circleResisArray < 0, 0, sourceArray)  

            if donutRadius>0:
                sourceArray=donut(sourceArray, donutRadius, subsetCenterRow, subsetCenterCol)

            sourceArray = npy.where(targetArray > 0, 0, sourceArray)    
            sourceSum = sourceArray.sum()
            if sourceSum == 0:
                print'ss0'
                continue
            # print 'center',centerRow,centerCol
            # print 'subset center', subsetCenterRow, subsetCenterCol

            
# test- small block in center
            if centerGround:
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
                
            if useClimate:
                climateArray= circ(climateBandArray, rowArray, colArray, subsetCenterRow, centerCol, radius)              
                sourceArrayNew=npy.zeros(sourceArray.shape,dtype='float64')
                
                # indices of valid targets
                tRows, tCols=npy.where(targetArray)
                if len(sourceArray) < 10:
                    print 'clim'
                    print climateArray
                targetTally=0
                for i in range(0,len(tRows)):
                    # print 'tRows[i]'
                    # print tRows[i]
                    # print 'targ'
                    # print targetArray[tRows[i],tCols[i]]
                    # each target point
                    tTarget=climateArray[tRows[i],tCols[i]] 
                    if tTarget==-9999:
                        continue
                    tDiffArray=npy.where(climateArray == -9999, 0, climateArray-tTarget) #target is cooler
                    if absClim:
                        tCutoffArray=npy.where(abs(tDiffArray)>=tCutoff,1,0)
                    else:
                        tCutoffArray=npy.where(tDiffArray>=tCutoff,1,0)
                        # for now: i=1 if tdiff > 1, otherwise 0
                    # print 'tTarget',tTarget
                    if len(sourceArray) < 10:
                        print 'tCutoffArray'
                        print tCutoffArray
                    sourceArrayTarget_i=npy.multiply(sourceArray,tCutoffArray)# sources for this target cell
                    if npy.max(sourceArrayTarget_i)>0:
                        targetTally+=1 # There's a valid source for this target. Increment target count for scaling current in weightTargOnly setting
                    
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

            if negTargets:
                sourceArray += targetArray
            print 'sourcesum, targsum',sourceSum,targetSum
# EXPERIMENT2 #nope, still causes solver failures 
            # targetArray = targetArray*sourceSum
            # sourceArray = sourceArray*targetSum

            circleHeader = get_subset_header(sourceArray, header, radius, centerRow, centerCol)


            sourceAsciiFile = path.join(scratchDir, 'source_r'+str(centerRow) + 'c' +str(centerCol)+'iter'+str(iter)+'.asc') 
            ascii_grid_writer(sourceAsciiFile, sourceArray, circleHeader, compress)

            groundAsciiFile = path.join(scratchDir, 'ground_r'+str(centerRow) + 'c' +str(centerCol)+'iter'+str(iter)+'.asc')            
            ascii_grid_writer(groundAsciiFile, groundArray, circleHeader, compress)

            outputFN = 'target_r'+str(centerRow) + 'c' +str(centerCol)+'.out'
            outputFile = path.join(scratchDir, outputFN)
            
            resisAsciiFile = path.join(scratchDir, 'resis_r'+str(centerRow) + 'c' +str(centerCol)+'iter'+str(iter)+'.asc') 

            options = setCircuitscapeOptions(resisAsciiFile,sourceAsciiFile,groundAsciiFile,outputFile)

            if calcVoltages:
                options['write_volt_maps']=True

            if calcNull:
                circleResisArray= npy.where(circleResisArray>0,1,circleResisArray)
            
            ascii_grid_writer(resisAsciiFile, circleResisArray, circleHeader, compress)

            configFN = 'config_target_b'+str(centerRow) + 'c' +str(centerCol)+'.ini'
            outConfigFile = path.join(scratchDir, configFN)
            writeCircuitscapeConfigFile(outConfigFile, options)  
            
            print '\nDone with prep'
            start_time0=elapsed_time(start_time0)
            print 'Ready to solve block row and col', centerRow, ',', centerCol, 'in band#',bandNum,' out of ',int(header['nrows']/blockSize)
            print 'Elapsed time so far: {0}'.format(datetime.datetime.now()-theStart)  
            solveInBand=True
            CSPATH = get_cs_path()                            
            start_time1 = time.clock()
            call_circuitscape(CSPATH, outConfigFile)
            start_time1 = elapsed_time(start_time1)
            start_time0 = time.clock()

            delete_data(configFN)
            configFN2 = 'target_r'+str(centerRow) + 'c' +str(centerCol)+'.ini'
            delete_data(path.join(scratchDir, configFN2))

            curMap = path.join(scratchDir, 'target_r'+str(centerRow) + 'c' +str(centerCol) + '_curmap.asc')
            
            if not path.exists(curMap):
                print "Can't find circuitscape output"
                exit(0)
                continue
            if weightTargOnly:
                if useClimate:
                    multiplier = targetTally
                else:
                    multiplier = targetSum
            elif useClimate:
                multiplier = targetSum
            elif noWeight:
                multiplier = 1
            else:    
                multiplier = sourceSum * targetSum 
            
            if logTrans and not noWeight:
                multiplier = math.log(multiplier+1)

            if subtractSources:
                currentArray = multiplier * (ascii_grid_reader(curMap, header['nodata'], 'float64') - abs(sourceArray))
            else:
                currentArray = multiplier * (ascii_grid_reader(curMap, header['nodata'], 'float64')) 


#fixme- may not be compatible with subsources, weights, polygons, etc:
#fixme: dont' mult by dists, but dist/radius?

            if fadeIn or fadeOut:
#fixme: could do these once (at least when current array is centered in band)
                grid = npy.indices((currentArray.shape))
                subsetRowArray = grid[0]
                subsetColArray = grid[1]
                del grid

                subsetDistArray = npy.sqrt(npy.multiply(subsetRowArray - subsetCenterRow, subsetRowArray- subsetCenterRow) + npy.multiply(subsetColArray-subsetCenterCol, subsetColArray-subsetCenterCol))           
                
                if fadeIn:
                    fadeArray=npy.where(subsetDistArray<fadeInDist,subsetDistArray/fadeInDist,1)
                    if fadeOut:
                        fadeArray=npy.where(subsetDistArray>fadeOutDist,(radius-subsetDistArray)/(radius-fadeOutDist),fadeArray)
                else:
                    fadeArray=npy.where(subsetDistArray>fadeOutDist,(radius-subsetDistArray)/(radius-fadeOutDist),1)
                
                currentArray=npy.where(currentArray>0,npy.multiply(currentArray,fadeArray),currentArray)
                del subsetRowArray, subsetColArray

            # if fadeOut:
                # fadeOutArray=npy.where(subsetDistArray>fadeOutDist,radius-subsetDistArray,1)
                # currentArray = npy.where(currentArray>0,npy.multiply(currentArray,fadeOutArray),currentArray)
                # 
                # ascii_grid_writer('c:\\temp\\curFade.asc', currentArray, circleHeader, compress)
                # blarg
# fixme: do fadeArray as multiplier- just a single one for fade in and out
            




            
            if maskDonut:
                currentArray=donut(currentArray, maskDonut, subsetCenterRow, subsetCenterCol)
# fixme: use numpyarraytoraster
            # temp
            # currentRaster = path.join(scratchDir,'fillblockcurrenttemp'+str(iter)+'.asc')
            # ascii_grid_writer(currentRaster, currentArray, circleHeader, compress)
            #TEMP- multiplying by targetsum to get around solver failures



            if calcVoltages:
                voltMapPath = path.join(scratchDir, 'target_r'+str(centerRow) + 'c' +str(centerCol) + '_voltmap.asc')
                voltMap = (ascii_grid_reader(voltMapPath, header['nodata'], 'float64'))
                voltMap = npy.where(circleResisArray<0,npy.nan,voltMap)
                max_vdiff = get_max_vdiff(voltMap,circleResisArray,options)
                max_vdiff=npy.where(npy.isnan(max_vdiff),0,max_vdiff)
                
                                                  

# temp fix for striping at edges
                max_vdiff[0,::]=0
                max_vdiff[:,0]=0
                max_vdiff[:,circleHeader['ncols']-1]=0
                max_vdiff[circleHeader['nrows']-1,:]=0
                # vDiffPath = path.join(scratchDir, 'target_r'+str(centerRow) + 'c' +str(centerCol) + '_VDIFF.asc')
                # ascii_grid_writer(vDiffPath, max_vdiff, circleHeader, compress)

                if fadeVoltages:
                    max_vdiff = npy.where(max_vdiff>0,npy.multiply(max_vdiff,subsetDistArray),max_vdiff)        
                    del fadeArray
                elif fadeIn or fadeOut:
                    max_vdiff = npy.where(max_vdiff>0,npy.multiply(max_vdiff,fadeArray),max_vdiff)
                    del fadeArray

                if weightTargOnly:
                    max_vdiff = targetSum * max_vdiff
                elif noWeight:
                    pass
                else:    
                    max_vdiff = sourceSum * targetSum * max_vdiff

                cumVdiffArray = addData(cumVdiffArray, max_vdiff, subsetCenterRow, centerCol, radius)
  
                # delete_data(voltMapPath)
            
            delete_data(curMap)
            
            delete_data(outConfigFile)
            delete_data(groundAsciiFile)
            delete_data(resisAsciiFile)
            delete_data(sourceAsciiFile)

            start_time1 = time.clock()

            # cumCurrentArray = addData(cumCurrentArray, currentArray, centerRow, centerCol, radius)
            cumCurrentArray = addData(cumCurrentArray, currentArray, subsetCenterRow, centerCol, radius)

#TEMP
            # bandCurrentArray = addData(bandCurrentArray, currentArray, subsetCenterRow, centerCol, radius)
#TEMP




# FIXME: do all raster stuff for each band, not each solve
            # cumCurrentRaster = addData_arcpy(cumCurrentRaster, currentRaster)
# FIXME: move file writing to band iteration
            # sumFile = os.path.join(scratchDir,'sumtemp'+str(iter)+'.tif')
            # cumCurrentRaster.save(sumFile)
            
            # delete_data(currentRaster)
            
            # sumFile = path.join(outputDir, 'curmap_tempsum_r'+str(centerRow) + 'c' +str(centerCol) + '_cutoff'+str(rCutoff) + 'block'+str(blockSize)+radiusText+subtr+donutText+'.tif')
            # ascii_grid_writer(sumFile, cumCurrentArray, header, compress)
            
            # delete_data(prevSumFile)
            # prevSumFile = sumFile


#TEMP
        # bandCurrentAsciiFile = path.join(outputDir, 'BANDCur_addData2_band'+str(bandNum) +'.asc') 
        # header2=header
        # header2['nrows']=min(radius*2+1,radius+centerRow+1,header['nrows'])
        # ascii_grid_writer(bandCurrentAsciiFile , bandCurrentArray, header2, compress)
        # arcpy.ASCIIToRaster_conversion(bandCurrentAsciiFile, path.join(outputDir, 'addd_BAND'+str(bandNum)), 
                               # "FLOAT")
#TEMP

#put below in fn
        
        if solveInBand:
            print 'Done with band #',bandNum,'.'
            print 'Writing temporary grid...'
            bandRows=min(radius*2+1,radius+centerRow+1)
            yMin= max(header['yllcorner'],header['yllcorner'] + ((header['nrows'] - centerRow - radius - 1) * header['cellsize']))
            LLC = arcpy.Point(header['xllcorner'],yMin)

            bandCurrentRaster = arcpy.NumPyArrayToRaster(cumCurrentArray,LLC,
                                             header['cellsize'],header['cellsize'],-9999)
            # bandCurrentFile = os.path.join(outputDir,'bandcur_band'+str(bandNum)+'.tif')
            # bandCurrentRaster.save(bandCurrentFile)

            cumCurrentRaster = addData_arcpy(cumCurrentRaster, bandCurrentRaster)
            cumCurrentFile = os.path.join(outputDir,'BAND'+str(bandNum)+'startBand'+str(startBand)+resisRasterBase+'_sumtemp.tif')
            try:
                cumCurrentRaster.save(cumCurrentFile)
            except:
                try:
                    print 'Error writing'
                    time.sleep(5)                
                    scratchDir2=os.path.join(projectDir,scratchDirBase+'a')
                    if not path.exists(scratchDir2):
                        os.mkdir(scratchDir2)
                    arcpy.env.scratchWorkspace = scratchDir2
                    arcpy.env.Workspace = scratchDir2

                    cumCurrentFile = os.path.join(outputDir,'BAND'+str(bandNum)+'startBand'+str(startBand)+resisRasterBase+'_sumtemp2.tif')
                    cumCurrentRaster.save(cumCurrentFile)
                except:
                    print 'Second error writing'
                    time.sleep(5)
                    scratchDir3=os.path.join(projectDir,scratchDirBase+'b')
                    if not path.exists(scratchDir3):
                        os.mkdir(scratchDir3)
                    arcpy.env.scratchWorkspace = scratchDir3
                    arcpy.env.Workspace = scratchDir3
                    cumCurrentFile = os.path.join(outputDir,'BAND'+str(bandNum)+'startBand'+str(startBand)+resisRasterBase+'_sumtemp3.tif')
                    cumCurrentRaster.save(cumCurrentFile)            
            
      
            delete_data(prevCumCurrentFile)
            
            prevCumCurrentFile = cumCurrentFile

            if calcVoltages:
                bandVdiffRaster = arcpy.NumPyArrayToRaster(cumVdiffArray,LLC,
                                                 header['cellsize'],header['cellsize'],-9999)
                # bandCurrentFile = os.path.join(outputDir,'bandcur_band'+str(bandNum)+'.tif')
                # bandCurrentRaster.save(bandCurrentFile)

                cumVdiffRaster = addData_arcpy(cumVdiffRaster, bandVdiffRaster)
                cumVdiffFile = os.path.join(outputDir,'BAND'+str(bandNum)+'startBand'+str(startBand)+resisRasterBase+'_vdiffSumtemp.tif')
                try:
                    cumVdiffRaster.save(cumVdiffFile)
                except:    
                    try:
                        print 'Error writing'
                        time.sleep(5)
                        cumVdiffFile = os.path.join(outputDir,'BAND'+str(bandNum)+'startBand'+str(startBand)+resisRasterBase+'_vdiffSumtemp2.tif')
                        cumVdiffRaster.save(cumVdiffFile)
                    except:
                        print 'Second error writing'
                        time.sleep(5)
                        cumVdiffFile = os.path.join(outputDir,'BAND'+str(bandNum)+'startBand'+str(startBand)+resisRasterBase+'_vdiffSumtemp3.tif')
                        cumVdiffRaster.save(cumVdiffFile)
                    
                delete_data(prevVdiffFile)
                
                prevVdiffFile = cumVdiffFile
        
        print 'Done with band #',bandNum,'.Elapsed time so far: {0}'.format(datetime.datetime.now()-theStart)  

    print 'solved thru centerrow',centerRow
    fileBase, fileExtension = os.path.splitext(resisRasterBase)
    if len(fileBase)>10:
        fileBase=fileBase[0:14]+'Trnc'
    if startBand > 0:
        startBandText = '_startBand'+str(startBand)
    else:
        startBandText=''
    if stopBand > 0:
        stopBandText = '_stopBand'+str(stopBand)
    else:
        stopBandText=''
    if startStripe > 0:
        startStripeText = '_startStripe'+str(startStripe)
    else:
        startStripeText=''
    if stopStripe > 0:
        stopStripeText = '_stopStripe'+str(stopStripe)
    else:
        stopStripeText=''

    for fl in glob.glob(path.join(outputDir,'*.xml')):
        try:
            os.remove(fl)
        except:
            pass
    for fl in glob.glob(path.join(outputDir,'*.tfw')):
        try:
            os.remove(fl)
        except:
            pass
    
    sumFile = path.join(outputDir, 'cur_sum_' + nullText+fileBase + '__'+srcText+climText+logText+donutText+radiusText+'bl'+str(blockSize)+'_cutoff'+str(rCutoff) +weightText+noWeightText+subtrText+centerText+startBandText+stopBandText+startStripeText+stopStripeText+straddleText+maskDonutText+negTargText+fadeInText+fadeOutText+bufferText+'.tif')
    print 'writing:\n',sumFile
    cumCurrentRaster.save(sumFile)
    delete_data(prevCumCurrentFile)
    if calcVoltages:
        sumFile = path.join(outputDir, 'vDiff_sum_' + nullText+fileBase + '__'+srcText+climText+logText+donutText+radiusText+'bl'+str(blockSize)+'_cutoff'+str(rCutoff) +weightText+noWeightText+subtrText+centerText+startBandText+stopBandText+startStripeText+stopStripeText+straddleText+maskDonutText+negTargText+fadeVoltText+bufferText+'.tif')
        print 'writing:\n',sumFile
        cumVdiffRaster.save(sumFile)
        delete_data(prevVdiffFile)

    delete_dir(scratchDir)
    for filename in glob.glob(os.path.join(outputDir,'band*.*')) :
        try:    
            os.remove( filename )
        except:
            pass
        
def get_max_vdiff(voltMap,circleResisArray,options):
    voltMap_l, voltMap_r = get_horiz_neighbors(voltMap)
    voltMap_u, voltMap_d = get_vert_neighbors(voltMap)

    if adjustVoltages:
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
    if options['connect_four_neighbors_only'] == False:
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
    # print'vm'
    # print voltMap
    # print 'vmd'
    # print voltMap_direction
    vdiff_direction = npy.absolute(voltMap - voltMap_direction)
# fixme: nans, not -9999 are in these maps
    # vdiff_direction = npy.where(voltMap == -9999, 0, npy.where(voltMap_direction == -9999,0,vdiff_direction))
    vdiff_direction = npy.where(voltMap == npy.nan, 0, npy.where(voltMap_direction == npy.nan,0,vdiff_direction))

    # vdiff_direction = npy.where(vdiff_direction < 0, 0, vdiff_direction)
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


def band(inRaster,header,centerRow, radius):
#FIXME put bounds on bandrows- minrow=0, maxrow=nrows...
    # bandRows=min(radius*2+1,radius+centerRow+1,header['nrows'])
    bandRows=1 + min(radius,centerRow) + min(header['nrows'] - (centerRow+1), radius)
    yMin= max(header['yllcorner'],header['yllcorner'] + ((header['nrows'] - centerRow - radius - 1) * header['cellsize']))
    LLC = arcpy.Point(header['xllcorner'],yMin)
    # arcpy.Raster(inRaster).save('c:\\temp\\bandInRas')
    bandArray = arcpy.RasterToNumPyArray(inRaster,LLC,"#",bandRows,-9999)
    # newRaster = arcpy.NumPyArrayToRaster(bandArray,LLC,
                                         # header['cellsize'],header['cellsize'],-9999)
    # newRaster.save('c:\\temp\\bandarray')
    return bandArray

def addData(cumCurrentArray, currentArray, centerRow, centerCol, radius):
    # print 'adding data for center row, col, rad'
    # print centerRow,centerCol,radius
    minRow = max(0,centerRow-radius)
    maxRow = min(minRow+cumCurrentArray.shape[0]-1, centerRow+radius)
    minCol = max(0,centerCol-radius)
    # fixme: check if next line needs mincol added like in maxrow line
    maxCol = min(cumCurrentArray.shape[1]-1, centerCol+radius)
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

def donut(array, donutRadius, centerRow, centerCol):        
        grid2 = npy.indices((array.shape))
        rowArray2 = grid2[0]
        colArray2 = grid2[1]
        distArray = npy.sqrt(npy.multiply(rowArray2 - centerRow, rowArray2 - centerRow) + npy.multiply(colArray2-centerCol, colArray2-centerCol))           
        array = npy.where(distArray >= donutRadius, array, 0)
        return array

def circ(array, rowArray, colArray, centerRow, centerCol, radius):
        distArray = npy.sqrt(npy.multiply(rowArray - centerRow, rowArray- centerRow) + npy.multiply(colArray-centerCol, colArray-centerCol))           
        arrayMasked = npy.where(distArray <= radius, array, -9999)
        del distArray
        startRow = max(centerRow - radius,0)
        endRow = min(centerRow + radius,array.shape[0]-1)
        startCol = max(centerCol - radius,0)
        endCol = min(centerCol + radius,array.shape[1]-1)
        circleArray = arrayMasked[startRow:endRow+1,startCol:endCol+1]
        del arrayMasked
        return circleArray

def center_block(array, blockSize, centerRow, centerCol):
        startRow = centerRow - ((blockSize-1)/2)
        endRow = centerRow + ((blockSize-1)/2)
        startCol = centerCol - ((blockSize-1)/2)
        endCol = centerCol + ((blockSize-1)/2)
        blockArray = npy.zeros(array.shape, dtype = 'float64') - 9999    
        blockArray[startRow:endRow+1,startCol:endCol+1] = array[startRow:endRow+1,startCol:endCol+1]
        return blockArray            

def center_mask(array, blockSize, centerRow, centerCol):
        print 'masking',centerRow,centerCol
        startRow = centerRow - ((blockSize-1)/2)
        endRow = centerRow + ((blockSize-1)/2)
        startCol = centerCol - ((blockSize-1)/2)
        endCol = centerCol + ((blockSize-1)/2)
        maskedArray = array
        maskedArray[startRow:endRow+1,startCol:endCol+1] = 0
        return maskedArray            

def get_center_block(array, blockSize, centerRow, centerCol):
        print 'getting center block',centerRow,centerCol
        startRow = centerRow - ((blockSize-1)/2)
        endRow = centerRow + ((blockSize-1)/2)
        startCol = centerCol - ((blockSize-1)/2)
        endCol = centerCol + ((blockSize-1)/2)
        blockArray = npy.where(array > 0, -9999, -9999)
        blockArray[startRow:endRow+1,startCol:endCol+1] = 1
        return blockArray            

def get_subset_header(array, fullHeader, radius, centerRow, centerCol):
    subsetHeader = {}
    llrow = min(centerRow + radius, fullHeader['nrows']-1)
    llcol = max(centerCol - radius, 0)
   
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
    # f = open(filename, 'r')
    # ncols = int(string.split(f.readline())[1])
    # nrows = int(string.split(f.readline())[1])
    # xllcorner = float(string.split(f.readline())[1])
    # yllcorner = float(string.split(f.readline())[1])
    # cellsize = float(string.split(f.readline())[1])
    # try:
        # [_ign, nodata] = string.split(f.readline())
        # try:
            # nodata = int(nodata)
        # except ValueError:
            # nodata = float(nodata)
    # except ValueError:
        # nodata = False
    # header = {}
    # header['ncols'] = ncols
    # header['nrows'] = nrows
    # header['xllcorner'] = xllcorner
    # header['yllcorner'] = yllcorner
    # header['cellsize'] = cellsize
    # header['nodata'] = -9999 if (nodata == False) else nodata     
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

def ascii_grid_writer(file_name, data, header, compress):
            """Writes rasters to ASCII grid or numpy formats."""     

            f = gzip.open(file_name+'.gz', 'w') if compress else open(file_name, 'w')    
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
    # try:
    
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
        
    # # Return GEOPROCESSING specific errors
    # except arcpy.ExecuteError:
        # lu.dashline(1)
        # print('****Failed in step 8. Details follow.****')
        # lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # # Return any PYTHON or system specific errors
    # except:
        # lu.dashline(1)
        # print('****Failed in step 8. Details follow.****')
        # lu.exit_with_python_error(_SCRIPT_NAME)


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
    options = {}

    options['data_type']='raster'
    options['version']='unknown'
    options['low_memory_mode']=False
    options['scenario']='advanced'
    options['habitat_file']='(Browse for a habitat map file)'
    options['habitat_map_is_resistances']=True
    options['point_file']=('(Browse for file with '
                          'locations of focal points or areas)')
    options['point_file_contains_polygons']=True
    options['connect_four_neighbors_only']=False
    options['connect_using_avg_resistances']=True
    options['use_polygons']=False
    options['polygon_file']='(Browse for a short-circuit region file)'
    options['source_file']='(Browse for a current source file)'
    options['ground_file']='(Browse for a ground point file)'
    options['ground_file_is_resistances']=True
    options['use_unit_currents']=False
    options['use_direct_grounds']=False
    options['remove_src_or_gnd']='rmvsrc'
    options['output_file']='(Choose a base name for output files)'
    options['write_cur_maps']=True
    options['write_cum_cur_map_only']=True
    options['log_transform_maps']=False
    options['write_volt_maps']=False
    options['solver']='cg+amg'
    options['compress_grids']=False
    options['print_timings']=False
    options['use_mask']=False
    options['mask_file']='None'
    options['use_included_pairs']=False
    options['included_pairs_file']='None'
    options['use_variable_source_strengths']=False
    options['variable_source_file']='None'
    options['write_max_cur_maps']=False
    options['set_focal_node_currents_to_zero']=True
    options['print_timings']=False
    options['output_file'] = outputFile
       
    options['remove_src_or_gnd']='rmvsrc'#'keepall'
    options['habitat_file'] = resisAsciiFile
    options['source_file'] = sourceAsciiFile
    options['ground_file']= groundAsciiFile

    return options

def writeCircuitscapeConfigFile(configFile, options):
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

    if options['ground_file_is_resistances']=='not entered':
        options['ground_file_is_resistances'] = False
    if options['point_file_contains_polygons']=='not entered':
        options['point_file_contains_polygons'] = False

    for option in sections:
        try:
            config.add_section(sections[option])
        except:
            pass
    for option in sections:
        config.set(sections[option], option, options[option])

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
    circleFlowBand()

# MULTISCALE
# assume 100m pix
# 1-3km (11/10/30)
# 3-10 (11/30/100)
# 10-20(11/100/200)
# 20-50 (11/200/500)
# projectDir = 'c:\\temp\\CF_stgr'
# outputDir = 'c:\\temp\\CF_stgr\csblockout4sub'
# scratchDir = 'c:\\temp\\CF_stgr\\scratch4'

# Tiling:
    # take large area, create overlapping tiles based on max radius. Run code
    # subset a strip each time. Height = radius*2+1. Start with first full strip, increment by blockSize

# beamtiles
# for tile, take 3x3 subset
# beam up, down, etc

# OR wedges.... could test by running 4x, with code that snuffs out sources outsid of each wedge
    # Doesn't seem to help, outer circles show up more


# adjacent tiles- beam up, beam down (and ne-sw etc, / sqrt 2.)
    # -could shift
    # -fewer calcs, even with some shifts, because tiles would be on order of 50km, not 5 pix
    # don't have block problems

# # VOLTAGE###
# dist matrix, get everything within radius
# fill radius with max diff
# BUT WOULD NEED TO DO MOVING WINDOW ACROSS ENTIRE MATRIX LIKE FOC STATS
# http://stackoverflow.com/questions/10996769/pixel-neighbors-in-2d-array-image-using-python
# https://geonet.esri.com/thread/16548
# ndimage.filters:
# http://gis.stackexchange.com/questions/34306/focal-statistics-with-variable-radius
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
        # sourceRaster = arcpy.sa.Con((arcpy.Raster(resisRaster) > rCutoff),0,1)

        
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
    # fadeindist of radius/2 also gets rid of blocking. But nice symmetry to radius- add 2 wedges, get a line.
    
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
# if useClimate
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
# Try sourceOnly- scale currents by # sources
# try scaling all equally- I guess the same as setting all points to be sources?
    #is this mathematically equiv to targ only?
    
# Go back to donut? model LDD (e.g. 10-50km)?
# Try a jelly donut. 50km circle, only recording middle 20km.
# To match W2W, needs to be radius of 1500 cells (!) with collection area of 750 cells.

# what about triangle fadeinfadeout where center is most important
  
# 12/15
# Put aside concerns about speed. can do parallel runs, coarsening and such. 

# proof of concept
    # ecology grid- does approach approximate pairwise?
    # yes. some open questions about block size effects.
    
# Null model to distinguish areas with high flow because of many src/target points vs flow pinchpoints
    
# try with seaboard
    # what happens if center ground is in a city?

# if not connecting nearby points, could match pairwise
# try two-scale? just outer half of large radius, then inner half with smaller block size?
    # would need to put sources at outer edge of inner half to match current contributed by outer area!
    
# compare outside only to exclude pairwise

# fadeinStart
# fadeInEnd

# fade in doing trick. test results against much finer run, like 15
# SPEED
# -test with large blocks 
# -PARALLEL- 4-8 at a time on one machine
# CLIMATE
# -Would need small blocks so have homog climate!
    # bl105=105*270m ~ 28km!
    # ..OR.. pixel-by-pixel strength pairing but still consolidated into blocks
        # would need to test how well it does



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

# # resisRasterBase = 'resis100km.asc'
# # resisRasterBase = 'ecolresis_100m.asc'

# # resisRasterBase = 'ecol_resis3_nodata2_10kbarriers2.asc'
# # resisRasterBase = 'EcologyResistances100k.asc'
# # sourceRasterBase = 'EcologySourcePoints.asc'
# # resisRasterBase = 'ecolres_halo.asc'
# # resisRasterBase = 'nlcdres_ne_NODATA.ASC'
# # resisRasterBase = 'resistance_surface_se.asc'
# # resisRasterBase = 'resistance_surface_se_nodata_270ma.asc'
# # resisRasterBase = 'tip_peninsula_resis_270m.ASC'
# # resisRasterBase = 'klamath2_r_md2.asc'
# # resisRasterBase = 'EcologyResistances.asc'
# resisRasterBase = '5x5R.asc'
# # resisRasterBase = 'ecol_resis3_nodata2_10kbarriers2.asc'
# resisRasterBase = 'ecol_resis3_nodata2.asc'

# useSourceRaster=False
# sourceRasterBase = '5x5_src.asc'
# sourceRasterBase = 'proof_points39.asc'

# useClimate=False
# tCutoff=3
# climateRasterBase ='5x5Temperature.asc'
# climateRasterBase = 'ecologytemperature.asc'

# projectDir = 'C:\\Dropbox\\Working\\Circuitscape_BraidedThruway\\CF_PROOF'
# outputDirBase = 'CF_5x5'
# outputDirBase = 'CF_climEcol'

# scratchDirBase = 'scratch'

# projectDir = 'C:\\dropbox\\Working\\Circuitscape_BraidedThruway\\CF_NE_NLCD2'
# outputDirBase = 'NY'
# # resisRasterBase = 'ny_nlcd_500m.asc'
# resisRasterBase = 'ny_nlcd_166m_water5_nodata.asc'
# # resisRasterBase = 'ny_nlcd_166m_water5_nodata_adks2.asc'



# projectDir = 'C:\\dropbox\\Working\\Circuitscape_BraidedThruway\\CF_NE_CLIM_Coarse'
# outputDirBase = 'adks450m'
# # resisRasterBase = 'ny_nlcd_500m.asc'
# resisRasterBase = 'resis_sq_adks_450m.ASC'
# resisRasterBase = 'resis_sq_adks_450m.tif'

# projectDir = 'C:\\temp'

# outputDirBase = 'neTest'
# # resisRasterBase = 'ny_nlcd_500m.asc'
# resisRasterBase = 'resistance.tif'
# useClimate=True
# tCutoff=0.1
# climateRasterBase ='tmean.tif'


# projectDir = 'D:\\GIS_DATA\\NACR\\McRae\\circuitscapeBraidedThruway\\CF_NE'
# outputDirBase = 'tile_test7'
# # resisRasterBase = 'ny_nlcd_500m.asc'
# resisRasterBase = 'resis_chesapeak2.tif'
# useClimate=True
# tCutoff=0.5
# climateRasterBase ='tmean_ne_chesapeak.tif'
# tile=1
# xxxx


# projectDir = 'C:\\dropbox\\Working\\Circuitscape_BraidedThruway\\CF_PROOF'
# # resisRasterBase = 'EcologyResistances.asc'
# resisRasterBase = '6x9R.asc'
# rCutoff=1
# useClimate=False
# resisRasterBase = 'ecol_resis3_nodata2_10kbarriers2.asc'
# resisRasterBase = 'ecol_resis3_nodata2.asc'


# projectDir = 'C:\\Dropbox\\Working\\Circuitscape_BraidedThruway\\CF_stgr'
# outputDir = 'C:\\Dropbox\\Working\\Circuitscape_BraidedThruway\\CF_stgr\\cf_stgr_srcs'
# scratchDir = 'C:\\Dropbox\\Working\\Circuitscape_BraidedThruway\\CF_stgr\\scratch'

# resisRasterBase = 'typh_r.ASC'
# sourceRasterBase = 'typh_sources.ASC'
# projectDir = 'c:\\temp\\CF_BAND'
# outputDir = 'c:\\temp\\CF_BAND\\csbandecolrVr100k_VDacross'
# scratchDir = 'c:\\temp\\CF_BAND\\scratchbandecolr'


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

    # circle resis array (SUBSET of radius*2 on a side, ND outside circle)
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
    # after this line: tCutoffArray=npy.where(tDiffArray>=tCutoff,1,0)
