new mosaic code
take all 'cur' files, add them
clip



Need to quantify regional flow. 
    look at 10km 20km focal stats.

COUNT of above average regional, local, etc
regional BOOST
multiply regional and local
degree of constriction
regional flow
Variance-

run 540m null and reg on cospatial. use cwd cutoff. RUNNING COSPAT
run 270m draft7 on 2 machines- MA and MN. DON'T clip to 100km! ideally would have 200km buff. 750 pixel overlap. NO random. VOLTAGES.

run 1080 flow on laptop. 

null current closely approximates mean hab vaLUE in 50km radius. 
suggests regional importance can be quickly estimated by summing current sources in analysis window
BUT, may be best to just run null. need to reclassify nodata resistances

# RUN FLOW AT 1080m
# RUN CS ACROSS PLATFORMS WITH DRAFT7
    # need to run with bigger overlap, also with no random

# create clip raster from nationalmap.gov states data and ecoregions. buffer out 270m or so. OUTSIDE IS ZERO NOT NODATA
ecoregions+1kmbuff + Duke_States_250mBuff. study area is already buffered out up to 10km from boundaries of non-CA ecoregions. CA ecoregions are buffered out about 1km.
re-buffer ecoregions by 1km all around?
intersected Duke_States_250mBuff with ecoregions buffered by 1km, including ID extension and northern WA. Got C:\DATADRIVE\DUKE_PNW_DATA\CIRCUITSCAPE_RESISTANCES_DUKE_PNW_CA\DUKE_FINAL_RESIS_RASTERS\PNW_study_area_poly_BHM.shp
 
# low-avg-high regional flow (5km)
# low-avg-high within each regional bracket, or at least top two
# LOW-AVG-HIGH WITHIN EACH ECOREGION

# count of high within ecoregion plus high within regional bracket, plus high within 5km?


startTile=0
endTile=8
projectDir = r'C:\DATADRIVE\DUKE_PNW_DATA\PNW_Omniscape_540m_Mockup'
outputDirBase = 'CS_540m_100km_Lim5_NULL'
extractMask =r'C:\DATADRIVE\DUKE_PNW_DATA\PNW_2014_DATA_201502\PNW_RESILIENCE_20150201.gdb\ECOREGIONS_ALL_ECOREG'
tasks = ['current']#,'flow']#,'voltage']
mosaicMethod="BLEND" # Either mean with 1-pixel overlap or blend with 5 seem to work. mean might be better
numTrimPixels=180
numQuantiles=20

if numTrimPixels > 0:
    print 'Rasters will be trimmed by',numTrimPixels,'pixels before mosaicking.'
    
import os
import arcpy
from arcpy.sa import *
arcpy.CheckOutExtension("spatial")
arcpy.env.overwriteOutput = True  

def mosaic():
    mosaicDir = os.path.join(projectDir,outputDirBase)
    if not os.path.exists(mosaicDir):
        os.mkdir(mosaicDir)

    
    for task in tasks:
        print 'Task=',task
        rasterString = ''
        for tile in range(startTile,endTile+1):
            print 'Tile:',tile
            
            tileDir = os.path.join(projectDir,outputDirBase+str(tile))
            if not os.path.exists(tileDir):
                print'Path does not exist:',tileDir
                continue
            
            for f in os.listdir(tileDir):
                fn,ext = os.path.splitext(f)


                if ext == '.tif' and ((fn[0]=='v' and task == 'voltage') or (fn[0]=='c' and task == 'current') or (fn[0]=='f' and task == 'flow')):
                    if tile==startTile:
                        mosFN = task + '_MOSB_' + (fn.replace('clip_'+str(tile),'')).replace('.','')
                    rasterPath= os.path.join(tileDir,f)
                    # rasterString = rasterString + '"'+rasterPath+'"'+";"  

                    if numTrimPixels > 0:
                        arcpy.env.snapraster=rasterPath
                        descData=arcpy.Describe(rasterPath)
                    
                        extent=descData.Extent
                        cellSize=descData.meanCellHeight
                        trimDistance = numTrimPixels * cellSize
                        print 'Trimming input raster by', trimDistance
                        # spatialReference=descData.spatialReference
                        arcpy.env.extent = extent
                        
                        arcpy.env.extent = arcpy.Extent(extent.XMin+trimDistance,extent.YMin+trimDistance,extent.XMax-trimDistance,extent.YMax-trimDistance)
                        outFile = 't'+f
                        trimRasterPath = os.path.join(tileDir,outFile)
                        conRas = arcpy.sa.Con(Raster(rasterPath) > -1, Raster(rasterPath))
                        conRas.save(trimRasterPath)
                        arcpy.env.extent = "MAXOF"     
                        rasterPath = trimRasterPath
                    
                    rasterString = rasterString + '"'+rasterPath+'"'+";"  
        
        outGDB= os.path.join(mosaicDir,outputDirBase+'.gdb')
        if not os.path.exists(outGDB):
            arcpy.CreateFileGDB_management(mosaicDir, outputDirBase+'.gdb')
        print "Mosaicking:"
        print rasterString
        print 'to:'
        print outGDB+'\\'+mosFN
        arcpy.MosaicToNewRaster_management(rasterString,outGDB,mosFN, "", "32_BIT_FLOAT", "#", "1", "BLEND", "MATCH")                        
        
        if extractMask is not None:
            print "Extracting by mask file"
            extractRasFN = mosFN+'clipB'
            extractRasPath =os.path.join(outGDB,extractRasFN)
            outRas = arcpy.sa.ExtractByMask(os.path.join(outGDB,mosFN), extractMask)
            outRas.save(extractRasPath)
            raster = extractRasPath
        else:
            raster = os.path.join(outGDB, mosFN)
        if numQuantiles is not None and task == 'current':
            quantilize(raster, numQuantiles)
        # arcpy.env.Workspace=mosaicDir
        # outInt = Int(Raster(os.path.join(outGDB,mosFN)))

        # rasterPath = os.path.join(outGDB,'INT_'+mosFN)
        # outInt.save(rasterPath)
        
def quantilize(raster,numQuantiles):
    # try:
        
        import numpy as npy
        if raster is not None:       
            inArray = arcpy.RasterToNumPyArray(raster,"#","#","#",-9999)
            interval = 100.0/numQuantiles
            breakSequence = seq(interval, 100-interval, interval)
            quantileArray = npy.zeros_like(inArray, dtype='int32')        
            quantileArray = npy.where(inArray > 0,numQuantiles,quantileArray)
                       
            inVector = inArray.flatten()
            ind=npy.where(inVector > 0) 
            inVector = inVector[ind] #Eliminate nodata and zeros

            # quantilize                      
            if len(inVector)==0:
                print 'Array is all zero. No results to process.'
                return
            print '\nPartitioning ' + str(len(inVector)) + ' non-zero values into ' + str(numQuantiles) + ' quantiles.'
            quantileBreaks = npy.percentile(inVector, breakSequence)

            # fixme: consider replacing below with an arcpy reclassify command. Would be faster and take less memory.
            # Not critical and remap table is fussy- can't be a string
            for i in range(numQuantiles-2,-1,-1):
                quantileArray = npy.where (inArray<quantileBreaks[i],i+1,quantileArray)
            
            
            
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
            # fileBase,ext=os.path.splitext(raster)
            outPath=raster+'_PCT'
            
            quantileRaster.save(outPath)
            print 'Saved quantilized raster to:',outPath
    # except:
        # print'failed to quantilize'
       
def seq(start, stop, step=1):
    n = int(round((stop - start)/float(step)))
    if n > 1:
        return([start + step*i for i in range(n+1)])
        
if __name__ == '__main__':
    mosaic()
    
print 'Done.'    