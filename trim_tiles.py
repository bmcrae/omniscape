# trim rasters- use un-extended tiles to do this!
# create clip raster from nationalmap.gov states data and ecoregions. buffer out 270m or so. OUTSIDE IS ZERO NOT NODATA
# low-avg-high regional flow (5km)
# low-avg-high within each regional bracket, or at least top two
# LOW-AVG-HIGH WITHIN EACH ECOREGION

# count of high within ecoregion plus high within regional bracket, plus high within 5km?

startTile=0
endTile=8
projectDir = r'C:\DATADRIVE\DUKE_PNW_DATA\PNW_Omniscape_540m_Mockup'
outputDirBase = 'CS_540m_100km_Lim5_NULL'
extractMask =r'C:\DATADRIVE\DUKE_PNW_DATA\PNW_2014_DATA_201502\PNW_RESILIENCE_20150201.gdb\ECOREGIONS_ALL_ECOREG'
numQuantiles = 100
tasks = ['current']#,'flow']#,'voltage']
mosaic_trimmed_rasters = True
numTrimPixels = 370

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
            # print tileDir
            # print os.listdir(tileDir)
            # files = [f for f in os.listdir(tileDir) if os.path.isfile(f)]
            # print 'FILES'
            # print files
            
            for f in os.listdir(tileDir):
                fn,ext = os.path.splitext(f)
                if ext == '.tif' and ((fn[0]=='v' and task == 'voltage') or (fn[0]=='c' and task == 'current') or (fn[0]=='f' and task == 'flow')):
                    if tile==startTile:
                        mosFN = task + '_MOSAIC_' + (fn.replace('clip_'+str(tile),'')).replace('.','')
                    rasterPath= os.path.join(tileDir,f)
                    if numTrimPixels > 0:
                        arcpy.env.snapraster=rasterPath
                        descData=arcpy.Describe(rasterPath)
                    
                        extent=descData.Extent
                        cellSize=descData.meanCellHeight
                        if trimtrimDistance = numTrimPixels * cellSize
                        print 'Trimming input raster by', trimDistance
                        # spatialReference=descData.spatialReference
                        arcpy.env.extent = extent
                        
                        arcpy.env.extent = arcpy.Extent(extent.XMin+trimDistance,extent.YMin+trimDistance,extent.XMax-trimDistance,extent.YMax-trimDistance)
                        outFile = 't'+f
                        outRasterPath= os.path.join(tileDir,outFile)
                        conRas = arcpy.sa.Con(Raster(rasterPath) > -1, Raster(rasterPath))
                        conRas.save(outRasterPath)
                        arcpy.env.extent = "MAXOF"
                    
                    # rasterString = rasterString + '"'+rasterPath+'"'+";"  
                    
        # outGDB= os.path.join(mosaicDir,outputDirBase+'.gdb')
        # if not os.path.exists(outGDB):
            # arcpy.CreateFileGDB_management(mosaicDir, outputDirBase+'.gdb')
        # print "Mosaicking:"
        # print rasterString
        # print 'to:'
        # print outGDB
        # print mosFN
        # arcpy.MosaicToNewRaster_management(rasterString,outGDB,mosFN, "", "32_BIT_FLOAT", "#", "1", "BLEND", "MATCH")                        
        
        # if extractMask is not None:
            # print "Extracting by mask file"
            # extractRasFN = mosFN+'clipB'
            # extractRasPath =os.path.join(outGDB,extractRasFN)
            # outRas = arcpy.sa.ExtractByMask(os.path.join(outGDB,mosFN), extractMask)
            # outRas.save(extractRasPath)
            # raster = extractRasPath
        # else:
            # raster = os.path.join(outGDB, mosFN)
        # if numQuantiles is not None and task == 'current':
            # quantilize(raster, numQuantiles)
        # # arcpy.env.Workspace=mosaicDir
        # # outInt = Int(Raster(os.path.join(outGDB,mosFN)))

        # # rasterPath = os.path.join(outGDB,'INT_'+mosFN)
        # # outInt.save(rasterPath)
        
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