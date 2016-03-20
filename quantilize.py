raster = r'C:\DATADRIVE\DUKE_PNW_DATA\PNW_Omniscape_180m\acur_R_d8_clpF_180m_src_Hab_d8_clp_180m__r278b31_targOnlycg_clp_artifactCorr_540m.tif'#'C:\DATADRIVE\DUKE_PNW_DATA\PNW_Omniscape_540m_Mockup\CS_540m_100km_Lim5\CS_540m_100km_Lim5.gdb\current_MOS_cur_draft6R540m0_lim925_srcRas_draft6Hab540m0__r185b25_distFn_targOnlycg_f25clipC'
extractMask =r'C:\DATADRIVE\DUKE_PNW_DATA\CIRCUITSCAPE_RESISTANCES_DUKE_PNW_CA\DUKE_FINAL_RESIS_RASTERS\PNW_study_area_poly_BHM.shp' #None to ignore, otherwise will cip to this

numQuantiles=100

import os 
import arcpy
import math
import ConfigParser
import string
import os.path as path

import sys
import numpy as npy
import shutil
import glob

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True
 
 
def quantilize(raster,numQuantiles):
    # try:
    
        if extractMask is not None:
            dir,file=os.path.split(raster)
            fileBase,ext = os.path.splitext(file)
                       
            print "Extracting by mask file"
            extractRasFN = fileBase+'clp.tif'
            extractRasPath = os.path.join(dir,extractRasFN)
            outRas = arcpy.sa.ExtractByMask(raster, extractMask)
            outRas.save(extractRasPath)
            # delete_data(os.path.join(outputDir,mosFN+'.tif'))
            raster = extractRasPath
   
    
    
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
            
            dir,file=os.path.split(raster)
            outPath=os.path.join(dir,'p'+file)
            
            quantileRaster.save(outPath)
            print 'Saved to:\n',outPath
            # except:
        # print'failed to quantilize'
       
def seq(start, stop, step=1):
    n = int(round((stop - start)/float(step)))
    if n > 1:
        return([start + step*i for i in range(n+1)])
        
if __name__ == '__main__':
    quantilize(raster,numQuantiles)
    
print 'Done.'    
        