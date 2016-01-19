projectDir = r'D:\GIS_DATA\NACR\McRae\Duke_PNW_Omniscape'# this is where all the input data are, and where out directory will be created.
inputDirBase = 'd6_540m_100km_NoLim'#'d6_540m_100km_5Lim'
outputDirBase = inputDirBase + 'AddedRasters'
numQuantiles=100
extractMask =r'D:\gis_data\NACR\McRae\Duke_PNW_Omniscape\PNW_study_area_poly_BHM.shp'
tasks = ['flow']#['cur']#,'flow']#,'volt']
    
import os
import shutil
import re
import arcpy
from arcpy.sa import *
arcpy.CheckOutExtension("spatial")
arcpy.env.overwriteOutput = True  

def mosaic():
    inputDir = os.path.join(projectDir,inputDirBase)
    outputDir = os.path.join(inputDir,outputDirBase)
    delete_dir(outputDir)
    if not os.path.exists(outputDir):
        os.mkdir(outputDir)

    
    for task in tasks:
        print 'Task=',task
        tileDir = os.path.join(projectDir,inputDirBase)
        if not os.path.exists(tileDir):
            print'Path does not exist:',tileDir
            continue

        for f in os.listdir(tileDir):
            fn,ext = os.path.splitext(f)

            if task == 'flow':
                if ext == '.tif' and (fn[0]=='f'):
                    flowRas = os.path.join(tileDir,f)
                    zeroFlowRas = arcpy.sa.Con(arcpy.sa.IsNull(flowRas),0,flowRas)
                    dir,file=os.path.split(flowRas)
                    zFRasPath = os.path.join(dir,'z'+file)
                    zeroFlowRas.save(zFRasPath)

        expression = ''
        
        count=0
        for f in os.listdir(tileDir):
            fn,ext = os.path.splitext(f)
                
            if ext == '.tif' and '_SB' in fn and ((fn[0]=='v' and task == 'volt') or (fn[0]=='c' and task == 'cur') or (fn[0]=='z' and task == 'flow')):
                rasterPath= os.path.join(tileDir,f)
                count=count+1
                if count==1:
                    print 'rasterPath is ',rasterPath
                    fn=re.sub(r'_SB\w+_f', '_f', fn)
                    mosFN = 'a'+fn.replace('.','')
                    outRaster = outputDir+'\\'+mosFN+'.tif'
                    delete_data(outRaster)
                    clipRaster = outputDir+'\\f'+mosFN+'clp.tif'         
                    delete_data(clipRaster)
                    expression = rasterPath
                else:
                    # rasterString = rasterString + '"'+rasterPath+'"'+";"  
                   
                    expression = expression + " + " +rasterPath 
        # expression = (expression + " + " + cwdRaster2 + " - " 
                      # + lcDist)

        print "Adding:"
        print expression
        print 'to:'
        print outRaster
        # outRas = Raster("input1") + Raster("input2") 

        arcpy.gp.SingleOutputMapAlgebra_sa(expression, outRaster)
        # statement = ('gp.SingleOutputMapAlgebra_sa(expression, '
                     # 'outRaster)')
        count = 0
        
        if extractMask is not None:
            print "Extracting by mask file"
            extractRasFN = mosFN+'clp.tif'
            extractRasPath =os.path.join(outputDir,extractRasFN)
            outRas = arcpy.sa.ExtractByMask(os.path.join(outputDir,mosFN+'.tif'), extractMask)
            outRas.save(extractRasPath)
            # delete_data(os.path.join(outputDir,mosFN+'.tif'))
            raster = extractRasPath
        else:
            raster = os.path.join(outputDir, mosFN + '.tif')
        if numQuantiles is not None and task == 'cur':
            quantilize(raster, numQuantiles)
        # arcpy.env.Workspace=inputDir
        # outInt = Int(Raster(os.path.join(outGDB,mosFN)))

        # rasterPath = os.path.join(outGDB,'INT_'+mosFN)
        # outInt.save(rasterPath)

        
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
            dir,file=os.path.split(raster)
            outPath=os.path.join(dir,'p'+file)
            
            quantileRaster.save(outPath)
            print 'Saved quantilized raster to:',outPath
    # except:
        # print'failed to quantilize'
       
def seq(start, stop, step=1):
    n = int(round((stop - start)/float(step)))
    if n > 1:
        return([start + step*i for i in range(n+1)])

        
def delete_dir(dir):
    try:
        if os.path.exists(dir):
            shutil.rmtree(dir)
        return
    except:
        return
        
if __name__ == '__main__':
    mosaic()
    
print 'Done.'    