projectDir = r'D:\GIS_DATA\NACR\McRae\CA_Mockups\HM_2016_0811'# this is where all the input data are, and where out directory will be created.

inputDirBase = 'rHMp10_s1minusHM_90m'#'d6_540m_100km_5Lim'
outputDirBase = inputDirBase + 'AddedRasters'
numQuantiles=100
extractMask = None #r'C:\DATADRIVE\DUKE_PNW_DATA\CIRCUITSCAPE_RESISTANCES_DUKE_PNW_CA\DUKE_FINAL_RESIS_RASTERS\PNW_study_area_poly_BHM.shp' #None to ignore, otherwise will cip to this
tasks = ['cur']#['cur','volt']
    
import os
import glob
import shutil
import re
import arcpy
from arcpy.sa import *
arcpy.CheckOutExtension("spatial")
arcpy.env.overwriteOutput = True  

def mosaic():
    inputDir = os.path.join(projectDir,inputDirBase)
    if not os.path.exists(inputDir):
        print("Error. Directory " + inputDir + " doesn't exist.")
        raw_input('Press Enter to exit.')
        exit(0)
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
                
            if ext == '.tif' and '_SB' in fn and '_PCT' not in fn and ((fn[0]=='v' and task == 'volt') or (fn[0]=='c' and task == 'cur') or (fn[0]=='z' and task == 'flow') or (fn[0]=='s' and task == 'sources') or (fn[0]=='t' and task == 'targets')):
                rasterPath= os.path.join(tileDir,f)
                count=count+1
                if count==1:
                    print 'rasterPath is ',rasterPath
                    fn=re.sub(r'_SB\w+.', '.', fn)
                    mosFN = 'a'+fn.replace('.','')
                    outRaster = outputDir+'\\'+mosFN+'.tif'
                    delete_data(outRaster)
                    clipRaster = outputDir+'\\f'+mosFN+'_clp.tif'         
                    delete_data(clipRaster)
                    expression = rasterPath
                else:
                    # rasterString = rasterString + '"'+rasterPath+'"'+";"  
                   
                    expression = expression + " + " +rasterPath 
        # expression = (expression + " + " + cwdRaster2 + " - " 
                      # + lcDist)

                      
        print 'Copying Python files from input to output directory:'
        for file in glob.glob(os.path.join(inputDir,'*.py')):
            print file                                                                                                                                        
            shutil.copy(file, outputDir)        
        
        
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
            extractRasFN = mosFN+'_clp.tif'
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
            print 'Saved quantilized raster to:'
            print outPath
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