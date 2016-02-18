projectDir = r'D:\Users\bmcrae\Duke_PNW_Omniscape\d8_omniscape' #r'C:\DATADRIVE\DUKE_PNW_DATA\PNW_OmniScape'# this is where all the input data are, and where out directory will be created.
resisRasterBase = 'R_d8_clpF_180m.tif'
blockSize=1

import os 
import arcpy
import math
import os.path as path
import numpy as npy
import shutil
arcpy.env.overwriteOutput = True

def get_header(filename):
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

raster=os.path.join(projectDir,resisRasterBase)
inArray = arcpy.RasterToNumPyArray(raster,"#","#","#",-9999)
bandArray=npy.zeros(inArray.shape, dtype = 'int32') 
stripeArray=bandArray.copy()
header = get_header(raster)
outDir,rasterName = os.path.split(raster)

bandRasterPath = os.path.join(outDir,resisRasterBase+'_bandsBlockSize'+str(blockSize)+'.tif')
stripeRasterPath = os.path.join(outDir,resisRasterBase+'_stripesBlockSize'+str(blockSize)+'.tif')
bandNum=0
stripeNum=0
print 'Processing',raster
for centerRow in range((blockSize-1)/2,header['nrows'],blockSize):
    bandNum+=1
    bandArray[centerRow-(blockSize-1)/2:centerRow+(blockSize-1)/2+1,:]=bandNum 
LLC = arcpy.Point(header['xllcorner'],header['yllcorner'])
bandRaster = arcpy.NumPyArrayToRaster(bandArray,LLC, header['cellsize'],header['cellsize'],-9999)    
bandRaster.save(bandRasterPath)
print 'Saved',bandRasterPath        

for centerCol in range((blockSize-1)/2, header['ncols'],blockSize):
    stripeNum+=1
    stripeArray[:,centerCol-(blockSize-1)/2:centerCol+(blockSize-1)/2+1]=stripeNum 
stripeRaster = arcpy.NumPyArrayToRaster(stripeArray,LLC, header['cellsize'],header['cellsize'],-9999)    
stripeRaster.save(stripeRasterPath)
print 'Saved',stripeRasterPath
