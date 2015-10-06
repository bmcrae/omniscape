startTile=0
endTile=8
projectDir = r'D:\GIS_DATA\NACR\McRae\circuitscapeBraidedThruway\CF_NE'
outputDirBase = 'cf_ne_tile025deg'

import os
import arcpy
from arcpy.sa import *
arcpy.CheckOutExtension("spatial")
arcpy.env.overwriteOutput = True  

mosaicDir = os.path.join(projectDir,outputDirBase)
if not os.path.exists(mosaicDir):
    os.mkdir(mosaicDir)

tasks = ['current','voltage']
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
            if ext == '.tif' and ((fn[0]=='v' and task == 'voltage') or (fn[0]=='c' and task == 'current')):
                if tile==startTile:
                    mosFN = task + '_MOSAIC_' + (fn.replace('clip_'+str(tile),'')).replace('.','')
                rasterPath= os.path.join(tileDir,f)
                rasterString = rasterString + '"'+rasterPath+'"'+";"  

    outGDB= os.path.join(mosaicDir,outputDirBase+'.gdb')
    if not os.path.exists(outGDB):
        arcpy.CreateFileGDB_management(mosaicDir, outputDirBase+'.gdb')
    print "Mosaicking:"
    print rasterString
    print 'to:'
    print outGDB
    print mosFN
    arcpy.MosaicToNewRaster_management(rasterString,outGDB,mosFN, "", "32_BIT_FLOAT", "#", "1", "MAXIMUM", "MATCH")                        
    
    # arcpy.env.Workspace=mosaicDir
    outInt = Int(Raster(os.path.join(outGDB,mosFN)))

    rasterPath = os.path.join(outGDB,'INT_'+mosFN)
    outInt.save(rasterPath)
    
