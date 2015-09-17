# ONLY CONNECT AREAS that are in different patches? so use a different rsrc layer for each target.


# #BHM NOTES
# Iterates through FIDs. If 10 iterations, will do 1st 10 points starting at whatever FID specified in tool
# need to add up results
# can limit CWD, could also mask out area beyond X km away from seed point.

# ---------------------------------------------------------------------------
# LandscapePermeability.py
# Created on: 11 October 2011
# Written by David Theobald, Colorado State University
# ---------------------------------------------------------------------------
import arcpy, os, sys
from arcpy import env
from arcpy.sa import *
arcpy.env.overwriteOutput=True

def delete_data(dataset):
    try:
        if arcpy.Exists(dataset):
            arcpy.Delete_management(dataset)
    except:
        pass
    return


inFC = sys.argv[1]      # input random points, ordered randomly!!!
rCW = sys.argv[2]
rNL0 = sys.argv[3]

arcpy.AddMessage("Define arguments complete")

maxCWD = int(sys.argv[4])
numIterations = int(sys.argv[5])
startIteration = int(sys.argv[6])
rOut = sys.argv[7]

ws = env.workspace

arcpy.AddMessage("Convert arguments complete")

# get list of points
lstFIDs = []
lstBlobs = []
rows = arcpy.SearchCursor( inFC )
cnt = 0
for row in rows:
    # if cnt/10 == int(cnt/10):
        # arcpy.AddMessage("Succesfully entered the loop for count #"+str(cnt))
    if cnt >= (startIteration + numIterations):
        break
    if cnt >= startIteration:
        lstFIDs.append(row.FID)   
        blob = row.getValue("Blob")
        lstBlobs.append(blob)
    cnt += 1


arcpy.AddMessage("FIDs: "+str(lstFIDs))

lstCentralityFiles=[]
prevSumRasterFile='None'
count=0    
for FID in lstFIDs:
    where_clause = '"FID" = ' + str(FID)
##    where_clause = "[FID] = " + str(FID)
    
    temppts = "c:\\temp\\temppts.shp"
    inFC2 = arcpy.Select_analysis(inFC, temppts, where_clause)

    rCB = CostBackLink(temppts, rCW, maxCWD)#, rOut + "d" + str(FID) ) #BHM can limit cwd here 999999999
    rFD = Con( rCB > 0, Power(2, (rCB - 1)))
    rFA = FlowAccumulation( rFD, rNL ) #BHM could alter rNL raster here to only include sources within given radius, or greater than given radius
    rFA2=Con(IsNull(rFA),0,rFA)
    centralityFile= rOut + "a" + str(FID) + '.tif'
    rFA2.save (centralityFile)   

    lstCentralityFiles.append(centralityFile)

    count += 1
    sumRasterFile =rOut+"sum"+str(count)+'.tif'
    if count==1:
        sumRaster = Raster(centralityFile)
    else:
        sumRaster=Int(Plus(Raster(prevSumRasterFile),Raster(centralityFile)))
        # outCon=Con((IsNull(sumRaster)),Raster()centralityFile,Con(IsNull(centralityFile),sumRaster,sumRaster+centralityFile))
    sumRaster.save(sumRasterFile)
    delete_data(prevSumRasterFile)
    delete_data(centralityFile)
    prevSumRasterFile=sumRasterFile
    arcpy.AddMessage("Iteration: " + " FID: " + str(FID)+". Saved "+sumRasterFile)

# prevSumRasterFile='None'
# count=0    
# for centralityFile in lstCentralityFiles:
    # count += 1
    # sumRasterFile =rOut+"sum"+str(count)+'.tif'
    # if count==1:
        # sumRaster = Raster(centralityFile)
    # else:
        # sumRaster=Int(Plus(Raster(prevSumRasterFile),Raster(centralityFile)))
        # # outCon=Con((IsNull(sumRaster)),Raster()centralityFile,Con(IsNull(centralityFile),sumRaster,sumRaster+centralityFile))
    # sumRaster.save(sumRasterFile)
    # delete_data(prevSumRasterFile)
    # prevSumRasterFile=sumRasterFile
    
    # print 'saved'+sumRasterFile
    
arcpy.AddMessage('DONE')
 
    

    
