# TO DO 7/22/16:

# comment code, functions
# move kluges into functions, e.g., trimming raster to numpy arrays
# check if there is easy way to do bands in chunks of stripes for memory footprint. Right now need entire band in numpy array

# put individual arcpy steps into functions to allow smoother migration away from arcpy later?


# GDAL looks like a huge pain. Reasonable for Circuitscape, but not to expect OmniScape users to deal with UNLESS COMPILED.
# #trying pytables. haing trouble with 'valid win32 application' in arcgis install, so trying in basic c:\python27\scripts install:
# pip install tables-3.2.3-cp27-cp27m-win_amd64.whl
# usual problems with not a supported wheel on this platform...
# bailin for now. try numpy memmap? problem is getting the raster INTO memmap... would need to read portions of it and save the portions to the memmap.
#maybe arcpy is the way to go at least for now.


# need nanmax
# Fix brightbug for sources and targets
# To eliminate arcpy, challenges are:
# 1) size of grids- manipulating and accumulating massive ones. 
# 2) cost distance and flow accumulation

# 4/19/16
# major problem: climate and non-climate differ for TO and Non-TO when have source ras
# climate not designed for source ras. uses target tally
# need to have a target tally that acts like target sum. 

# Add text for artifact correction
