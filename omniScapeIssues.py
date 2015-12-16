how to do center of mass (maybe best not to do this, keep flow going to area in center of block):
import numpy as np
img = np.array([[4, 2, 1, 0.0],[0,0,0,0],[5,0,0,0]])
print img

x = range(0, img.shape[1])
y = range(0, img.shape[0])

(X,Y) = np.meshgrid(x,y)

x_coord = (X*img).sum() / img.sum().astype("float")
y_coord = (Y*img).sum() / img.sum().astype("float")

print x_coord
print y_coord


add fade to FA?

Note: running with total fade (distance=radius) at blocksize 15 closely apprximates results from blocksize=1, at least with 
disteq. Actually, bs15 results without disteq apprximate results from bs1 disteq anyway. Interesting question about 
how total fade actually affects results though. short distance movements will be short-shrifted.


run one more analogous to the totfade but with blocksize fade.
what about source strength declines iwth square of distance, so that when combined with totfade, it comes out in the wash?
radius-dist/radius


priorities
clean up 

Add option to do max cur along with sum?

voltages
        # fixme: not completed. Aim is to 
        # see how much improvement is possible if r reduced to 1.
        # vdiff * (r-1)/r... if r=1, no improvement possible. if r=2, half of voltage could be reduced. etc.
        # need r's, will include r from 3 pixels. but try to calc how much v would drop if that one pixel were restored.
        # HOW MUCH WOULD V DROP ACROSS PIXEL IF R>1
        
        # But know current, so v=ir to get r?
        10 amps going through cell. max drop = x. assume

limit cur by cwd? would require numpytoras, then con.

# maxcur, maxflow to adjust for components DONE I THINK
# HOW TO DEAL WITH components/strengths
# Don't connect anything > radius
# if rad = 100km
# water=400
# 1km water = 400km cwd

Coastal issue
noweight doesn't help
dividing current by focal sum doesn't really help- tip of peninsula lights up, but only because it has lots of nodata around it. shorelines don't really.

what if sources and targs in diff components? does that mess up multiplier?
    need to run components separately?
    remove nodata?
    look at current flowing to ground to determine nsources? 
    then divide current by this number so that sums to 1, adjust nsources accordingly to get multiplier:
    nsources=nsources*maxcur
    
    or...
    sourceCorrection= maxcur
    if noweight, divide by sourceCorrection
    if targonly, ""
    else, do nothing (washes out)
    
what if center ground is in nodata?
    rare since water isn't nodata
what if targets are not all in same component?
    assume rare for now

uses NANs instead of -999
# drop fade from flow?
# try flow without rand... would be nice for results to be consistent across runs

#try difference between run with some sources missing... does pattern look right? 
#Yes.

# do a donut on flow? just show long distance? 
    # probably still get spaghetti. and more elegant to have same settings as current
    
# try to use inverse of FA as resistance. Could consolidate lines?
    # it does consolidate, but not very predictable- suggest just doing major lines instead.


# subtract out sources from FA grid? could eliminate whiskers.
    # this line: rFA = arcpy.sa.Plus(arcpy.sa.FlowAccumulation( rFD, sourceRasFA ), sourceRasFA) #Note! added in source strengths 8/15/15 
    #doesn't really help, can just get rid of these later anyway
    # don't do int before band saves for flow, unless speed at issue
# save int and intln flow files


# Flow mapping
# break out values into 10 classes
# grow each class by value. Nah.
# Have separate line layers for each class? Easier to have   one on top of the other, transparancy etc.


# go with fade... setting whole block to average of edges has its own artifacts.
# Based on 25km comparisons, full block size fade distance works much better than half. Matches BS1 data well.


# What if using distfn- does this mean targets in block that are not in center cell will have less weight?
# either ignore or set center block source modifiers to 1
# OK- targets are already 1, modify only works on sources.


# write and read numpy grids?
# adddata for FA?

