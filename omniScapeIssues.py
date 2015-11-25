priorities
clean up 
maxcur, maxflow to adjust for components
voltages

HOW TO DEAL WITH components/strengths
Don't connect anything > radius
if rad = 100km
water=400
1km water = 400km cwd


Costal issue
noweight doesn't help
dividing current by focal sum doesn't really help- tip of peninsula lights up, but only because it has lots of nodata around it. shorelines don't really.

what if sources and targs in diff components? does that mess up multiplier?
    need to run components separately?
    remove nodata?
    look at current flowing to ground to determine nsources? 
    then divide current by this number, adjst nsources accordingly to get multiplier.

what if center ground is in nodata?
what if targets are not all in same component?

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

