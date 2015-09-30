
import numpy as np
import scipy as sp
import scipy.stats as stats

import copy
import sys
def quantile(y, k=4):
    """
    Calculates the quantiles for an array
    Parameters
    ----------
    y : array
        (n,1), values to classify
    k : int
        number of quantiles
    Returns
    -------
    implicit  : array
                (n,1), quantile values
    Examples
    --------
    >>> x = np.arange(1000)
    >>> quantile(x)
    array([ 249.75,  499.5 ,  749.25,  999.  ])
    >>> quantile(x, k = 3)
    array([ 333.,  666.,  999.])
    >>>
    Note that if there are enough ties that the quantile values repeat, we
    collapse to pseudo quantiles in which case the number of classes will be
    less than k
    >>> x = [1.0] * 100
    >>> x.extend([3.0] * 40)
    >>> len(x)
    140
    >>> y = np.array(x)
    >>> quantile(y)
    array([ 1.,  3.])
    """
    w = 100. / k
    p = np.arange(w, 100 + w, w)
    if p[-1] > 100.0:
        p[-1] = 100.0
    q = np.array([stats.scoreatpercentile(y, pct) for pct in p])
    return np.unique(q)


y = [1,2,3,4,5,6,7,8,9,10,10,11]
k=4
x=quantile(y,k)
print x




import numpy as npy
inArray=npy.random.rand(20,20)
inArray[0:3] = -9999
inArray[4:8] = 0

def seq(start, stop, step=1):
    n = int(round((stop - start)/float(step)))
    if n > 1:
        return([start + step*i for i in range(n+1)])
     

     
numQuantiles = 100
interval = 100/numQuantiles
breakSequence = seq(interval, 100-interval, interval)

# set up quantile array
quantileArray = npy.zeros_like(inArray, dtype='int32')        
quantileArray = npy.where(inArray > 0,numQuantiles,quantileArray)
           
inVector = inArray.flatten()
ind=npy.where(inVector > 0) 
inVector = inVector[ind] #Eliminate nodata and zeros

# quantilize                      
print 'Partitioning ' + str(len(inVector)) + ' non-zero values into ' + str(numQuantiles) + ' quantiles.'
quantileBreaks = npy.percentile(inVector, breakSequence)

for i in range(numQuantiles-2,-1,-1):
    print 'quantile:',breakSequence[i]
    print 'break: ',quantileBreaks[i]
    quantileArray = npy.where (inArray<quantileBreaks[i],i,quantileArray)
    
quantileArray = npy.where(inArray < 0,-9999,quantileArray)
quantileArray = npy.where(inArray == 0,0,quantileArray)
print quantileArray

















asldkfj

# breakArray2 =[20,40,60,80]
# print breakArray2
quantileBreaks = npy.percentile(vec, breakSeq)
print quantileBreaks
print len(quantileBreaks)
for i in range(numQuantiles-2,-1,-1):
    print quantileBreaks[i]
brea
for qbreak in quantileBreaks:
    print qbreak

# 5
# 20, 40, 60, 80
# 10
# 10, 20 30 40 50 60 70 80 90
# 100
# 1--99

# interval = 100/numQuantiles
# breakArray = arange(interval, 100-interval, interval)

# 100/numQuant, 2*

# test = arange(0,numQuantiles-1,
# [start,] stop[, step,][, dtype])

def calc_resil_quantile_by_zone(resilArray, zoneArray, numQuantiles):
    try:
        zones = get_unique(zoneArray)
        
        resilClassArray = zeros(resilArray.shape, dtype='int32')
        resilClassArray = npy.where(resilArray < 0, -9999, resilClassArray)
        pctDone = 0
        x = 0
        print '0 percent done'
        for zone in zones:              
            # print'zone'+str(zone) + ' out of ' + str(len(zones))
            resilArrayZone = npy.where(zoneArray == zone, resilArray, -9999)
            resilVector = resilArrayZone.flatten()
            
            ind=where(resilVector > -1) 
            resilVector = resilVector[ind] #Eliminate nodata
            if len(resilVector) == 0:
                print 'zone=',zone
                print 'non'
                resilClassArrayZone =  resilArrayZone # set to nodata
            else:
                [qb1,qb2,qb3,qb4]=percentile(resilVector, [20,40,60,80])
                # print 'Quintile Breaks: ' + str([qb1,qb2,qb3,qb4])
                pctDone = report_pct_done(x, len(zones), pctDone)
                x = x+1
                
                resilClassArrayZone = npy.where (resilArrayZone,5,5)
                resilClassArrayZone = npy.where (resilArrayZone<qb4,4,resilClassArrayZone)
                resilClassArrayZone = npy.where (resilArrayZone<qb3,3,resilClassArrayZone)
                resilClassArrayZone = npy.where (resilArrayZone<qb2,2,resilClassArrayZone)
                resilClassArrayZone = npy.where (resilArrayZone<qb1,1,resilClassArrayZone)
                resilClassArrayZone = npy.where (resilArrayZone<0,-9999,resilClassArrayZone) 
           
            resilClassArray = npy.where(zoneArray == zone, resilClassArrayZone, resilClassArray)
            del resilClassArrayZone
            del resilArrayZone
            gc.collect()
        return resilClassArray
    except:
        pass