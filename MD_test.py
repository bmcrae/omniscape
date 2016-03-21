# ---run with arcgis win32 python... that's the only one I have scipy setup for

import numpy as np
import scipy.spatial

#putting non-zero values into vector (kind of thing we'd do with dist function in circ() function)
v3=np.array([[1,3,4,6],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
print v3
rows,cols = np.where(v3)
test=(v3[rows,cols]).flatten()
print test



# IS THIS IT??? Yes.
# http://stackoverflow.com/questions/21003272/difference-between-all-1d-points-in-array-with-python-diff
# a = np.random.random(700)
# print 'a0 a1 a2',a[0],a[1],a[2]
# print a[0]-a[1]
# print a[0]-a[2]
# diff=np.subtract.outer(a,a)[np.tril_indices(a.shape[0],k=-1)]
# print 'diff',diff
# print diff.size

Z1 = np.random.random(700)
Z2 = np.random.random(700)
Var1=np.var(Z1)
Var2=np.var(Z2)

diff1 = np.subtract.outer(Z1,Z1)[np.tril_indices(Z1.shape[0],k=-1)]
diff2 = np.subtract.outer(Z2,Z2)[np.tril_indices(Z2.shape[0],k=-1)]
diff1Sq = np.multiply(diff1,diff1)
diff2Sq = np.multiply(diff2,diff2)
diff1SqDivVar = np.divide(diff1Sq,Var1)
diff2SqDivVar = np.divide(diff2Sq,Var2)

SEucDist = np.sqrt(diff1SqDivVar + diff2SqDivVar)
 
print 'Manual calc'
print SEucDist
print 'Manual mean',np.mean(SEucDist)


Var = [Var1,Var2]
Z=np.zeros((700,2),dtype='float64')
Z[:,0]=Z1
Z[:,1]=Z2

Y1 = scipy.spatial.distance.pdist(Z, 'seuclidean', V=Var)#=None)
Y2 = scipy.spatial.distance.pdist(Z, 'seuclidean', V=None)#=None)
print 'Y1',Y1
print 'Y2',Y2
print 'average Y1, Y2',np.mean(Y1),np.mean(Y2)
print Y1.size

# a = np.array(range(5, 10))
# b = np.array(range(1, 6))

# res = a[np.newaxis, :] - b[:, np.newaxis]
# print 'a',a
# print 'b',b
# print(res)

blarg

# use numpy?
# http://stackoverflow.com/questions/26076576/subtract-all-pairs-of-values-from-two-arrays


# for Var1, get all (x1-y1)
# Diffxy[x,y] = Var1[x]-Var1[y]


# http://stackoverflow.com/questions/17936587/in-numpy-find-euclidean-distance-between-each-pair-from-two-arrays

v1=(np.array([[1,2,3,4],[0,0,0,0],[0,0,0,0],[0,0,0,0]])).transpose()
v2=(np.array([1,3,4,6])).transpose()
v3=np.array([[1,3,4,6],[0,0,0,0],[0,0,0,0],[0,0,0,0]])

test=v1-v3
print v1
print v3

print 'test'
print test
print np.multiply(test,test)
blar
Var = np.array([1.0,3.0])
# Var=Var.transpose()

print v1-v2
print np.multiply((v1-v2),(v1-v2))
distances = (v1-v2)^2
distances = distances.sum(axis=-1)
print 'dist',distances



X=np.zeros((2,4),dtype='int32')
X[0,:]=v1
X[1,:]=v2
X = X.transpose()
print X
Y = scipy.spatial.distance.pdist(X, 'seuclidean', V=Var)#=None)

print Y
