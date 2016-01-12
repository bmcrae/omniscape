import re

string1 = 'pubet_habitat_270m_sb245_eb65_f20.tif'
string2 = 'one.two.hello.four.five.six.seven'
s1=re.sub(r'_sb\w+_f', '_f', string1)

s2=re.sub(r'^one\.two\.\w+\.four', '', string2)
print s1
print s2

fn = 'puget_habitat270m'
fn2=re.sub(r'^_h.*m', '', fn)
print fn2
blar

import numpy as np
img = np.array([[4, 2, 1, 0.0],[0,0,0,0],[5,0,0,0]])
print img

duh = 3
test=np.where(img >= duh,img,-1)
print 'test'
print test

test=np.where(img <= duh,img,-1)
print test

x = range(0, img.shape[1])
y = range(0, img.shape[0])

(X,Y) = np.meshgrid(x,y)

x_coord = (X*img).sum() / img.sum().astype("float")
y_coord = (Y*img).sum() / img.sum().astype("float")

print x_coord
print y_coord