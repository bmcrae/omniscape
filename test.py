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