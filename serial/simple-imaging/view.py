import matplotlib
import matplotlib.pylab

import numpy

img = numpy.fromfile('out.image', dtype=float)
theta= 0.08
l=300000

img_s = theta*l
#img_s =24000
img = img.reshape(img_s,img_s)

matplotlib.pylab.imshow(img[6000:18000,6000:18000])

matplotlib.pylab.show()
