################################################################################
#          KERNELCUT - software for image segmentation / graph clustering      #
#          http://vision.csd.uwo.ca/code                                       #
################################################################################
#!/usr/bin/python
import numpy as np
from scipy import sparse


def addsmoothness(g, weight, image, connecttype, gridvariance):
    if weight < 1e-100:
        return
    pixelshifts = np.asarray([[1,0],[0,1],[1,1],[1,-1],])
    if not (connecttype in [4,8]):
        print "Grid connect type should be 4 or 8"
        exit(-2)
    [h, w, c] = image.shape
    for i in range(connecttype/2):
        shift_x = pixelshifts[i,0]
        shift_y = pixelshifts[i,1]
        imageshifted = image[max(0,(0+shift_y)):min(h,(h+shift_y)),
                             max(0,(0+shift_x)):min(w,(w+shift_x)),:]
        nodeids = np.reshape(range(h*w), (h,w))
        if shift_x>0:
            imagecrop = image[:, 0:-shift_x, :]
        else:
            imagecrop = image[:, -shift_x:, :]
        if shift_y>0:
            imagecrop = imagecrop[0:-shift_y, :, :]
        else:
            imagecrop = imagecrop[-shift_y:, :, :]
        if shift_x>0:
            nodeidscrop = nodeids[:, 0:-shift_x]
        else:
            nodeidscrop = nodeids[:, -shift_x:]
        if shift_y>0:
            nodeidscrop = nodeidscrop[0:-shift_y, :]
        else:
            nodeidscrop = nodeidscrop[-shift_y:, :]
        diff = imageshifted.astype(float) - imagecrop.astype(float)
        diff = diff*diff
        edgeweight = np.exp( - np.sum(diff,axis=2) / 2 / gridvariance) / np.sqrt(shift_x*shift_x + shift_y*shift_y) * weight
        structure = np.zeros((3,3))
        structure[shift_y+1, shift_x+1] = 1
        g.add_grid_edges(nodeidscrop, weights=edgeweight, structure=structure, symmetric=True)
        
def maskoncolorimage(colorimage, mask, bkgcolor=[255,255,255]):
    maskedimage = np.copy(colorimage)
    for i in range(3):
        onechannel = maskedimage[:,:,i]
        onechannel[mask==0] = bkgcolor[i]
        maskedimage[:,:,i] = onechannel
    return maskedimage    

def computeGridVariance(image, connecttype):
    pixelshifts = np.asarray([[1,0],[0,1],[1,1],[1,-1],])
    if not (connecttype in [4,8]):
        print "Grid connect type should be 4 or 8"
        exit(-2)
    [h, w, c] = image.shape
    counter = int(0)
    totalvariance = 0
    for i in range(connecttype/2):
        shift_x = pixelshifts[i,0]
        shift_y = pixelshifts[i,1]
        imageshifted = image[max(0,(0+shift_y)):min(h,(h+shift_y)),
                             max(0,(0+shift_x)):min(w,(w+shift_x)),:]
        if shift_x>0:
            imagecrop = image[:, 0:-shift_x, :]
        else:
            imagecrop = image[:, -shift_x:, :]
        if shift_y>0:
            imagecrop = imagecrop[0:-shift_y, :, :]
        else:
            imagecrop = imagecrop[-shift_y:, :, :]
        diff = imageshifted.astype(float) - imagecrop.astype(float)
        totalvariance += np.sum(diff*diff)
        counter += diff.size / 3
    return totalvariance / counter

# perturb the knn in 'image grid', this is a trick to reach out to large
# neighborhood in the color space without knn search of large K
def perturbknn(neighbors, height, width):
    knnidx_c = np.remainder(neighbors,width);
    knnidx_r = (neighbors-knnidx_c)/width;
    knnidx_r = knnidx_r + np.random.randint(-2, 3, size=neighbors.shape)
    knnidx_c = knnidx_c + np.random.randint(-2, 3, size=neighbors.shape)
    knnidx_r[knnidx_r<0] = 0
    knnidx_r[knnidx_r>=height] = height-1
    knnidx_c[knnidx_c<0] = 0
    knnidx_c[knnidx_c>=width] = width-1
    return knnidx_r * width + knnidx_c




