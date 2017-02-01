################################################################################
#          KERNELCUT - software for image segmentation / graph clustering      #
#          http://vision.csd.uwo.ca/code                                       #
################################################################################
import numpy as np
import time
from scipy import sparse
from skimage.color import rgb2lab
from skimage.util import random_noise
from sklearn.neighbors import *
from segutil import *
from matplotlib import pyplot as plt
import maxflow

def NormalizedCutEnergy(A, clustering):
    if isinstance(A, np.ndarray):
        d = np.sum(A, axis=1)
        #print "degree vector"
        #print d
        #print d.shape
        #print type(d)
    elif isinstance(A, sparse.csr.csr_matrix):
        d = A.sum(axis=1)
        #print "degree vector"
        #print d
        #print d.shape
        #print type(d)
    maxclusterid = np.max(clustering)
    #print "max cluster id is: ", maxclusterid
    nassoc_e = 0;
    num_cluster = 0;
    for k in range(maxclusterid+1):
        #print k
        # binary indicators (N-by-1) for cluster k
        S_k = np.array(clustering == k,dtype=np.float)
        #print S_k
        if 0 == np.sum(clustering==k):
             continue # skip empty cluster
        num_cluster = num_cluster + 1
        if isinstance(A, np.ndarray):
            nassoc_e = nassoc_e + np.dot( np.dot(np.transpose(S_k),  A) , S_k) / np.dot(np.transpose(d), S_k)
        elif isinstance(A, sparse.csr.csr_matrix):
            nassoc_e = nassoc_e + np.dot(np.transpose(S_k), A.dot(S_k)) / np.dot(np.transpose(d), S_k)
            nassoc_e = nassoc_e[0,0]
    #print "number of clusters: ", num_cluster
    ncut_e = num_cluster - nassoc_e
    return ncut_e
    
def KernelBound(A, K, current_clustering):
    N = current_clustering.size
    unaries = np.zeros((N, K), dtype=np.float)
    d = A.sum(axis=1)
    for i in range(K):
        #print i
        S_i = np.array(current_clustering == i, dtype=np.float)
        volume_s_i = np.dot(np.transpose(d), S_i)
        volume_s_i = volume_s_i[0,0]
        #print volume_s_i
        temp = np.dot(np.transpose(S_i), A.dot(S_i)) / volume_s_i / volume_s_i
        temp = temp * d
        #print temp.shape
        temp2 = temp + np.reshape( - 2 * A.dot(S_i) / volume_s_i, (N,1))
        #print type(temp2)
        unaries[:,i] = temp2.flatten()
    return unaries

def ConstructKNNGraph(image, KNN_K=400, xyscale=0, samplerate=5):
    print "computing k nearest neighbors..."
    (height, width, channel) = image.shape
    # construct KNN graph in LAB color space. Find KNN_K neighbors and sample the neighbors
    #add gaussian noise to image (Images typically come with quantization artifcat,
    #which makes KNN graph not connected if noise not added)
    noisyimage = random_noise(image,mode='gaussian',var=0.0002);
    labimage = rgb2lab(noisyimage)
    # normalize so that each channel has unit variance
    labimage[:,:,0] /= np.std(labimage[:,:,0])
    labimage[:,:,1] /= np.std(labimage[:,:,1])
    labimage[:,:,2] /= np.std(labimage[:,:,2])
    # compute K nearest neighbor
    X = np.reshape(labimage,(height * width, channel))
    if xyscale > 1e-100:
        print "KNN over color + XY feature"
        rgbxy = np.zeros((height*width, 5), dtype=np.float)
        rgbxy[:,0:3] = X
        rgbxy[:,3] = np.repeat(range(height),width) * xyscale
        rgbxy[:,4] = np.tile(range(width),height).flatten() * xyscale
        X = np.copy(rgbxy)
    start = time.clock()
    tree = KDTree(X, leaf_size=50)
    neighbors = tree.query(X, k=KNN_K,return_distance=False)
    print time.clock() - start, "seconds process time"
 
    # sample the neighbors
    neighbors = neighbors[:,0::5]
    neighbors = perturbknn(neighbors, height, width)
    row = np.repeat(range(height*width), KNN_K/samplerate)
    col = neighbors.flatten()
    data = np.ones(height*width*KNN_K/samplerate, dtype=np.float)
    affinity_matrix = sparse.csr_matrix((data, (row, col)), shape=(height*width,height*width),dtype=np.float) 
    affinity_matrix = (affinity_matrix + affinity_matrix.transpose(True)) /2
    return affinity_matrix

def KernelCut(image, opt):
    [height, width, channel] = image.shape
    top = opt.box[0]
    left = opt.box[1]
    boxheight = opt.box[2]
    boxwidth = opt.box[3]
    connecttype = 8
    gridvariance = computeGridVariance(np.asarray(image), connecttype)
    #print "Grid variance is: ", gridvariance
    
    # intial segmentation
    initsegmentation = np.zeros([height, width],dtype=int)
    initsegmentation[(top-1):(top+boxheight-2),(left-1):(left+boxwidth-2)] = 1
    # hard constraints
    hardconstraints = np.zeros([height, width], dtype=int)
    hardconstraints[(top-1):(top+boxheight-2),(left-1):(left+boxwidth-2)] = -1
    print "hardconstraintsflag", opt.hardconstraintsflag
    plt.figure(2)
    plt.imshow(maskoncolorimage(image, initsegmentation))
    plt.show(block=False)
    
    affinity_matrix = ConstructKNNGraph(image, opt.KNN_K, opt.xyscale)

    energy = NormalizedCutEnergy(affinity_matrix, initsegmentation.flatten())
    print "normalized cut energy is ", energy
    
    segmentation = np.copy(initsegmentation)
    itr_num = 1
    while True:
        itr_num = itr_num +1
        start = time.time()
        g = maxflow.Graph[float](height*width, height*width*4)
        nodeids = g.add_grid_nodes((height,width))
        addsmoothness(g, opt.weight_smoothness, image, connecttype, gridvariance)
        unaries = KernelBound(affinity_matrix, 2, segmentation.flatten())
        capsource = np.reshape(unaries[:,0],(height,width))
        capsink = np.reshape(unaries[:,1],(height,width))
        # hard constraints?
        if opt.hardconstraintsflag:
            capsink[hardconstraints==0] = 1e+100
            capsource[hardconstraints==1] = 1e+100
        
        g.add_grid_tedges(nodeids, capsink, capsource)
        g.maxflow()
        newsegmentation = np.asarray(g.get_grid_segments(nodeids), dtype=np.int)
        print "# of pixels with new label", np.sum(newsegmentation.flatten()!=segmentation.flatten())
        if np.sum(newsegmentation.flatten()!=segmentation.flatten()) >10:
            segmentation = np.copy(newsegmentation)
            energy = NormalizedCutEnergy(affinity_matrix, segmentation.flatten())
            end = time.time()
            print "normalized cut energy is ", energy
        else:
            print "converged"
            break
        if itr_num>opt.MAX_ITERATION:
            print "max number of iteration"
            break;
        #print end - start
    return segmentation

