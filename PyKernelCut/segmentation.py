################################################################################
#          KERNELCUT - software for image segmentation / graph clustering      #
#          http://vision.csd.uwo.ca/code                                       #
#          https://github.com/meng-tang/KernelCut                              #
#          Meng Tang (mtang73@uwo.ca) (contact author)                         #
#          Dmitrii Marin (dmitrii.a.marin@gmail.com)                           #
#          Ismail Ben Ayed (ismail.benayed@etsmtl.ca)                          #
#          Yuri Boykov (yuri@csd.uwo.ca)                                       #
################################################################################
#!/usr/bin/python 

import sys, getopt
import os, time
import skimage.io as io
import numpy as np
from optparse import OptionParser
import timeit

from segutil import *
from kernelcut import *
from scipy import sparse
from matplotlib import pyplot as plt


def parsearguments():
    parser = OptionParser()
    parser.add_option("-i", "--image", dest="imagename",
                     help="name of the image, can also be an URL")
    parser.add_option("-d", "--dir", dest="imagedir",
                     help="directory that contains image, defalt is .py file path", type="string", default=os.path.dirname(os.path.realpath(__file__)))
    parser.add_option("-b", "--box", dest="box",
                     help="bounding box (top, left, height, width)",type="int", nargs=4)
    parser.add_option("--hard", dest="hardconstraintsflag",
                     help="enforce hard constraints, default off",action="store_true",default=False)
    parser.add_option("-k", "--knn", dest="KNN_K",
                     help="K nearest neighbors, defualt 400",type="int",default=400)
    parser.add_option("-s", "--smoothness", dest="weight_smoothness",
                     help="weight of smoothness term, defualt 0.0001", type="float",default=0.0001)
    parser.add_option("--maxitr", dest="MAX_ITERATION",
                     help="maximum number of iteration, defualt 80", type="int",default=80)
    parser.add_option("--xyscale", dest="xyscale",
                     help="scale of x,y coordinate for KNN graph, defualt 0", type="float",default=0)

    (opt, args) = parser.parse_args()
    print 'image directory is ', opt.imagedir
    print 'image name is ', opt.imagename
    return opt

def main():
    opt = parsearguments()
    if "http" in opt.imagename:
        image = io.imread(opt.imagename)
    else:
        image = io.imread(os.path.join(opt.imagedir, opt.imagename))
    (height, width, channel) = image.shape
    plt.figure(1)
    plt.imshow(image)
    plt.show(block=False)

    segmentation = KernelCut(image, opt)
    plt.figure(3)
    plt.imshow(maskoncolorimage(image, segmentation))
    plt.show()
    return

if __name__ == "__main__":
    main()
