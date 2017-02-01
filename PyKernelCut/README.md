## KernelCut in Python ##
Python implementation of KernelCut for binary segmentation is provided.
Example Usage:
```{r, engine='bash'}
$ python segmentation.py -h
Usage: segmentation.py [options]

Options:
  -h, --help            show this help message and exit
  -i IMAGENAME, --image=IMAGENAME
                        name of the image, can also be an URL
  -d IMAGEDIR, --dir=IMAGEDIR
                        directory that contains image, defalt is .py file path
  -b BOX, --box=BOX     bounding box (top, left, height, width)
  --hard                enforce hard constraints, default off
  -k KNN_K, --knn=KNN_K
                        K nearest neighbors, defualt 400
  -s WEIGHT_SMOOTHNESS, --smoothness=WEIGHT_SMOOTHNESS
                        weight of smoothness term, defualt 0.0001
  --maxitr=MAX_ITERATION
                        maximum number of iteration, defualt 80
  --xyscale=XYSCALE     scale of x,y coordinate for KNN graph, defualt 0

$ python segmentation.py -i 124084.jpg -b 10 20 300 300 -k 100
$ python segmentation.py -i 0_5_5303.bmp -b 50 100 200 200 -s 0
$ python segmentation.py -i 0_5_5303.bmp -b 50 100 200 200
$ python segmentation.py -i 314016.jpg -b 10 80 300 300 --hard

```
Note that To use PyKernelCut several dependencies (skimage, scipy, [PyMaxflow](https://github.com/pmneila/PyMaxflow)) have to be installed.
