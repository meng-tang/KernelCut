This is the code for the paper:

	"Normalized Cut Meets MRF"
	Meng Tang, Dmitrii Marin, Ismail Ben Ayed, Yuri Boykov
	In European Conference on Computer Vision (ECCV), Amsterdam, the Netherlands, October, 2016

The CORE of our algorithm is linearization or unary bound for Normalized Cut (NC).<br />
Simple implementation of such linearization is given in a FEW lines in "matlab/KernelBound.m". <br />
The function **[ unaries ] = KernelBound( A, K, current_clustering)** simply takes affinity, cluster number and current clustering and gives unary terms.<br />
Example of optimizing NC or AA (avearge association) ONLY is in "matlab/syntheticclustering.m". Below is sample result with NC:<br />
<span><img src="matlab/NC_init.png" alt="" width="400"/>
<img src="matlab/NC_clustering.png" alt="" width="400"/>
<img src="matlab/NC_energy.png" alt="" width="400"/></span>

## Motion Segmentation using KernelCut ##
Input image frames: directory "motionsegmentation/ducks01/images"  
Initial Strokes for the first frame: directory "motionsegmentation/ducks01/seedsmulti"  

Build dependency libraries (maxflow and easybmp)  
```{r, engine='bash'}
cd libs
make all
```
Build main program
```{r, engine='bash'}
cd ../kernelcut
make main
cd ../
```
Download executable for optical flow
```{r, engine='bash'}
wget http://lmb.informatik.uni-freiburg.de/resources/binaries/pami2010Linux64.zip
unzip pami2010Linux64.zip -d libs/LDOF
```
Compute optical flow
```{r, engine='bash'}
chmod +x libs/LDOF/ldof
matlab -nojvm -nosplash -nodisplay -r "cd motionsegmentation/ducks01; computeopticalflow; exit()"
```
Compute KNN graph for joint LAB + XY + M space
```{r, engine='bash'}
matlab -nojvm -nosplash -nodisplay -r "cd motionsegmentation/ducks01; getsubpixelimages; exit();"
matlab -nojvm -nosplash -nodisplay -r "cd motionsegmentation/ducks01; computeknn; exit()"
```
(Visualization of KNN graph is by clicking on image pixel, simply run motionsegmentation/visualizeknnbyclick.m)

Go to motionsegmentation/motion.sh, change codepath, and run script
```{r, engine='bash'}
chmod +x ./motionsegmentation/motion.sh
./motionsegmentation/motion.sh
```
Output segmentations are in the directory "motionsegmentation/ducks01/output".

(note that if initialized from seeds, the colors has to be of the following: {white,red,blue,green,black,navy})

## Image Clustering using Kernel Cut and Spectral Cut##
Download images:  

    cd imageclustering
    wget http://people.csail.mit.edu/torralba/code/spatialenvelope/spatial_envelope_256x256_static_8outdoorcategories.zip
    unzip spatial_envelope_256x256_static_8outdoorcategories.zip -d images
    mv images/spatial_envelope_256x256_static_8outdoorcategories/* images/
    
Download feature extractor:  

    wget http://vision.stanford.edu/projects/objectbank/MATLAB_release.zip
    unzip MATLAB_release.zip -d ./
    mv MATLAB_release object-bank
Change the variable "modelpath" in getFeatureOB.m accordingly and move the file:  

    mv getFeatureOB.m object-bank/code/partless/
    
Extract features and compute Gaussian kernel:  

    matlab -nojvm -nosplash -nodisplay -r "computegaussiankernel; exit()"
Compute eigenvalues:  

    matlab -nojvm -nosplash -nodisplay -r "eigen_labelme; exit()"
Prepare ground truth labels:  

    matlab -nojvm -nosplash -nodisplay -r "preparelabels; exit()"
Spectral clustering:  

    matlab -nojvm -nosplash -nodisplay -r "sc_labelme; exit()"
Spectral cut:  

    matlab -nojvm -nosplash -nodisplay -r "sc_labelme; exit()"
Kernel k-means:  

    matlab -nojvm -nosplash -nodisplay -r "kkm_labelme; exit()"
Kernel cut:  

    matlab -nojvm -nosplash -nodisplay -r "kcut_labelme; exit()"
These matlab scripts will report NMI (normalized mutual information) values for clustering obtained.

