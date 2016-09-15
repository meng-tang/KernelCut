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
```
Download executable for optical flow
```{r, engine='bash'}
wget http://lmb.informatik.uni-freiburg.de/resources/binaries/pami2010Linux64.zip
unzip pami2010Linux64.zip -d libs/LDOF
```
Compute optical flow
```{r, engine='bash'}
cd ../
matlab -nojvm -nosplash -nodisplay -r "cd motionsegmentation/ducks01; computeopticalflow"
```
Compute KNN graph for joint LAB + XY + M space
\#Convert from RGB space to LAB space
matlab -nojvm -nosplash -nodisplay -r "cd motionsegmentation/ducks01; getsubpixelimages"
```{r, engine='bash'}
matlab -nojvm -nosplash -nodisplay -r "cd motionsegmentation/ducks01; computeknn"
```
Go to motionsegmentation/motion.sh, change codepath, and run script
```{r, engine='bash'}
chmod +x ./motionsegmentation/motion.sh
./motionsegmentation/motion.sh
```
Output segmentations are directory "motionsegmentation/ducks01/output".
