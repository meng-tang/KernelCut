Input image frames: directory "motionsegmentation/ducks01/images"  
Initial Strokes for the first frame: directory "motionsegmentation/ducks01/seedsmulti"  

<ol start="0">
  <li>build dependency libraries (maxflow and easybmp)
```{r, engine='bash'}
cd libs
make all
```</li>
  <li>item 4</li>
</ol>

1.  build dependency libraries (maxflow and easybmp)  
```{r, engine='bash'}
cd libs
make all
```
2.  build main program
```{r, engine='bash'}
cd ../kernelcut
make main
```
3.  download executable for optical flow
```{r, engine='bash'}
wget http://lmb.informatik.uni-freiburg.de/resources/binaries/pami2010Linux64.zip
unzip pami2010Linux64.zip -d libs/LDOF
```
4.  compute optical flow
```{r, engine='bash'}
cd ../
matlab -nojvm -nosplash -nodisplay -r "cd motionsegmentation/ducks01; computeopticalflow"
```
5.  convert from RGB space to LAB space
```{r, engine='bash'}
matlab -nojvm -nosplash -nodisplay -r "cd motionsegmentation/ducks01; getsubpixelimages"
```
6.  compute KNN graph for joint LAB + XY + M space
```{r, engine='bash'}
matlab -nojvm -nosplash -nodisplay -r "cd motionsegmentation/ducks01; computeknn"
```
7.  go to motionsegmentation/motion.sh, change codepath, and run script
```{r, engine='bash'}
chmod +x ./motionsegmentation/motion.sh
./motionsegmentation/motion.sh
```
Output segmentations are directory "motionsegmentation/ducks01/output".
