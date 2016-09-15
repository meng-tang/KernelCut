Input image frames: directory "motionsegmentation/ducks01/images"
Initial Strokes for the first frame: directory "motionsegmentation/ducks01/seedsmulti"
Output segmentation: directory "motionsegmentation/ducks01/output"
1. build dependency libraries (maxflow and easybmp)
cd libs
make all
2. build main program
cd ../kernelcut
make main
4. download executable for optical flow
wget http://lmb.informatik.uni-freiburg.de/resources/binaries/pami2010Linux64.zip
unzip pami2010Linux64.zip -d libs/LDOF
3. compute optical flow
cd ../
matlab -nojvm -nosplash -nodisplay -r "cd motionsegmentation/ducks01; computeopticalflow"
4. convert from RGB space to LAB space
matlab -nojvm -nosplash -nodisplay -r "cd motionsegmentation/ducks01; getsubpixelimages"
5. compute KNN graph for joint LAB + XY + M space
matlab -nojvm -nosplash -nodisplay -r "cd motionsegmentation/ducks01; computeknn"
6. go to motionsegmentation/motion.sh, change codepath, and run script
chmod +x ./motionsegmentation/motion.sh
./motionsegmentation/motion.sh
