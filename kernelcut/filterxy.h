#ifndef _H_FILTERXY_
#define _H_FILTERXY_
#include <math.h>
#include <vector>
#include "basicutil.h"
#include "permutohedral.h"

Table2D<float> filterrgb(const Table2D<vector<float> > image, Table2D<bool> ROI, double * sigmargb);
Table2D<float> filterrgbxy(const Table2D<vector<float> > image, Table2D<bool> ROI, double * sigmargb, double sigmax, double sigmay);

Table2D<float> filterrgb(const Table2D<vector<float> > image, Table2D<bool> ROI, double * sigmargb){
	int img_w = image.getWidth();
	int img_h = image.getHeight();
	Table2D<float> returnv(img_w,img_h,0);
	float * pos;
	float * val;
	float * out;
	val = new float[img_w*img_h];
	out = new float[img_w*img_h];

	pos = new float[img_w*img_h*3];
	for(int i=0;i<img_w;i++){
		for(int j=0;j<img_h;j++){
			int idx = i*img_h + j;
			pos[idx*3+0] = (float)(image[i][j][0]) / sigmargb[0];
			pos[idx*3+1] = (float)(image[i][j][1]) / sigmargb[1];
			pos[idx*3+2] = (float)(image[i][j][2]) / sigmargb[2];
			if(ROI[i][j])
				val[idx] = 1.0;
			else
				val[idx] = 0;
		}
	}
	PermutohedralLattice::filter(pos, 3, val, 1, img_w*img_h, out);
	delete [] pos;
	delete [] val;
	for(int i=0;i<img_w;i++){
		for(int j=0;j<img_h;j++){
			int idx = i*img_h + j;
			returnv[i][j] = out[idx]/15 ;
		}
	}
	delete [] out;
	//outv(returnv.getMax());
	//outv(returnv.getMin());
	//outv(returnv.getMean());
	return returnv;
}

Table2D<float> filterrgbxy(const Table2D<vector<float> > image, Table2D<bool> ROI, double * sigmargb, double sigmax, double sigmay){
	int img_w = image.getWidth();
	int img_h = image.getHeight();
	Table2D<float> returnv(img_w,img_h,0);
	float * pos;
	float * val;
	float * out;
	val = new float[img_w*img_h];
	out = new float[img_w*img_h];
	pos = new float[img_w*img_h*5];
	for(int i=0;i<img_w;i++){
		for(int j=0;j<img_h;j++){
			int idx = i*img_h + j;
			pos[idx*5+0] = (float)(image[i][j][0]) / sigmargb[0];
			pos[idx*5+1] = (float)(image[i][j][1]) / sigmargb[1];
			pos[idx*5+2] = (float)(image[i][j][2]) / sigmargb[2];
			pos[idx*5+3] = ((float)(i)) / sigmax;
			pos[idx*5+4] = ((float)(j)) / sigmay;
			if(ROI[i][j])
				val[idx] = 1.0;
			else
				val[idx] = 0;
		}
	}
	PermutohedralLattice::filter(pos, 5, val, 1, img_w*img_h, out);
	delete [] pos;
	delete [] val;
	for(int i=0;i<img_w;i++){
		for(int j=0;j<img_h;j++){
			int idx = i*img_h + j;
			returnv[i][j] = out[idx] / 35;
		}
	}
	delete [] out;
	return returnv;
}
#endif
