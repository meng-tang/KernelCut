#ifndef _KNN_H_
#define _KNN_H_
#include <algorithm>

void perturbconnection(Table2D<int> & knntable,int img_w,int img_h){
	int knn_w = knntable.getWidth();
	int knn_h = knntable.getHeight();
	for(int i=0;i<knn_w;i++){
	    for(int j=0;j<knn_h;j++){
			int x = (knntable[i][j]) / img_h;
			int y = (knntable[i][j]) % img_h;
			x = x+rand()%9-4;
			y = y+rand()%9-4;
			x = min(x,img_w-1);x=max(x,0);
			y = min(y,img_h-1);y=max(y,0);
			knntable[i][j] = x*img_h + y ;
		}
	}
}

Table2D<double> zeroonekernel(const Table2D<int> & knntable, Table2D<bool> ROI, int KNN_K){
	int img_w = ROI.getWidth();
	int img_h = ROI.getHeight();
	Table2D<double> returnv(img_w,img_h,0);
	for(int i=0;i<img_w;i++){
	    for(int j=0;j<img_h;j++){
			int p_idx = j+i*img_h;
		    for(int k=0;k<KNN_K;k++){
			    int q_idx = knntable[p_idx][k];
				if(ROI[q_idx/img_h][q_idx%img_h]){
				    returnv[i][j] += 1.0;
				}
			}
		}
	}
	for(int i=0;i<img_w;i++){
	    for(int j=0;j<img_h;j++){
			int q_idx = j+i*img_h;
			if(ROI[i][j]){
				for(int k=0;k<KNN_K;k++){
					int p_idx = knntable[q_idx][k];
					returnv[p_idx/img_h][p_idx%img_h] += 1.0;
				}
			}
		}
	}
	return returnv;
}
#endif
