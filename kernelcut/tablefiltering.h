#ifndef _TABLE_FILTERING_H_
#define _TABLE_FILTERING_H_
#include <vector>
#include "gpufiltering.h"
#include "Table2D.h"

Table2D<float> table_filter(const Table2D<vector<float> > & img, const Table2D<bool> & ROI, const Table2D<float> & sigmas_table, const char * filtermode);
Table2D<float> table_filter(const Table2D<vector<float> > & img, const Table2D<bool> & ROI, float sigma);

Table2D<float> table_filter(const Table2D<vector<float> > & img, const Table2D<bool> & ROI, const Table2D<float> & sigmas_table, const char * filtermode)
{
    int img_w = img.getWidth();
    int img_h = img.getHeight();
    int n,N = img_w*img_h;
    // intialization
    float * X = new float[N*3];
    float * values = new float[N];
    float * sigmas = new float[N];
    for(int i=0;i<img_w;i++){
        for(int j=0;j<img_h;j++){
            n = i*img_h + j;
            X[n*3] = img[i][j][0];
            X[n*3+1] = img[i][j][1];
            X[n*3+2] = img[i][j][2];
            if(ROI[i][j])  values[n] = 1.0;
            else values[n] = 0;
            sigmas[n] = sigmas_table[i][j];    
        }
    }
    Table2D<float> outputtable(img_w,img_h,0);
    if((strcmp(filtermode,"p")==0)||(strcmp(filtermode,"q")==0)||(strcmp(filtermode,"m")==0))
    {
        float * output = gpufiltering(X, values, sigmas, N, filtermode[0]);
        for(int i=0;i<img_w;i++){
            for(int j=0;j<img_h;j++){
                n = i*img_h + j;
                outputtable[i][j] += output[n];
            }
        }
        delete [] output;
    }
    else if(strcmp(filtermode,"pq")==0)
    {
        float * output_p = gpufiltering(X, values, sigmas, N, 'p');
        float * output_q = gpufiltering(X, values, sigmas, N, 'q');
        for(int i=0;i<img_w;i++){
            for(int j=0;j<img_h;j++){
                n = i*img_h + j;
                outputtable[i][j] += output_p[n] + output_q[n];
            }
        }
        delete [] output_p;
        delete [] output_q;
    }
    delete [] X;
    delete [] values;
    delete [] sigmas;
    return outputtable;
}

Table2D<float> table_filter(const Table2D<vector<float> > & img, const Table2D<bool> & ROI, float sigma)
{
    int img_w = img.getWidth();
    int img_h = img.getHeight();
    Table2D<float> sigmas_table(img_w,img_h,sigma);
    return table_filter(img, ROI, sigmas_table, "p");
}

#endif
