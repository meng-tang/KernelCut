#include <stdlib.h>
#include <stdio.h>
#include <iostream>
using namespace std;
#include "dirent.h"

double arrayMean(double * vs, int n);
vector<Vect3D> getvectordata(const Table2D<Vect3D> & floatimg, Table2D<bool> ROI);

double arrayMean(double * vs, int n){
    double sum=0;
    for(int i=0;i<n;i++)
        sum += vs[i];
    return sum / n;
}

int countFilesInDirectory(const char * dir){
    DIR *dpdf;
    struct dirent *epdf;
    dpdf = opendir(dir);
    if(dpdf == NULL) return -1;
    int filecount=0;
    while(epdf = readdir(dpdf)){
        char filename[10];
        strcpy(filename,epdf->d_name);
        if(strlen(filename)<=2) continue;
        filecount++;
    }
    return filecount;
}

vector<Vect3D> getvectordata(const Table2D<Vect3D> & floatimg, Table2D<bool> ROI){
    int N = countintable(ROI,true);
    vector<Vect3D> data(N);
    int img_w = floatimg.getWidth();
    int img_h = floatimg.getHeight();
    int id=0;
    for(int i=0;i<img_w;i++){
        for(int j=0;j<img_h;j++){
            if(ROI[i][j])
                data[id++] = floatimg[i][j];
        }
    }
    return data;
}
