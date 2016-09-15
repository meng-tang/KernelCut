#include <stdio.h>
#include <iostream>
#include <ctime>
#include <string>

#include <graph.h>
#include <EasyBMP.h>
#include <dirent.h>
#include "basicutil.h"
#include "Filtering.h"
#include "filterxy.h"
#include "tablefiltering.h"
#include "PPBCadp.h"

#define MODE_ADPKERNEL 4
#define INIT_RANDSEEDS 3
#define INIT_BOX 1
#define INIT_STROKES 5

#define MAPPING_MODE_CONST 1
#define MAPPING_MODE_IDENTITY 2
#define MAPPING_MODE_POWER 3
#define MAPPING_MODE_TRUNCATION 4

#define N_IMG 9

using namespace std;
void kernelsegmentation(int argc,char *argv[], int init_mode,double * returnv=NULL,double * energyv=NULL);

// global variables
int kmeans_mode = MODE_ADPKERNEL;
int mapping_mode = MAPPING_MODE_CONST;

void getRGBSpatialStds(const Table2D<RGB> &img,double & std_r,double & std_g,double & std_b, const Table2D<bool> & ROI, bool countequal=true);

int main(int argc, char * argv[])
{
    clock_t start = clock();
	srand( (unsigned)time( NULL ) );
	if(1==argc){
        printf("Adaptive kernel K-means\n");
        printf("Usage:--------------------------------------------\n");
        printf("ppbcadp imagename w_smooth mapping theta1 theta2\n");
        printf("mapping can be const or identity or power or log\n");
        printf("theta1 and theta2 are parameters for mappings\n");
        printf("End of Usage--------------------------------------\n");
        exit(-1);
    }
	int init_mode = INIT_BOX;
	if(strcmp(argv[1],"all")==0){
	    // read directory
	    DIR *dpdf;
        struct dirent *epdf;
        dpdf = opendir("/home/mtang73/dataset/selected/smallcut/");
        if (dpdf == NULL) return -1;
        int imgid=0;
        double returnvs[N_IMG];
        double energies[N_IMG];
        while (epdf = readdir(dpdf)){
            char filename[10];
            strcpy(filename,epdf->d_name);
            if(strlen(filename)<=2) continue;
            char * shortname = (char *) malloc(strlen(filename)-3);
            strncpy(shortname,filename,strlen(filename)-4);
            shortname[strlen(filename)-4]='\0';
            printf("image %d : %s\n",imgid,shortname);
            // run for each image
            argv[1] = shortname;
	        kernelsegmentation(argc,argv,init_mode,returnvs+imgid,energies+imgid);
            free(shortname);
            imgid++;
        }
        for(int i=0;i<N_IMG;i++)
            printf("%.2f \n",returnvs[i]);
        printf("average error: %f\n", averagev(returnvs, N_IMG));
	}
	else
	    kernelsegmentation(argc,argv,init_mode);
	clock_t finish = clock();
	double time_to_now = (double)(finish-start)/CLOCKS_PER_SEC;
	cout<<"evaluation time "<<time_to_now<<" seconds!"<<endl;
	return -1;
}

void kernelsegmentation(int argc,char *argv[],int init_mode,double * returnv, double * energyv)
{
    // which kmeans method to run
	const char * mode_str="adpkernel";;
	kmeans_mode=MODE_ADPKERNEL;
    printf("adaptive kernel segmentation! argc =%d \n",argc);
    
    // set root directory and database directory
	char * destdir,	* sourcedir;
	sourcedir = "/home/mtang73/dataset/selected/";
	destdir = "/home/mtang73/kernel/output/adaptive/whard/";
	
	// set energy parameter
	double w_smooth = atof(argv[2]);
    if(0==strcmp(argv[3],"const"))  mapping_mode = MAPPING_MODE_CONST;
    else if(0==strcmp(argv[3],"identity"))  mapping_mode = MAPPING_MODE_IDENTITY;
    else if(0==strcmp(argv[3],"power"))  mapping_mode = MAPPING_MODE_POWER;
    else if(0==strcmp(argv[3],"truncation")) mapping_mode = MAPPING_MODE_TRUNCATION;
    else {printf("mode not valid!\n");exit(-1);}
    
	double theta1 = atof(argv[4]); // parameters for adaptive kernel
	double theta2 = atof(argv[5]);
	
	// read image
	char * imgname = argv[1];
	printf("argc = %d image %s\n",argc,imgname);
	Image image = Image((sourcedir+string("images/")+string(imgname) + string(".bmp")).c_str(), imgname, 8.0, 8);
    image.print();
    
    Table2D<RGB> rgbimg = image.img;
    int img_w = rgbimg.getWidth();
	int img_h = rgbimg.getHeight();
	
    //image.sigma_square = 1000000000;
	//image.computesmoothnesscost();
    
    // read subpixel images
    Table2D<vector<float> > vector_img(img_w,img_h,vector<float>(3));
	Table2D<double> column_img;
	readbinfile(column_img, (sourcedir+string("subpixelimages/")+string(imgname) + string(".bin")).c_str(), 1, img_h*img_w*3);
	int idx = 0;
	for(int c=0;c<3;c++){
		for(int i=0;i<img_w;i++){
			for(int j=0;j<img_h;j++){
				vector_img[i][j][c] = column_img[0][idx++];
			}
		}
	}

    // init labeling and hard constraints
	Table2D<Label> init_labeling;
	Table2D<Label> hard_constraints(img_w,img_h,NONE);
	
	if(INIT_BOX == init_mode){
		init_labeling = getinitlabeling(loadImage<RGB>((sourcedir+string("templates/")+string(imgname) + string(".bmp")).c_str()),0);
		for(int i=0;i<img_w;i++){
			for(int j=0;j<img_h;j++){
				if(init_labeling[i][j]==BKG) hard_constraints[i][j] = BKG;
				else hard_constraints[i][j] = NONE;
			}
		}
	}
	else if(INIT_RANDSEEDS==init_mode){
        cout<<"INIT_RANDSEEDS"<<endl;
        Table2D<RGB> seedimg = loadImage<RGB>((sourcedir+string("randomseeds/")+string(imgname) + string(".bmp")).c_str());
        init_labeling = getinitlabelingFB(seedimg, red, blue);
    }
    else if(INIT_STROKES==init_mode){
        cout<<"INIT_STROKES"<<endl;
        Table2D<RGB> seedimg = loadImage<RGB>((sourcedir+string("strokes/")+string(imgname) + string(".bmp")).c_str());
        init_labeling = getinitlabelingFB(seedimg, red, blue);
        outv(countintable(init_labeling,OBJ));
        outv(countintable(init_labeling,BKG));
        for(int i=0;i<img_w;i++){
			for(int j=0;j<img_h;j++){
				if(init_labeling[i][j]==BKG) hard_constraints[i][j] = BKG;
				else hard_constraints[i][j] = NONE;
			}
		}
    }

	// adaptive kernel using density mapping theory
	Table2D<double> d_k;
	readbinfile(d_k, (sourcedir+string("dk/")+string(imgname) + string(".bin")).c_str(), img_w, img_h);
	for(int i=0;i<img_w;i++){
        for(int j=0;j<img_h;j++){
	    	d_k[i][j] = max(d_k[i][j],0.2);
		}
	}
    // physical densities
    Table2D<double> densities=table_filter(vector_img,Table2D<bool>(img_w,img_h,true), d_k,"p");
	for(int i=0;i<img_w;i++){
        for(int j=0;j<img_h;j++){
			densities[i][j] = densities[i][j] /(3.1415926*4.0/3.0*d_k[i][j]*d_k[i][j]*d_k[i][j]);
		}
	}
	outv(densities.getMax());
	outv(densities.getMin());
	outv(densities.getMean());
	double mediandensity = tablemedian(densities);
	outv(mediandensity);
    
	Table2D<double> adpsigmas(img_w,img_h,0);
	double max_density = densities.getMax();
	//printf("density mapping: min(density,%.3f)\n",threshold);
	printf("density mapping: 200* (density/3)^(1.0/ %.3f)\n",gamma);
	printf("sigma: %.3f * pow(d'/d,0.333)\n",gamma);
    // density mapping
    double olddensity,newdensity;
    for(int i=0;i<img_w;i++){
        for(int j=0;j<img_h;j++){
            olddensity = densities[i][j];
            switch(mapping_mode){
                case MAPPING_MODE_CONST:
                    newdensity = theta1;
                    break;
                case MAPPING_MODE_TRUNCATION:
                    newdensity = min((theta2/max_density*theta1)*olddensity,theta2);
                    break;
                case MAPPING_MODE_POWER:
                    newdensity = theta1 * pow(olddensity / mediandensity, 1.0/theta2);
                    break;
                case MAPPING_MODE_IDENTITY:
                    newdensity = olddensity;
                    break;
                default:
                    break;
            }
            adpsigmas[i][j] = pow(newdensity/olddensity,0.333);
        }
    }
	outv(adpsigmas.getMax());
    outv(adpsigmas.getMin());
    outv(adpsigmas.getMean());
    
	// ground truth
	Table2D<int> gtimg = loadImage<RGB>((sourcedir+string("groundtruth/")+string(imgname) + string(".bmp")).c_str());
	Table2D<Label> gt = getinitlabeling(gtimg,255,0);
	if(strcmp(argv[argc-1],"gt")==0)
		init_labeling=gt;

    //hard_constraints.reset(NONE);
	PPBCadp ppbc = PPBCadp(image, 1.0, w_smooth, &vector_img,adpsigmas,hard_constraints);
	ppbc.setpara(-0.2,0.2,0.008,20,false);
	int bo_itr = 0,pbo_itr = 0;
	double ppbc_energy = 1e+10;
	while(++bo_itr){
		ppbc.setinitlabeling(init_labeling);
		BreakPoint best_bp = ppbc.explorepara(0);
		best_bp.original_e = ppbc.computeenergy(best_bp.solution);
		best_bp.print();
		//return;
		if((best_bp.original_e<ppbc_energy-0.1)&&(best_bp.ssize!=0)&&(best_bp.ssize!=img_w*img_h))
		{
			init_labeling = best_bp.solution;
			ppbc_energy = best_bp.original_e;
			outv(ppbc_energy);
			//savebinarylabeling(rgbimg,init_labeling,(destdir+string(imgname)+"_ppbc_"+pch::to_string(ppbc_itr)+string(".bmp")).c_str());
		}
		else
			break;
		if(bo_itr>20)
			break;
	}
	savebinarylabeling(rgbimg,init_labeling,(destdir+string(imgname)+string("_")+string(argv[3])+string("_")+string(argv[4])+string("_")+string(argv[5])+string("_s")+pch::to_string(w_smooth)+string("_bo_e")+pch::to_string((int)ppbc_energy)+string(".bmp")).c_str());
    
    ppbc.showflag = true;
    ppbc.computeenergy(init_labeling);
    ppbc.showflag = false;
    
	while(++pbo_itr){
		if(pbo_itr>10)
			break;
		ppbc.setinitlabeling(init_labeling);
		ppbc.explore();
		//ppbc.savesolutions(image,string(destdir));
		BreakPoint best_bp = ppbc.SelectBestBP();
		best_bp.print();
		
		//return;
		if((best_bp.original_e<ppbc_energy-0.1)&&(best_bp.ssize!=0)&&(best_bp.ssize!=img_w*img_h))
		{
			init_labeling = best_bp.solution;
			ppbc_energy = best_bp.original_e;
			outv(ppbc_energy);
			//savebinarylabeling(rgbimg,init_labeling,(destdir+string(imgname)+"_pbo_"+pch::to_string(pbo_itr)+string(".bmp")).c_str());
		}
		else
			break;
	}

	int solution_obj_size = countintable(init_labeling,OBJ);
	cout<<"foreground ratio"<<(double)solution_obj_size/img_w/img_h<<"ratio"<<endl;
	int box_size = countintable(hard_constraints,NONE);
	double errorrate = geterrorcount(init_labeling,gt)/(double)(box_size);
	cout<<"error"<<errorrate<<"error"<<endl;
	//cout<<"fmeasure"<<fmeasure(init_labeling, gt, OBJ)<<"fmeasure\n";
	//cout<<"jaccard"<<jaccard(init_labeling, gt, OBJ)<<"jaccard\n";
	if(returnv!=NULL) *returnv = errorrate;
	if(energyv!=NULL) *energyv = ppbc_energy;
	//return;
	savebinarylabeling(rgbimg,init_labeling,(destdir+string(imgname)+string("_")+string(argv[3])+string("_")+string(argv[4])+string("_")+string(argv[5])+string("_s")+pch::to_string(w_smooth)+string("_pbo_e")+pch::to_string((int)ppbc_energy)+string(".bmp")).c_str());
	
	ppbc.showflag = true;
    ppbc.computeenergy(init_labeling);
	return;
}

void getRGBSpatialStds(const Table2D<RGB> &img,double & std_r,double & std_g,double & std_b, const Table2D<bool> & ROI, bool countequal)
{
	int node_id = 0;
	int img_w = img.getWidth();
	int img_h = img.getHeight();
	double sigma_square_sum_r = 0;
	double sigma_square_sum_g = 0;
	double sigma_square_sum_b = 0;
	int sigma_square_count_r = 0;
	int sigma_square_count_g = 0;
	int sigma_square_count_b = 0;
	Point kernelshifts [] = {Point(1,0),Point(0,1)};
	for (int y=0; y<img_h; y++) // adding edges (n-links)
	{
		for (int x=0; x<img_w; x++) 
		{ 
			Point p(x,y);
			for(int i=0;i<2;i++)
			{
				Point q = p + kernelshifts[i];
				if(img.pointIn(q)&&ROI[x][y])
				{
					if((countequal)||((!countequal)&&(img[p].r!=img[q].r))){
						sigma_square_sum_r += pow(double(img[p].r)-double(img[q].r),2.0);
						sigma_square_count_r ++;
					}
					if((countequal)||((!countequal)&&(img[p].g!=img[q].g))){
						sigma_square_sum_g += pow(double(img[p].g)-double(img[q].g),2.0);
						sigma_square_count_g ++;
					}
					if((countequal)||((!countequal)&&(img[p].b!=img[q].b))){
						sigma_square_sum_b += pow(double(img[p].b)-double(img[q].b),2.0);
						sigma_square_count_b ++;
					}
				}
			}
		}
	}
	std_r =  sqrt(sigma_square_sum_r/sigma_square_count_r);
	std_g =  sqrt(sigma_square_sum_g/sigma_square_count_g);
	std_b =  sqrt(sigma_square_sum_b/sigma_square_count_b);
}
