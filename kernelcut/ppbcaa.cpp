#include <stdio.h>
#include <iostream>
#include <ctime>
#include <string>

#include <graph.h>
#include <EasyBMP.h>
#include "dirent.h"
#include "basicutil.h"
#include "knn.h"
#include "PPBCaa.h"
#define MODE_KNN 5
#define INIT_BOX 1
#define INIT_RANDBOX 2
#define INIT_RANDSEEDS 3
#define INIT_GT 4

using namespace std;
void knnsegmentation(int argc,char *argv[], int init_mode, double * returnv=NULL,double * energyv=NULL);

// global variables
int kmeans_mode = MODE_KNN;
bool showflag = true;
#define N_IMG 50
int main(int argc, char * argv[])
{
    clock_t start = clock();
	srand( (unsigned)time( NULL ) );
	int init_mode = INIT_BOX;
	if(strcmp(argv[1],"all")==0){
	    // read directory
	    DIR *dpdf;
        struct dirent *epdf;
        dpdf = opendir("/home/mtang73/dataset/GrabCutclean/images/");
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
	        knnsegmentation(argc,argv,init_mode, returnvs+imgid,energies+imgid);
            free(shortname);
            imgid++;
        }
		printf("energies:\n");
		for(int i=0;i<N_IMG;i++)
            printf("%.1f \n",energies[i]);
		printf("errors:\n");
		for(int i=0;i<N_IMG;i++)
            printf("%.3f \n",returnvs[i]);
		printf("average error: %f\n", averagev(returnvs, N_IMG));
	}
	else{
	    double seed_e, seed_error;
	    knnsegmentation(argc,argv, INIT_BOX, &seed_error, &seed_e);
		printf("%.2f\t%.4f\n",seed_e,seed_error);
	}
	clock_t finish = clock();
	double time_to_now = (double)(finish-start)/CLOCKS_PER_SEC;
	printf("evaluation time %.5f seconds\n",time_to_now);
	return -1;
}

void knnsegmentation(int argc,char *argv[],int init_mode, double * returnv,double * energyv)
{
    // which kmeans method to run
	const char * mode_str="ppbc";
    if(showflag) printf("knn segmentation! argc =%d \n",argc);
    
    // set root directory and database directory
	char * destdir,	* sourcedir;
	sourcedir = "/home/mtang73/dataset/GrabCutclean/";
	destdir = "/home/mtang73/kernel/output/ppbc/aa/";
	
	// set energy parameter
	double w_smooth;
	double binsize = 8.0; // for histogram or smoothed histogram
	
	// read image
	char * imgname = argv[1];
	if(showflag) printf("argc = %d image %s\n",argc,imgname);
    Table2D<RGB> rgbimg = loadImage<RGB>((sourcedir+string("images/")+string(imgname) + string(".bmp")).c_str());
    if(showflag) cout<<rgbimg.getWidth()<<' '<<rgbimg.getHeight()<<endl;
    int img_w = rgbimg.getWidth();
	int img_h = rgbimg.getHeight();
	Image image = Image((sourcedir+string("images/")+string(imgname) + string(".bmp")).c_str(),imgname,8,8); // rgb image
	//image.sigma_square = 1000000000;
	//image.computesmoothnesscost();
	
	int KNN_K = atoi(argv[2]);
	w_smooth = atof(argv[3]);
    // init labeling and hard constraints
	Table2D<Label> init_labeling;
	Table2D<Label> hard_constraints(img_w,img_h,NONE);
	
	if(INIT_BOX == init_mode){
		init_labeling = getinitlabeling(loadImage<RGB>((sourcedir+string("templates/")+string(imgname) + string("_t.bmp")).c_str()),0);
		for(int i=0;i<img_w;i++){
			for(int j=0;j<img_h;j++){
				if(init_labeling[i][j]==BKG) hard_constraints[i][j] = BKG;
				else hard_constraints[i][j] = NONE;
			}
		}
	}
	else if(INIT_RANDBOX==init_mode){
		if(showflag) cout<<"INIT_RANDBOX"<<endl;
		init_labeling = Table2D<Label>(img_w,img_h,BKG);
		double wmin = (rand()%10000/10000.0)*0.8; // 0-0.8
		double wmax = wmin + (1-wmin) * (1 - (rand()%10000/10000.0)*0.8 );
		double hmin = (rand()%10000/10000.0)*0.8;
		double hmax = hmin + (1-hmin) * (1 - (rand()%10000/10000.0)*0.8 );
		//wmin = (rand()%1000/1000.0)*0.5;wmax = wmin+(rand()%1000/1000.0/2+0.5)*(1.0-wmin);
		//hmin = (rand()%1000/1000.0)*0.5;hmax = hmin+(rand()%1000/1000.0/2+0.5)*(1.0-hmin);
		for(int i=img_w*wmin;i<img_w*wmax;i++){
			for(int j=img_h*hmin;j<img_h*hmax;j++){
					init_labeling[i][j] = OBJ;
			}
		}
	}
	else if(INIT_RANDSEEDS==init_mode){
        if(showflag) cout<<"INIT_RANDSEEDS"<<endl;
        Table2D<RGB> seedimg = loadImage<RGB>((sourcedir+string("randomseeds/")+string(imgname) + string(".bmp")).c_str());
        init_labeling = getinitlabelingFB(seedimg, red, blue);
    }
    
	//char knnfile[100] = "/scratch/mtang73/GrabCutcleanknn400/";
	char knnfile[100] = "/scratch/mtang73/Grabcutcleanknn3200/";
	strcat(knnfile,imgname);
	strcat(knnfile,".bin");
	if(showflag) printf("%s\n",knnfile);

	// read knn table
	Table2D<int> temp_knntable; 
	readbinfile(temp_knntable,knnfile, KNN_K,img_w * img_h);
	Table2D<int> knntable(img_w*img_h,KNN_K); // index from zero, KNN_K rows
	for(int i=0;i<KNN_K;i++){
		for(int j=0;j<img_w*img_h;j++){
			knntable[j][i] = temp_knntable[i][j]-1;
		}
	}
	temp_knntable.resize(1,1);

	// full knn
	readbinfile(temp_knntable,knnfile, 400,img_w * img_h);
	Table2D<int> knntable2(img_w*img_h,400); // index from zero, KNN_K rows
	for(int i=0;i<400;i++){
		for(int j=0;j<img_w*img_h;j++){
			knntable2[j][i] = temp_knntable[i][j]-1;
		}
	}
	temp_knntable.resize(1,1);
	PPBCaa ppbc2 = PPBCaa(image, 5.0, 0.96, &knntable2, 400,hard_constraints);
	ppbc2.setpara(-0.1,0.1,0.002,20,false);


	//perturbconnection(knntable,img_w,img_h);

	// ground truth
	Table2D<int> gtimg = loadImage<RGB>((sourcedir+string("groundtruth/")+string(imgname) + string(".bmp")).c_str());
	Table2D<Label> gt = getinitlabeling(gtimg,255,0);
	double gt_obj_ratio = countintable(gt,OBJ)/(double)(img_w*img_h);
	if( (strcmp(argv[argc-1],"gt")==0) || (INIT_GT == init_mode) )
		init_labeling=gt;

	PPBCaa ppbc = PPBCaa(image, 5.0, w_smooth, &knntable, KNN_K,hard_constraints);
	ppbc.setpara(-0.1,0.1,0.002,20,false);

	int bo_itr = 0,pbo_itr=0;double ppbc_energy = 1e+10;
	while(++bo_itr){
		ppbc.setinitlabeling(init_labeling);
		BreakPoint best_bp = ppbc.explorepara(0);
		//best_bp.print();
		best_bp.original_e = ppbc.computeenergy(best_bp.solution);
		if((best_bp.original_e<ppbc_energy-0.1)&&(best_bp.ssize!=0)&&(best_bp.ssize!=img_w*img_h))
		{
			init_labeling = best_bp.solution;
			ppbc_energy = best_bp.original_e;
			outv(ppbc_energy);
			//savebinarylabeling(rgbimg,init_labeling,(destdir+string(imgname)+"_ppbc_"+pch::to_string(ppbc_itr)+string(".bmp")).c_str());
		}
		else
			break;
		if(bo_itr>20) break;
	}
	savebinarylabeling(rgbimg,init_labeling,(destdir+string(imgname)+string("_bo_k")+pch::to_string(8*KNN_K)+string("_e")+pch::to_string((int)ppbc_energy)+string(".bmp")).c_str());
	while(++pbo_itr){
		ppbc.setinitlabeling(init_labeling);
		ppbc.explore();
		BreakPoint best_bp = ppbc.SelectBestBP();
		//best_bp.print();
		if((best_bp.original_e<ppbc_energy-0.1)&&(best_bp.ssize!=0)&&(best_bp.ssize!=img_w*img_h))
		{
			init_labeling = best_bp.solution;
			ppbc_energy = best_bp.original_e;
			outv(ppbc_energy);
			//savebinarylabeling(rgbimg,init_labeling,(destdir+string(imgname)+"_ppbc_"+pch::to_string(ppbc_itr)+string(".bmp")).c_str());
		}
		else
			break;
		if(pbo_itr>8) break;
	}
	
	int solution_obj_size = countintable(init_labeling,OBJ);
	if(showflag) cout<<"foreground ratio"<<(double)solution_obj_size/img_w/img_h<<"ratio"<<endl;
	int box_size = countintable(hard_constraints,NONE);
	double errorrate = geterrorcount(init_labeling,gt)/(double)(box_size);
	double myjaccard = jaccard(init_labeling, gt, OBJ);
	if(showflag) cout<<"error"<<errorrate<<"error"<<endl;
	if(showflag) cout<<"fmeasure"<<fmeasure(init_labeling, gt, OBJ)<<"fmeasure\n";
	if(showflag) cout<<"jaccard"<<myjaccard<<"jaccard\n";
	savebinarylabeling(rgbimg,init_labeling,(destdir+string(imgname)+string("_pbo_k")+pch::to_string(8*KNN_K)+string("_e")+pch::to_string((int)ppbc_energy)+string(".bmp")).c_str());
	//savebinarylabeling(rgbimg, init_labeling,(destdir+string(imgname) +"_"+string(mode_str)+"_k"+pch::to_string(KNN_K)+"_s"+pch::to_string(w_smooth)+string(".bmp")).c_str());
	outv(ppbc2.computeenergy(init_labeling));
	if(returnv!=NULL) *returnv = errorrate;
	if(energyv!=NULL) *energyv = ppbc_energy;
}

