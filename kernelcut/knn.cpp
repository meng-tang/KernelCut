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
#include "knn.h"
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
        dpdf = opendir("/home/mtang73/dataset/GrabCut/images/");
		//dpdf = opendir("/home/mtang73/dataset/interactive-gt/images/");
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
        printf("average error: %f\n", averagev(returnvs, N_IMG));
		for(int i=0;i<N_IMG;i++)
            printf("%.2f \n",returnvs[i]);
	}
	else{
	    double seed_e, seed_error;
	    knnsegmentation(argc,argv, INIT_BOX, &seed_error, &seed_e);
		printf("%.2f\t%.4f\n",seed_e,seed_error);
		
		/*printf("from box, GT, rand (50 times)\n");
		double box_e, box_error;
	    knnsegmentation(argc,argv, INIT_BOX, &box_error, &box_e);
		printf("%.2f\t%.4f\n",box_e,box_error);
		double gt_e, gt_error;
	    knnsegmentation(argc,argv, INIT_GT, &gt_error, &gt_e);
		printf("%.2f\t%.4f\n",gt_e,gt_error);
		double rand_e[50], rand_error[50];
		for(int i=0;i<50;i++){
			knnsegmentation(argc,argv, INIT_RANDBOX, rand_error+i, rand_e+i);
			printf("%.2f\t%.4f\n",rand_e[i],rand_error[i]);
		}*/
	}
	clock_t finish = clock();
	double time_to_now = (double)(finish-start)/CLOCKS_PER_SEC;
	cout<<"evaluation time "<<time_to_now<<" seconds!"<<endl;
	return -1;
}

void knnsegmentation(int argc,char *argv[],int init_mode, double * returnv,double * energyv)
{
    // which kmeans method to run
	const char * mode_str="knn";
    if(showflag) printf("knn segmentation! argc =%d \n",argc);
    
    // set root directory and database directory
	char * destdir,	* sourcedir;
	sourcedir = "/home/mtang73/dataset/GrabCut/";
	destdir = "/home/mtang73/kernel/output/aa/2009knn/";
	
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
    
	//char knnfile[100] = "/scratch/mtang73/GrabCutcroppedknn400/";
	//char knnfile[100] = "/scratch/mtang73/GrabCutcroppedknn1600pruned/";
	//char knnfile[100] = "/scratch/mtang73/interactive-gtknn400/";
	char knnfile[100] = "/scratch/mtang73/GrabCutcleanknn400/";
	strcat(knnfile,imgname);
	strcat(knnfile,".bin");
	if(showflag) printf("%s\n",knnfile);

	Table2D<int> knntable; 
	readbinfile(knntable,knnfile, KNN_K,img_w * img_h);
	//perturbconnection(knntable,img_w,img_h);

	// ground truth
	Table2D<int> gtimg = loadImage<RGB>((sourcedir+string("groundtruth/")+string(imgname) + string(".bmp")).c_str());
	Table2D<Label> gt = getinitlabeling(gtimg,255,0);
	double gt_obj_ratio = countintable(gt,OBJ)/(double)(img_w*img_h);
	if( (strcmp(argv[argc-1],"gt")==0) || (INIT_GT == init_mode) )
		init_labeling=gt;
	double current_flow = 1e+10;
	Table2D<double> colorweights = Table2D<double>(img_w,img_h,1.0);
	
	int imgid=0;
	while(1)
	{
		GraphType * g;
		g = new GraphType(/*estimated # of nodes*/ img_w*img_h, /*estimated # of edges*/ 4*img_w*img_h); 
	    g->add_node(img_w*img_h);    // adding nodes
		
		if(w_smooth>1e-10)
			addsmoothnessterm(g,image,w_smooth,Table2D<bool>(img_w,img_h,true),false);
		Table2D<double> capsource = Table2D<double>(img_w,img_h,0), capsink = Table2D<double>(img_w,img_h,0);

		Table2D<double> w_data_OBJ = zeroonekernel(knntable, getROI(init_labeling,OBJ),KNN_K);
		double obj_size  = countintable(init_labeling, OBJ);
		double obj_sum = w_data_OBJ.sum(getROI(init_labeling, OBJ));
		for(int i=0;i<img_w;i++){
			for(int j=0;j<img_h;j++){
				capsink[i][j] = (- 2*w_data_OBJ[i][j]/obj_size + obj_sum / obj_size / obj_size)*5 ;
			}
		}

		Table2D<double> w_data_BKG = zeroonekernel(knntable, getROI(init_labeling,BKG),KNN_K);
		double bkg_size = countintable(init_labeling, BKG);
		double bkg_sum = w_data_BKG.sum(getROI(init_labeling, BKG));
		for(int i=0;i<img_w;i++){
			for(int j=0;j<img_h;j++){
				capsource[i][j] = (- 2*w_data_BKG[i][j]/bkg_size + bkg_sum / bkg_size / bkg_size)*5;
			}
		}
		
		Table2D<double> diff = capsource - bkg_sum / bkg_size / bkg_size *5 - ( capsink - obj_sum / obj_size / obj_size *5 );
		outv(diff.getMax());
		outv(diff.getMin());
		outv(diff.getMean());
		savetableascolorimage(diff*(1600), (destdir+string(imgname) + string("_id") + pch::to_string(imgid) + string("_probs.bmp")).c_str());
		Table2D<double> diff2 = capsource - capsink;
		outv(diff2.getMax());
		outv(diff2.getMin());
		outv(diff2.getMean());
		savetableascolorimage(diff2*(1600), (destdir+string(imgname) + string("_id") + pch::to_string(imgid) + string("_unarys.bmp")).c_str());
		//return;
		
		double obj_ratio = countintable(init_labeling, OBJ)/ (double)(img_w*img_h);
		for(int i=0;i<img_w;i++){
			for(int j=0;j<img_h;j++){
				if(OBJ==hard_constraints[i][j])
					capsource[i][j] = 1e20;
				else if(BKG==hard_constraints[i][j])
					capsink[i][j] = 1e20;
				g->add_tweights(i+j*img_w,capsource[i][j],capsink[i][j]);
			}
		}

		double flow = g -> maxflow();
		if(showflag) outv(flow);
		
		Table2D<Label> m_labeling(img_w,img_h);
		if(!getlabeling(g,m_labeling))
		{
			cout<<"trivial solution!"<<endl;
			break; // trivial solution
		}
		delete g;

		if(init_labeling==m_labeling)
		{
			if(showflag) cout<<"labeling converged!"<<endl;
			break;
		}
		if(current_flow-flow<0.1)
		{
			if(showflag) cout<<"energy increase!"<<endl;
			break;
		}
		init_labeling = m_labeling;
		current_flow = flow;
		imgid++;
		savebinarylabeling(rgbimg,init_labeling,(destdir+string(imgname)+"_"+string(mode_str)+"_"+pch::to_string(imgid)+string(".bmp")).c_str());
		if(imgid>20){
			if(showflag) printf("max iterations.\n");
			break;
		}
	}

	int solution_obj_size = countintable(init_labeling,OBJ);
	if(showflag) cout<<"foreground ratio"<<(double)solution_obj_size/img_w/img_h<<"ratio"<<endl;
	int box_size = countintable(hard_constraints,NONE);
	double errorrate = geterrorcount(init_labeling,gt)/(double)(box_size);
	double myjaccard = jaccard(init_labeling, gt, OBJ);
	if(showflag) cout<<"error"<<errorrate<<"error"<<endl;
	if(showflag) cout<<"fmeasure"<<fmeasure(init_labeling, gt, OBJ)<<"fmeasure\n";
	if(showflag) cout<<"jaccard"<<myjaccard<<"jaccard\n";
	
	savebinarylabeling(rgbimg, init_labeling,(destdir+string(imgname) +"_"+string(mode_str)+"_k"+pch::to_string(KNN_K)+"_s"+pch::to_string(w_smooth)+string(".bmp")+"_e"+pch::to_string(errorrate)).c_str());
	
	if(INIT_RANDBOX==init_mode){
	    Table2D<Label> fliplabeling = init_labeling;
	    for(int i=0;i<img_w;i++){
		    for(int j=0;j<img_h;j++){
			    if(init_labeling[i][j]==OBJ)
				    fliplabeling[i][j] = BKG;
			    else
				    fliplabeling[i][j] = OBJ;
		    }
	    }
	    double errorrate2 = geterrorcount(fliplabeling,gt)/(double)(img_w*img_h);
	    if(errorrate2<errorrate)
		    errorrate = errorrate2;
    }
	if(returnv!=NULL) *returnv = errorrate;
	if(energyv!=NULL) *energyv = current_flow;
}

