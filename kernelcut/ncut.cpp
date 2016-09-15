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
#define MODE_NCUT 6
#define INIT_BOX 1
#define INIT_RANDBOX 2
#define INIT_GT 3
#define INIT_RANDSEEDS 4

using namespace std;

void getRGBSpatialStds(Table2D<vector<float> > img,double & std_r,double & std_g,double & std_b);
void ncutsegmentation(int argc,char *argv[], int init_mode,double * returnv=NULL,double * energyv=NULL);

// global variables
int kmeans_mode = MODE_NCUT;
bool showflag = true;

#define N_IMG 50

int main(int argc, char * argv[])
{
    clock_t start = clock();
	srand( (unsigned)time( NULL ) );
	if(argc==1){
		printf("Usage: -------------------------------\n");
		printf("Ncut with KNN: imgname KNN_K w_smooth\n");
		printf("Ncut with Gaussian (RGB): imgname sigmaR sigmaG sigmaB w_smooth\n");
		printf("Ncut with Gaussian (RGBXY): imgname sigmaR sigmaG sigmaB sigmaX sigmaY w_smooth\n");
		printf("End of usage: ------------------------\n");
		return -1;
	}
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
	        ncutsegmentation(argc,argv,init_mode,returnvs+imgid,energies+imgid);
            free(shortname);
            imgid++;
        }
        printf("average error: %f\n", averagev(returnvs, N_IMG));
		for(int i=0;i<N_IMG;i++)
            printf("%.3f \n",returnvs[i]);
	}
	else
	    ncutsegmentation(argc,argv,init_mode);
	clock_t finish = clock();
	double time_to_now = (double)(finish-start)/CLOCKS_PER_SEC;
	cout<<"evaluation time "<<time_to_now<<" seconds!"<<endl;
	return -1;
}

void ncutsegmentation(int argc,char *argv[],int init_mode, double * returnv,double * energyv)
{
    // which kmeans method to run
	const char * mode_str="ncut";
    printf("knn segmentation! argc =%d \n",argc);
    char * graphmode; // gaussianrgb,gaussianrgbxy,knn

    // set root directory and database directory
	char * destdir,	* sourcedir, * LABdir;
	sourcedir = "/home/mtang73/dataset/GrabCutclean/";
	//sourcedir = "/home/mtang73/dataset/interactive-gt/";
	destdir = "/home/mtang73/kernel/output/ncut/2009gaussian/";
	//destdir = "/home/mtang73/kernel/output/ncut/interactive-gt/";
	
	// set energy parameter
	double w_smooth;
	double binsize = 8.0; // for histogram or smoothed histogram
	double sigmas[3];double sigmax;double sigmay;
	int KNN_K;
	switch(argc){
		case 4:
			graphmode = "knn";
			KNN_K = atoi(argv[2]);
			w_smooth = atof(argv[3]);
			break;
		case 6:
			graphmode = "gaussianrgb";
			sigmas[0]=atof(argv[2]);sigmas[1]=atof(argv[3]);sigmas[2]=atof(argv[4]);
			w_smooth = atof(argv[5]);
			break;
		case 8:
			graphmode = "gaussianrgbxy";
			sigmas[0]=atof(argv[2]);sigmas[1]=atof(argv[3]);sigmas[2]=atof(argv[4]);
			sigmax=atof(argv[5]);sigmay=atof(argv[6]);
			w_smooth = atof(argv[7]);
			break;
		default:
			break;
	}
	// read image
	char * imgname = argv[1];
	printf("argc = %d image %s\n",argc,imgname);
    Table2D<RGB> rgbimg = loadImage<RGB>((sourcedir+string("images/")+string(imgname) + string(".bmp")).c_str());
    cout<<rgbimg.getWidth()<<' '<<rgbimg.getHeight()<<endl;
    int img_w = rgbimg.getWidth();
	int img_h = rgbimg.getHeight();
	Image image = Image((sourcedir+string("images/")+string(imgname) + string(".bmp")).c_str(),imgname,8,8); // rgb image
	
	image.sigma_square = 1000000000;
	image.computesmoothnesscost();

	Table2D<vector<float> > vector_img(img_w,img_h,vector<float>(3));
	// read subpixel images
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
    
    if(0==strcmp(graphmode,"gaussianrgb")){
        getRGBSpatialStds(vector_img,sigmas[0],sigmas[1],sigmas[2]);
        outv(sigmas[0]);
        outv(sigmas[1]);
        outv(sigmas[2]);
        sigmas[0] = sigmas[0]*atof(argv[2]);
        sigmas[1] = sigmas[1]*atof(argv[3]);
        sigmas[2] = sigmas[2]*atof(argv[4]);
        outv(sigmas[0]);
        outv(sigmas[1]);
        outv(sigmas[2]);
    }
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
        cout<<"INIT_RANDSEEDS"<<endl;
        Table2D<RGB> seedimg = loadImage<RGB>((sourcedir+string("randomseeds/")+string(imgname) + string(".bmp")).c_str());
        init_labeling = getinitlabelingFB(seedimg, red, blue);
    }

	Table2D<int> knntable;
	Table2D<double> weights;
	if(strcmp(graphmode,"knn")==0){
		//char knnfile[100] = "/scratch/mtang73/GrabCutcroppedknn1600pruned/";
		//char knnfile[100] = "/scratch/mtang73/GrabCutknn400/";
		//char knnfile[100] = "/scratch/mtang73/GrabCutcleanknn400/";
		char knnfile[100] = "/scratch/mtang73/interactive-gtknn400/";
		strcat(knnfile,imgname);
		strcat(knnfile,".bin");
		printf("%s\n",knnfile);
		readbinfile(knntable,knnfile, KNN_K,img_w * img_h);
		//perturbconnection(knntable,img_w,img_h);
		weights = zeroonekernel(knntable, Table2D<bool>(img_w,img_h,true), KNN_K);
	}
	else if(strcmp(graphmode,"gaussianrgb")==0){
		weights = filterrgb(vector_img, Table2D<bool>(img_w,img_h,true), sigmas);
		outv(weights.getMean());
		outv(weights.getMax());
		outv(weights.getMin());
	}
	else if(strcmp(graphmode,"gaussianrgbxy")==0){
		weights = filterrgbxy(vector_img, Table2D<bool>(img_w,img_h,true), sigmas,sigmax,sigmay);
		/*int x = rand()%img_w;
		int y = rand()%img_h;
		outv(weights[x][y]);
		vector<float> I_p = vector_img[x][y];
		float wp = 0;
		for(int i=0;i<img_w;i++){
        	for(int j=0;j<img_h;j++){
				vector<float> I_q = vector_img[i][j];
				float diff = pow(I_p[0]-I_q[0],2) + pow(I_p[1]-I_q[1],2) + pow(I_p[2]-I_q[2],2);
				float diffxy = pow(x-i,2) + pow(y-j,2);
				wp += exp(-0.5*diff/5.0/5.0-0.5*diffxy/sigmax/sigmax); 
			}
		}
		outv(wp);
		outv(weights[x][y]/wp);	*/
	}
	// ground truth
	Table2D<int> gtimg = loadImage<RGB>((sourcedir+string("groundtruth/")+string(imgname) + string(".bmp")).c_str());
	Table2D<Label> gt = getinitlabeling(gtimg,255,0);
	if( (strcmp(argv[argc-1],"gt")==0) || (INIT_GT == init_mode) )
		init_labeling=gt;
	double current_flow = 1e+10;

	double offset = -2;
	for(int i=0;i<img_w;i++){
		for(int j=0;j<img_h;j++){
			offset += 1 / weights[i][j];
		}
	}
	offset = offset * 10000;
	int imgid=0;
	while(1)
	{
		GraphType * g;
		g = new GraphType(/*estimated # of nodes*/ img_w*img_h, /*estimated # of edges*/ 4*img_w*img_h); 
	    g->add_node(img_w*img_h);    // adding nodes
		
		if(w_smooth>1e-10)
			addsmoothnessterm(g,image,w_smooth,Table2D<bool>(img_w,img_h,true),false);
		Table2D<double> capsource = Table2D<double>(img_w,img_h,0), capsink = Table2D<double>(img_w,img_h,0);

		Table2D<double> w_data_OBJ;
		if(strcmp(graphmode,"knn")==0){
			w_data_OBJ = zeroonekernel(knntable, getROI(init_labeling,OBJ),KNN_K);
		}
		else if(strcmp(graphmode,"gaussianrgb")==0){
			w_data_OBJ = filterrgb(vector_img, getROI(init_labeling, OBJ), sigmas);
		}
		else if(strcmp(graphmode,"gaussianrgbxy")==0){
			w_data_OBJ = filterrgbxy(vector_img, getROI(init_labeling, OBJ), sigmas,sigmax,sigmay);
		}
		double obj_size  = weights.sum(getROI(init_labeling,OBJ));
		double obj_sum = w_data_OBJ.sum(getROI(init_labeling, OBJ));
		for(int i=0;i<img_w;i++){
			for(int j=0;j<img_h;j++){
				double wp = weights[i][j];
				capsink[i][j] = (1.0 / wp - 2*w_data_OBJ[i][j]/obj_size + obj_sum / obj_size / obj_size * wp)*10000 ;
			}
		}

		Table2D<double> w_data_BKG;
		if(strcmp(graphmode,"knn")==0){
			w_data_BKG = zeroonekernel(knntable, getROI(init_labeling,BKG),KNN_K);
		}
		else if(strcmp(graphmode,"gaussianrgb")==0){
			w_data_BKG = filterrgb(vector_img, getROI(init_labeling, BKG), sigmas);
		}
		else if(strcmp(graphmode,"gaussianrgbxy")==0){
			w_data_BKG = filterrgbxy(vector_img, getROI(init_labeling, BKG), sigmas,sigmax,sigmay);
		}
		double bkg_size = weights.sum(getROI(init_labeling,BKG));
		double bkg_sum = w_data_BKG.sum(getROI(init_labeling, BKG));
		for(int i=0;i<img_w;i++){
			for(int j=0;j<img_h;j++){
				double wp = weights[i][j];
				capsource[i][j] = (1.0 / wp - 2*w_data_BKG[i][j]/bkg_size + bkg_sum / bkg_size / bkg_size * wp)*10000;
			}
		}
		
		for(int i=0;i<img_w;i++){
			for(int j=0;j<img_h;j++){
				if(OBJ==hard_constraints[i][j])
					capsource[i][j] = 1e20;
				else if(BKG==hard_constraints[i][j])
					capsink[i][j] = 1e20;
				g->add_tweights(i+j*img_w,capsource[i][j],capsink[i][j]);
			}
		}
		// weighted kernel kmeans energy
		double wkk_e = capsink.sum(getROI(init_labeling,OBJ)) + capsource.sum(getROI(init_labeling,BKG));
		// noralized cut energy
		double ncut_e = ( w_data_OBJ.sum(getROI(init_labeling, BKG))/ weights.sum(getROI(init_labeling, BKG)) 
				+ w_data_OBJ.sum(getROI(init_labeling, BKG))/ weights.sum(getROI(init_labeling, OBJ)) ) * 10000;
		printf("ncut %.4f (wkk_e - offset %.4f)\n",ncut_e,wkk_e - offset);
		double flow = g -> maxflow();
		outv(flow);
		
		Table2D<Label> m_labeling(img_w,img_h);
		if(!getlabeling(g,m_labeling))
		{
			cout<<"trivial solution!"<<endl;
			break; // trivial solution
		}
		delete g;

		if(init_labeling==m_labeling)
		{
			cout<<"labeling converged!"<<endl;
			break;
		}
		if(current_flow-flow<0.1)
		{
			cout<<"energy increase!"<<endl;
			break;
		}
		init_labeling = m_labeling;
		current_flow = flow;
		imgid++;
		//savebinarylabeling(rgbimg,init_labeling,(destdir+string(imgname)+"_"+string(mode_str)+"_"+pch::to_string(imgid)+string(".bmp")).c_str());
		if(imgid>25){
			printf("max iterations.\n");
			break;
		}
	}

	int solution_obj_size = countintable(init_labeling,OBJ);
	cout<<"foreground ratio"<<(double)solution_obj_size/img_w/img_h<<"ratio"<<endl;
	int box_size = countintable(hard_constraints,NONE);
	double errorrate = geterrorcount(init_labeling,gt)/(double)(box_size);
	double myjaccard = jaccard(init_labeling, gt, OBJ);
	cout<<"error"<<errorrate<<"error"<<endl;
	cout<<"fmeasure"<<fmeasure(init_labeling, gt, OBJ)<<"fmeasure\n";
	cout<<"jaccard"<<myjaccard<<"jaccard\n";
	if(returnv!=NULL) *returnv = errorrate;
	if(energyv!=NULL) *energyv = current_flow;
	
	return;
	if(strcmp(graphmode,"knn")==0){
		savebinarylabeling(rgbimg, init_labeling,(destdir+string(imgname) +"_"+string(mode_str)+"_k"+pch::to_string(KNN_K)+"_s"+pch::to_string(w_smooth)+string(".bmp")).c_str());
	}
	else if(strcmp(graphmode,"gaussianrgb")==0){
		savebinarylabeling(rgbimg, init_labeling,(destdir+string(imgname) +"_grgb"+"_w"+pch::to_string(atof(argv[2]))+"_s"+pch::to_string(w_smooth)+string(".bmp")).c_str());
	}
	else if(strcmp(graphmode,"gaussianrgbxy")==0){
		savebinarylabeling(rgbimg, init_labeling,(destdir+string(imgname) +"_grgbxy"+"_s"+pch::to_string(w_smooth)+string(".bmp")).c_str());
	}
	
	return;
}

void getRGBSpatialStds(Table2D<vector<float> > img,double & std_r,double & std_g,double & std_b)
{
	int img_w = img.getWidth();
	int img_h = img.getHeight();
	double sigma_square_sum_r = 0;
	double sigma_square_sum_g = 0;
	double sigma_square_sum_b = 0;
	int sigma_square_count_r = 0;
	int sigma_square_count_g = 0;
	int sigma_square_count_b = 0;
	Point kernelshifts [] = {Point(1,0),Point(0,1)};
	for (int y=0; y<img_h; y++)
	{
		for (int x=0; x<img_w; x++) 
		{ 
			Point p(x,y);
			for(int i=0;i<2;i++)
			{
				Point q = p + kernelshifts[i];
				if(img.pointIn(q))
				{
					sigma_square_sum_r += pow(float(img[p][0])-float(img[q][0]),2.0);
					sigma_square_count_r ++;
					sigma_square_sum_g += pow(float(img[p][1])-float(img[q][1]),2.0);
					sigma_square_count_g ++;
					sigma_square_sum_b += pow(float(img[p][2])-float(img[q][2]),2.0);
					sigma_square_count_b ++;
				}
			}
		}
	}
	std_r =  sqrt(sigma_square_sum_r/sigma_square_count_r);
	std_g =  sqrt(sigma_square_sum_g/sigma_square_count_g);
	std_b =  sqrt(sigma_square_sum_b/sigma_square_count_b);
}
