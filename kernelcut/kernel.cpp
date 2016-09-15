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

#define MODE_KERNEL 3
#define MODE_KERNELXY 5

#define INIT_RANDSEEDS 3
#define INIT_BOX 1
#define INIT_STROKES 5

#define N_IMG 15

using namespace std;
void kernelsegmentation(int argc,char *argv[], int init_mode,double * returnv=NULL,double * energyv=NULL);

// global variables
int kmeans_mode = MODE_KERNEL;
void getRGBSpatialStds(const Table2D<RGB> &img,double & std_r,double & std_g,double & std_b, const Table2D<bool> & ROI, bool countequal=true);

int main(int argc, char * argv[])
{
    clock_t start = clock();
	srand( (unsigned)time( NULL ) );
	int init_mode = INIT_BOX;
	if(strcmp(argv[1],"all")==0){
	    // read directory
	    DIR *dpdf;
        struct dirent *epdf;
        dpdf = opendir("/home/mtang73/dataset/selected/images/");
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
            printf("%.2f \n",energies[i]);
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
	const char * mode_str;
	if(strcmp("kernel",argv[2])==0)
	    {kmeans_mode=MODE_KERNEL;mode_str="RGB";}
	else if(strcmp("kernelxy",argv[2])==0)
	    {kmeans_mode=MODE_KERNELXY;mode_str="RGBXY";}
	else{
	    printf("please specify mode!\n");
	    return;
	}
    printf("kernel segmentation! argc =%d \n",argc);
    
    // set root directory and database directory
	char * destdir,	* sourcedir;
	//sourcedir = "/home/mtang73/dataset/GrabCutcropped/";
	//LABdir = "/home/mtang73/dataset/GrabCutcropped/images/";
	//sourcedir = "/home/mtang73/dataset/cooperativecut/";
	sourcedir = "/home/mtang73/dataset/selected/";
	destdir = "/home/mtang73/kernel/output/aa/selected/";
	
	// set energy parameter
	double w_smooth;
	double sigmas[3]; // for kernel in color space
	double sigmax, sigmay; // for kernel in XY domain
	if(MODE_KERNEL==kmeans_mode)
	{
	    sigmas[0] = atof(argv[3]);sigmas[1] = atof(argv[4]);sigmas[2] = atof(argv[5]);
	    w_smooth = atof(argv[6]);
		//getRGBSpatialStds(image.img,sigmas[0],sigmas[1],sigmas[2], Table2D<bool>(img_w,img_h,true));
    }
	else if(MODE_KERNELXY==kmeans_mode)
	{
	    sigmas[0] = atof(argv[3]);sigmas[1] = atof(argv[4]);sigmas[2] = atof(argv[5]);
		sigmax = atof(argv[6]);sigmay = atof(argv[7]);
	    w_smooth = atof(argv[8]);
	}
	
	// read image
	char * imgname = argv[1];
	printf("argc = %d image %s\n",argc,imgname);
    Table2D<RGB> rgbimg = loadImage<RGB>((sourcedir+string("images/")+string(imgname) + string(".bmp")).c_str());
    cout<<rgbimg.getWidth()<<' '<<rgbimg.getHeight()<<endl;
    int img_w = rgbimg.getWidth();
	int img_h = rgbimg.getHeight();
	
    Image image = Image((sourcedir+string("images/")+string(imgname) + string(".bmp")).c_str(), imgname, 8, 8);
    image.print();
    //image.sigma_square = 1000000000;
	//image.computesmoothnesscost();
    
    Table2D<vector<float> > vector_img(img_w,img_h,vector<float>(3));
    for(int i=0;i<img_w;i++){
        for(int j=0;j<img_h;j++){
            vector_img[i][j][0] = float(image.img[i][j].r);
            vector_img[i][j][1] = float(image.img[i][j].g);
            vector_img[i][j][2] = float(image.img[i][j].b);
        }
    }
    
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
        Table2D<RGB> seedimg = loadImage<RGB>((sourcedir+string("seeds/")+string(imgname) + string(".bmp")).c_str());
        init_labeling = getinitlabelingFB(seedimg, red, blue);
        for(int i=0;i<img_w;i++){
			for(int j=0;j<img_h;j++){
				if(OBJ == init_labeling[i][j]) hard_constraints[i][j] = OBJ;
				else if(BKG == init_labeling[i][j]) hard_constraints[i][j] = BKG;
			}
		}
    }
    
	// ground truth
	Table2D<int> gtimg = loadImage<RGB>((sourcedir+string("groundtruth/")+string(imgname) + string(".bmp")).c_str());
	Table2D<Label> gt = getinitlabeling(gtimg,255,0);
	if(strcmp(argv[argc-1],"gt")==0)
		init_labeling=gt;
		
	double current_flow = 1e+10;
	Table2D<double> colorweights = Table2D<double>(img_w,img_h,1.0);
	
	hard_constraints.reset(NONE);
	int imgid=0;
	while(1)
	{
		GraphType * g;
		g = new GraphType(/*estimated # of nodes*/ img_w*img_h, /*estimated # of edges*/ 4*img_w*img_h); 
	    g->add_node(img_w*img_h);    // adding nodes
		
		if(w_smooth>1e-10)
			addsmoothnessterm(g,image,w_smooth,Table2D<bool>(img_w,img_h,true),false);
		Table2D<double> capsource = Table2D<double>(img_w,img_h,0), capsink = Table2D<double>(img_w,img_h,0);
        if((MODE_KERNEL==kmeans_mode)||(MODE_KERNELXY==kmeans_mode)){
			Table2D<double> w_data_OBJ;
			if(MODE_KERNEL==kmeans_mode) w_data_OBJ = filterrgb(vector_img, getROI(init_labeling, OBJ), sigmas);
			if(MODE_KERNELXY==kmeans_mode) w_data_OBJ = filterrgbxy(vector_img, getROI(init_labeling, OBJ), sigmas,sigmax, sigmay);
			double obj_size  = colorweights.sum(getROI(init_labeling, OBJ));
			double obj_sum = w_data_OBJ.sum(getROI(init_labeling, OBJ));
			for(int i=0;i<img_w;i++){
				for(int j=0;j<img_h;j++){
					capsink[i][j] = colorweights[i][j] - 2*w_data_OBJ[i][j]/obj_size + obj_sum / obj_size / obj_size*colorweights[i][j] ;
				}
			}

			Table2D<double> w_data_BKG;
			if(MODE_KERNEL==kmeans_mode) w_data_BKG = filterrgb(vector_img, getROI(init_labeling, BKG), sigmas);
			if(MODE_KERNELXY==kmeans_mode) w_data_BKG = filterrgbxy(vector_img, getROI(init_labeling, BKG), sigmas,sigmax, sigmay);
			double bkg_size = colorweights.sum(getROI(init_labeling, BKG));
			double bkg_sum = w_data_BKG.sum(getROI(init_labeling, BKG));
			for(int i=0;i<img_w;i++){
				for(int j=0;j<img_h;j++){
					capsource[i][j] = colorweights[i][j] - 2*w_data_BKG[i][j]/bkg_size + bkg_sum / bkg_size / bkg_size *colorweights[i][j]  ;
				}
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
		if(imgid>20){
			printf("max iterations.\n");
			break;
		}
	}

	int solution_obj_size = countintable(init_labeling,OBJ);
	cout<<"foreground ratio"<<(double)solution_obj_size/img_w/img_h<<"ratio"<<endl;
	int box_size = countintable(hard_constraints,NONE);
	double errorrate = geterrorcount(init_labeling,gt)/(double)(box_size);
	cout<<"error"<<errorrate<<"error"<<endl;
	//cout<<"fmeasure"<<fmeasure(init_labeling, gt, OBJ)<<"fmeasure\n";
	//cout<<"jaccard"<<jaccard(init_labeling, gt, OBJ)<<"jaccard\n";
	if(returnv!=NULL) *returnv = errorrate;
	if(energyv!=NULL) *energyv = current_flow;
	//return;
	savebinarylabeling(rgbimg, init_labeling,(destdir+string(imgname) +"_"+string("RGBXY")+pch::to_string((int)sigmax)+/*"_s"+pch::to_string(w_smooth)+*/string(".bmp")).c_str());
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
