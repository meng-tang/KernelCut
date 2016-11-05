// Basic utilities
#ifndef _BASICUTIL_H__
#define _BASICUTIL_H__
#include <time.h>
#include <string>
#include <algorithm>
#include <vector>

#include <graph.h>
#include "ezi/Table2D.h"
#include "ezi/Image2D.h"
#include "Image.h"

#include "dirent.h"



double arrayMean(double * vs, int n);
vector<Vect3D> getvectordata(const Table2D<Vect3D> & floatimg, Table2D<bool> ROI);

#define outv(A) cout << #A << ": " << (double)(A) << endl; // output a Variable
#define outs(S) cout<<S<<endl; // output a String

#define INFTY 1e+20
#define EPS 1e-20 

#include <sstream>

namespace pch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}

template<typename T>
inline vector<T> VectorMultiply(const vector<T> & v1, double d);

template<typename T>
inline vector<T> VectorDevide(const vector<T> & v1, double d);

template<typename T>
inline vector<T> VectorSum(const vector<T> & v1, const vector<T> & v2);

template<typename T>
inline vector<T> VectorMinus(const vector<T> & v1, const vector<T> & v2);

template<typename T>
inline T VectorDotProduct(const vector<T> & v1, const vector<T> & v2);

enum Label {NONE=0, OBJ=1, BKG=2, UNKNOWN=3}; 
// UNKNOWN is introduced for monotonic parametric maxflow
// In this case, NONE means out of region of interest.
typedef Graph<double,double,double> GraphType;

// get current date and time as a string
const std::string currentDateTime();
// count certain element in table
template<typename T>
int countintable(const Table2D<T> & table, T t);
template<typename T>
int tablediffcount(const Table2D<T> & t1, const Table2D<T> & t2);
template<typename T>
double fmeasure(const Table2D<T> & solution, const Table2D<T> & gt, T obj_t);
template<typename T>
double jaccard(const Table2D<T> & solution, const Table2D<T> & gt, T obj_t);
template<typename T>
void replaceintable(Table2D<T> & table, T oldt, T newt);

template<typename T>
Table2D<bool> getROI(const Table2D<T> & table, T t);

// count certain element in ROI in table
template<typename T>
int countintableROI(const Table2D<T> & table, T t, const Table2D<bool> & ROI);
// labeling.l ---> OBJ
// any other label ---> BKG
Table2D<Label> replacelabeling(const Table2D<Label> labeling, Label l);
// get labeling according to OBJ color
Table2D<Label> getinitlabeling(const Table2D<int> & initimg, int OBJcolor);
Table2D<Label> getinitlabelingFB(const Table2D<RGB> & initimg, RGB OBJcolor, RGB BKGcolor);
Table2D<int> getinitlabelingMULTI(const Table2D<RGB> & initimg, RGB colors[], int numColor);
// save binary labeling as B&W image
void savebinarylabelingBW(const Table2D<Label> & labeling,string outname);
// save binary labeling
void savebinarylabeling(const Table2D<RGB> & img, const Table2D<Label> & labeling,string outname,bool BW = false);
void savecontour(const Table2D<RGB> & img, const Table2D<Label> & labeling,string outname,bool BW = false);
// save multi labeling
template<typename T>
void savemultilabeling(Table2D<T> & labeling,const char * savedname,RGB * colors,Table2D<RGB> rgbimg);
// get labeling from maxflow instances
bool getlabeling(GraphType * g, Table2D<Label> & m_labeling);
// add smoothness term to the graph
// lambda is the weight of the smoothness term
// ROI is the region of interest
void addsmoothnessterm(GraphType * g, const Image & image, double lambda, 
	const Table2D<bool> & ROI, bool bordersmooth = false);
void savesmoothnessterm(const Image & image, const char * dest_file);
// error rate
double geterrorrate(Table2D<Label> & m_labeling,Table2D<int> & gtimg, int boxsize, int gtOBJcolor=0);
// smoothness cost - weight not applied
double getsmoothnesscost(const Image & image, const Table2D<Label> & m_labeling, bool bordersmoothness = false);
// distance transform based for L2 distance
// positive for object and negative for background
Table2D<double> getDistanceTransform(Table2D<Label> & labeling);

vector<int> getrandomvector(int n);
vector<Point> getrandomvector2dim(int n);
Table2D<Label> complementlabel(Table2D<Label> table);

Table2D<double> getoneDprob(const Table2D<RGB> & rgbimg,int temp_gap=1);

template <class T>
bool readbinfile(Table2D<T> & table, const char * filename, int w, int h);

template <class T>
bool writebinfile(const Table2D<T> & table, const char * filename);

double tablemedian(const Table2D<double> & table);

double tablemedian(const Table2D<double> & table)
{
	int w = table.getWidth();
	int h = table.getHeight();
	vector<double> v(w*h);
	for(int i=0;i<w;i++){
		for(int j=0;j<h;j++){
			v[i*h+j] = table[i][j];
		}	
	}
	std::sort(v.begin(),v.end());
	return v[w*h/2];
}
template<typename T>
inline vector<T> VectorMultiply(const vector<T> & v1, double d)
{
	int n=v1.size();
	vector<T> v(n);
	T returnv = 0;
	for(int i=0;i<n;i++)
		v[i] = v1[i]*d;
	return v;
}

template<typename T>
inline vector<T> VectorDevide(const vector<T> & v1, double d)
{
	int n=v1.size();
	vector<T> v(n);
	T returnv = 0;
	for(int i=0;i<n;i++)
		v[i] = v1[i]/d;
	return v;
}

template<typename T>
inline vector<T> VectorSum(const vector<T> & v1, const vector<T> & v2)
{
	int n=v1.size();
	if(v1.size()!=v2.size())
		cout<<"not the same length"<<endl;
	vector<T> v(n);
	T returnv = 0;
	for(int i=0;i<n;i++)
		v[i] = v1[i] + v2[i];
	return v;
}

template<typename T>
inline vector<T> VectorMinus(const vector<T> & v1, const vector<T> & v2)
{
	int n=v1.size();
	if(v1.size()!=v2.size())
		cout<<"not the same length"<<endl;
	vector<T> v(n);
	T returnv = 0;
	for(int i=0;i<n;i++)
		v[i] = v1[i] - v2[i];
	return v;
}

template<typename T>
inline T VectorDotProduct(const vector<T> & v1, const vector<T> & v2)
{
	int n=v1.size();
	T returnv = 0;
	for(int i=0;i<n;i++)
		returnv = returnv + v1[i]*v2[i];
	return returnv;
}
Table2D<Label> complementlabel(Table2D<Label> table)
{
	int img_w = table.getWidth();
	int img_h = table.getHeight();
	Table2D<Label> complement(img_w,img_h,NONE);
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			if(table[i][j]==OBJ)
				complement[i][j]=BKG;
			else if(table[i][j]==BKG)
				complement[i][j]=OBJ;
		}
	}
	return complement;
}
void savetableasgrayimage(Table2D<double> table, const char * imgname);
void savetableasgrayimage(Table2D<double> table, const char * imgname)
{
	int img_w = table.getWidth();
	int img_h = table.getHeight();
	double tablemax = table.getMax();
	double tablemin = table.getMin();
	Table2D<int> tmp(img_w,img_h);
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			int grayvalue = (table[i][j]-tablemin)/(tablemax-tablemin)*255;
			//int grayvalue = table[i][j]*10;
			if(grayvalue>255) grayvalue=255;
			if(grayvalue<0) grayvalue=0;
			tmp[i][j] = grayvalue;
		}
	}
	if( saveImage(tmp, imgname))
	    cout<<"saved into "<<imgname<<endl;
}

void savetableascolorimage(Table2D<double> table, const char * imgname)
{
	int img_w = table.getWidth();
	int img_h = table.getHeight();
	double tablemax = table.getMax();
	double tablemin = table.getMin();
	Table2D<RGB> tmp(img_w,img_h);
	// matlab jet colormap: colormap(jet); cmap = colormap; for i=1:64 disp([num2str(cmap(i,1)) ',' num2str(cmap(i,2)) ',' num2str(cmap(i,3)) ',' ]) end
	double cmap[64][3] = {0,0,0.5625,
0,0,0.625,
0,0,0.6875,
0,0,0.75,
0,0,0.8125,
0,0,0.875,
0,0,0.9375,
0,0,1,
0,0.0625,1,
0,0.125,1,
0,0.1875,1,
0,0.25,1,
0,0.3125,1,
0,0.375,1,
0,0.4375,1,
0,0.5,1,
0,0.5625,1,
0,0.625,1,
0,0.6875,1,
0,0.75,1,
0,0.8125,1,
0,0.875,1,
0,0.9375,1,
0,1,1,
0.0625,1,0.9375,
0.125,1,0.875,
0.1875,1,0.8125,
0.25,1,0.75,
0.3125,1,0.6875,
0.375,1,0.625,
0.4375,1,0.5625,
0.5,1,0.5,
0.5625,1,0.4375,
0.625,1,0.375,
0.6875,1,0.3125,
0.75,1,0.25,
0.8125,1,0.1875,
0.875,1,0.125,
0.9375,1,0.0625,
1,1,0,
1,0.9375,0,
1,0.875,0,
1,0.8125,0,
1,0.75,0,
1,0.6875,0,
1,0.625,0,
1,0.5625,0,
1,0.5,0,
1,0.4375,0,
1,0.375,0,
1,0.3125,0,
1,0.25,0,
1,0.1875,0,
1,0.125,0,
1,0.0625,0,
1,0,0,
0.9375,0,0,
0.875,0,0,
0.8125,0,0,
0.75,0,0,
0.6875,0,0,
0.625,0,0,
0.5625,0,0,
0.5,0,0};
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
		    int idx = table[i][j] + 32;
		    idx = min(idx,63);idx=max(idx,0);
			tmp[i][j].r = cmap[idx][0]*255;tmp[i][j].g = cmap[idx][1]*255;tmp[i][j].b = cmap[idx][2]*255;
		}
	}
	if( saveImage(tmp, imgname))
	    cout<<"saved into "<<imgname<<endl;
}
const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    tstruct = *localtime(&now);
	int h = tstruct.tm_hour;
	int m = tstruct.tm_min;
	int s = tstruct.tm_sec;
	string time_s;
	char buff[4];
	sprintf( buff, "%d", h );
	time_s += buff;
	sprintf( buff, "%d", m );
	time_s += buff;
	sprintf( buff, "%d", s );
	time_s += buff;
    return time_s;
}

template<typename T>
int countintable(const Table2D<T> & table, T t)
{
	int table_w = table.getWidth();
	int table_h = table.getHeight();
	int tsize = 0;
	for(int y=0; y<table_h; y++) 
	{
		for(int x=0; x<table_w; x++) 
		{ 
			if(table[x][y]==t) // certain element t
				tsize++;
		}
	}
	return tsize;
}

template<typename T>
int tablediffcount(const Table2D<T> & t1, const Table2D<T> & t2)
{
	int table_w = t1.getWidth();
	int table_h = t1.getHeight();
	int diffcount = 0;
	for(int y=0; y<table_h; y++) 
	{
		for(int x=0; x<table_w; x++) 
		{ 
			if(!(t1[x][y]==t2[x][y])) // certain element t
				diffcount++;
		}
	}
	return diffcount;
}

template<typename T>
double fmeasure(const Table2D<T> & solution, const Table2D<T> & gt, T obj_t)
{
	int w = gt.getWidth();
	int h = gt.getHeight();
	int tpcount = 0;
	for(int y=0; y<h; y++) 
	{
		for(int x=0; x<w; x++) 
		{ 
			if(gt[x][y]==obj_t && solution[x][y]==obj_t) // certain element t
				tpcount++;
		}
	}
	int solution_size = countintable(solution,obj_t);
	int gt_size = countintable(gt,obj_t);
	double precision = (double)tpcount / gt_size;
	double recall = (solution_size != 0 ? (double)tpcount / solution_size:1.0);
	return 2*precision*recall / (precision + recall);
}

template<typename T>
double jaccard(const Table2D<T> & solution, const Table2D<T> & gt, T obj_t)
{
	int w = gt.getWidth();
	int h = gt.getHeight();
	int unioncount = 0;
	int interseccount = 0;
	for(int y=0; y<h; y++) 
	{
		for(int x=0; x<w; x++) 
		{ 
			if(gt[x][y]==obj_t && solution[x][y]==obj_t) // certain element t
				interseccount++;
			if(gt[x][y]==obj_t || solution[x][y]==obj_t) // certain element t
				unioncount++;
		}
	}
	
	return (double)interseccount / unioncount;
}

template<typename T>
void replaceintable(Table2D<T> & table, T oldt, T newt)
{
	int table_w = table.getWidth();
	int table_h = table.getHeight();
	for(int y=0; y<table_h; y++) 
	{
		for(int x=0; x<table_w; x++) 
		{ 
			if(table[x][y]==oldt) // certain element t
				table[x][y]=newt;
		}
	}
}

template<typename T>
Table2D<bool> getROI(const Table2D<T> & table, T t)
{
	int table_w = table.getWidth();
	int table_h = table.getHeight();
	Table2D<bool> ROI(table_w,table_h,false);
	for(int y=0; y<table_h; y++) 
	{
		for(int x=0; x<table_w; x++) 
		{ 
			if(table[x][y]==t) // certain element t
				ROI[x][y]=true;
		}
	}
	return ROI;
}

template<typename T>
int countintableROI(const Table2D<T> & table, T t, const Table2D<bool> & ROI)
{
	int table_w = table.getWidth();
	int table_h = table.getHeight();
	int tsize = 0;
	for(int y=0; y<table_h; y++) 
	{
		for(int x=0; x<table_w; x++) 
		{ 
			if((table[x][y]==t)&&(ROI[x][y])) // certain element t
				tsize++;
		}
	}
	return tsize;
}

// labeling.l ---> OBJ
// any other label ---> BKG
Table2D<Label> replacelabeling(const Table2D<Label> labeling, Label l)
{
	Table2D<Label> newlabeling(labeling.getWidth(),labeling.getHeight(),BKG);
	for(int i=0;i<labeling.getWidth();i++)
	{
		for(int j=0;j<labeling.getHeight();j++)
		{
			if(labeling[i][j]==l)
				newlabeling[i][j]=OBJ;
		}
	}
	return newlabeling;
}

Table2D<Label> getinitlabeling(const Table2D<int> & initimg, int OBJcolor)
{
	int img_w = initimg.getWidth();
	int img_h = initimg.getHeight();
	Table2D<Label> initlabeling(img_w,img_h);
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			if(initimg[i][j]==OBJcolor)
				initlabeling[i][j] = OBJ;
			else
				initlabeling[i][j] = BKG;
		}
	}
	return initlabeling;
}

Table2D<Label> getinitlabeling2(const Table2D<RGB> & initimg)
{
	int img_w = initimg.getWidth();
	int img_h = initimg.getHeight();
	Table2D<Label> initlabeling(img_w,img_h);
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			if(initimg[i][j]==white)
				initlabeling[i][j] = BKG;
			else
				initlabeling[i][j] = OBJ;
		}
	}
	return initlabeling;
}

Table2D<Label> getinitlabeling(const Table2D<int> & initimg, int OBJcolor, int BKGcolor)
{
	int img_w = initimg.getWidth();
	int img_h = initimg.getHeight();
	Table2D<Label> initlabeling(img_w,img_h);
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			if(initimg[i][j]==OBJcolor)
				initlabeling[i][j] = OBJ;
			else if(initimg[i][j]==BKGcolor)
				initlabeling[i][j] = BKG;
			else
				initlabeling[i][j] = NONE;
		}
	}
	return initlabeling;
}

int geterrorcount(const Table2D<Label> & labeling, const Table2D<Label> & gt)
{
	int img_w = labeling.getWidth();
	int img_h = labeling.getHeight();
	int errorcount = 0;
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			if((gt[i][j]==OBJ)&&(labeling[i][j]!=OBJ))
				errorcount++;
			else if((gt[i][j]==BKG)&&(labeling[i][j]!=BKG))
				errorcount++;
		}
	}
	return errorcount;
}

Table2D<Label> getinitlabelingFB(const Table2D<RGB> & initimg, RGB OBJcolor, RGB BKGcolor)
{
	int img_w = initimg.getWidth();
	int img_h = initimg.getHeight();
	Table2D<Label> initlabeling(img_w,img_h);
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			if(initimg[i][j]==OBJcolor)
				initlabeling[i][j] = OBJ;
			else if(initimg[i][j]==BKGcolor)
				initlabeling[i][j] = BKG;
			else
				initlabeling[i][j] = NONE;
		}
	}
	return initlabeling;
}

Table2D<int> getinitlabelingMULTI(const Table2D<RGB> & initimg, RGB colors[], int numColor)
{
	int img_w = initimg.getWidth();
	int img_h = initimg.getHeight();
	Table2D<int> initlabeling(img_w,img_h,-1);
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
		    for(int c=0;c<numColor;c++){
		        if(initimg[i][j]==colors[c])
		            initlabeling[i][j] = c;
		    }
		}
	}
	return initlabeling;
}

void savebinarylabelingBW(const Table2D<Label> & labeling,string outname)
{
	int img_w = labeling.getWidth();
	int img_h = labeling.getHeight();
	Table2D<RGB> tmp(img_w,img_h);
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			if(labeling[i][j]==BKG)
				tmp[i][j] = white;
			else if(labeling[i][j]==OBJ)
				tmp[i][j] = black;
			else
				tmp[i][j] = gray;
		}
	}
	if(saveImage(tmp, outname.c_str()))
	    cout<<"saved into: "<<outname<<endl;
}

void savebinarylabeling(const Table2D<RGB> & img, const Table2D<Label> & labeling,string outname,bool BW)
{
	int img_w = labeling.getWidth();
	int img_h = labeling.getHeight();
	Table2D<RGB> tmp(img_w,img_h);
	//BW = true;
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			if(labeling[i][j]==OBJ)
			{
				if(BW) tmp[i][j] = black;
				else tmp[i][j] = img[i][j];
			}
			else
				tmp[i][j] = white;
		}
	}
	if(saveImage(tmp, outname.c_str()))
	    cout<<"saved into: "<<outname<<endl;
}

void savecontour(const Table2D<RGB> & img, const Table2D<Label> & labeling,string outname,bool BW)
{
	int img_w = labeling.getWidth();
	int img_h = labeling.getHeight();
	Table2D<RGB> tmp = img;
	for(int i=1;i<img_w-1;i++)
	{
		for(int j=1;j<img_h-1;j++)
		{
			if((labeling[i][j]==OBJ)&&(labeling[i-1][j]==BKG))
				tmp[i][j] = red;
			if((labeling[i][j]==OBJ)&&(labeling[i+1][j]==BKG))
				tmp[i][j] = red;
			if((labeling[i][j]==OBJ)&&(labeling[i][j-1]==BKG))
				tmp[i][j] = red;
			if((labeling[i][j]==OBJ)&&(labeling[i][j+1]==BKG))
				tmp[i][j] = red;
		}
	}
	saveImage(tmp, outname.c_str());
	cout<<"saved into: "<<outname<<endl;
}

template<typename T>
void savemultilabeling(Table2D<T> & labeling,const char * outname,RGB * colors,Table2D<RGB> rgbimg)
{
	int img_w = labeling.getWidth();
	int img_h = labeling.getHeight();
	Table2D<RGB> tmp(img_w,img_h);
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
		    if(labeling[i][j]>=0)
			    tmp[i][j] = colors[labeling[i][j]];
			else
			    tmp[i][j] = black;
		}
	}
	saveImage(tmp, outname);
	cout<<"saved into: "<<outname<<endl;
}

bool getlabeling(GraphType * g, Table2D<Label> & m_labeling)
{
	int img_w = m_labeling.getWidth();
	int img_h = m_labeling.getHeight();
	int n=0;
	m_labeling.reset(NONE);
	int sumobj =0, sumbkg = 0;
	for (int y=0; y<img_h; y++) 
	{
		for (int x=0; x<img_w; x++) 
		{ 
			n = x+y*img_w;
			if(g->what_segment(n) == GraphType::SOURCE)
			{
				m_labeling[x][y]=OBJ;
				sumobj++;
			}
			else if(g->what_segment(n) == GraphType::SINK)
			{
				m_labeling[x][y]=BKG;
				sumbkg++;
			}
			//else
				//exit(-1);
		}
	}
	if((sumobj==0)||(sumbkg==0))
		return false;
	else
		return true;
}

// add smoothness term to the graph
// lambda is the weight of the smoothness term
// ROI is the region of interest
void addsmoothnessterm(GraphType * g, const Image & image, double lambda,
	const Table2D<bool> & ROI, bool bordersmooth)
{
	int img_w = image.img_w;
	int img_h = image.img_h;
	// number of neighboring pairs of pixels
	int numNeighbor = image.pointpairs.size();
	// n-link - smoothness term
	int node_id1 =0, node_id2 =0;
	
	for(int i=0;i<numNeighbor;i++)
	{
		PointPair pp = image.pointpairs[i];
		node_id1 = pp.p1.x+pp.p1.y*img_w;
		node_id2 = pp.p2.x+pp.p2.y*img_w;
		// if the two points are inside region of interest.
		if(ROI[pp.p1]&&ROI[pp.p2])
		{
			//double v = fn(dI(image.img[pp.p1],image.img[pp.p2]),lambda,image.sigma_square)/(pp.p1-pp.p2).norm();   
			double v = lambda*image.smoothnesscosts[i];
			g->add_edge(node_id1,node_id2,v,v);
		}
	}
	if(bordersmooth)
	{
		for(int i=0;i<image.img_w;i++)
		{
			for(int j=0;j<image.img_h;j++)
				if((i==0)||(i==img_w-1)||(j==0)||(j==img_h-1)) // border
					g->add_tweights(i+j*img_w,0,lambda);
		}
	}
}

// save smoothness term as an image
void savesmoothnessterm(const Image & image, const char * dest_file)
{
	int img_w = image.img_w;
	int img_h = image.img_h;
	// number of neighboring pairs of pixels
	int numNeighbor = image.pointpairs.size();
	// n-link - smoothness term
	int node_id1 =0, node_id2 =0;
	Table2D<double> smoothnesscosts(img_w,img_h,0);
	for(int i=0;i<numNeighbor;i++)
	{
		PointPair pp = image.pointpairs[i];
		double v = image.smoothnesscosts[i];
		smoothnesscosts[pp.p1] += v;
		smoothnesscosts[pp.p2] += v;
	}
	savetableasgrayimage(smoothnesscosts, dest_file);
}

double geterrorrate(Table2D<Label> & m_labeling,Table2D<int> & gtimg, int boxsize, int gtOBJcolor)
{
	double errorrate = 0 ;
	int errornum = 0;
	for(int j=1;j<gtimg.getHeight()-1;j++)
	{
		for(int i=1;i<gtimg.getWidth()-1;i++)
		{
			if((gtimg[i][j]==gtOBJcolor)&&(m_labeling[i][j]==OBJ))
				errornum++;
			else if((gtimg[i][j]==(255-gtOBJcolor))&&(m_labeling[i][j]==BKG))
				errornum++;
		}
	}
	cout<<"errornum / boxsize "<<errornum<<' '<<boxsize<<endl;
	errorrate = (double)errornum / boxsize;
	return errorrate;
}

double getsmoothnesscost(const Image & image, const Table2D<Label> & m_labeling, bool bordersmoothness)
{
	int img_w = image.img_w;
	int img_h = image.img_h;
	double smoothenergy = 0;
	// number of neighboring pairs of pixels
	int numNeighbor = image.pointpairs.size();
	for(int i=0;i<numNeighbor;i++)
	{
		PointPair pp = image.pointpairs[i];
		if((m_labeling[pp.p1]!=NONE)&&(m_labeling[pp.p2]!=NONE))
		{
			if(m_labeling[pp.p1]!=m_labeling[pp.p2])
				smoothenergy += image.smoothnesscosts[i];
		}
		if((m_labeling[pp.p1]==NONE)||(m_labeling[pp.p2]==NONE))
			exit(-1);
	}
	if(bordersmoothness)
	{
		for(int i=0;i<image.img_w;i++)
		{
			for(int j=0;j<image.img_h;j++)
				if((i==0)||(i==img_w-1)||(j==0)||(j==img_h-1)) // border
					if(m_labeling[i][j]==OBJ)
						smoothenergy = smoothenergy +1;
		}
	}
	return smoothenergy;
}

////////////////////////////////////////////////////////////////////////
// getDistanceTransform is a function that computes distance 
// transform. It implements forward-backward pass algorithm with
// windows that approximate Euclidean Metric.
Table2D<double> getDistanceTransform(Table2D<Label> & labeling)
{
	int img_w = labeling.getWidth();
	int img_h = labeling.getHeight();
	Table2D<double> dt(img_w,img_h,INFTY); // distance transform
	//find all the edges
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			if((i==0)||(i==img_w-1)||(j==0)||(j==img_h-1)) // border
			{
				if(labeling[i][j]==OBJ)
					dt[i][j] = 0;
			}
			else // non-border
			{
				if((labeling[i][j]==OBJ)&&
				((labeling[i+1][j]==BKG)||(labeling[i-1][j]==BKG)||(labeling[i][j+1]==BKG)||(labeling[i][j-1]==BKG)))
					dt[i][j] = 0;
			}
		}
	}
	//forward pass
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			Point p(i,j);
			Point q;
			q=p+Point(-1,-1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.414+dt[q]);
			q=p+Point(0,-1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.0+dt[q]);
			q=p+Point(1,-1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.414+dt[q]);
			q=p+Point(-1,0);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.0+dt[q]);
		}
	}
	//backward pass
	for(int i=img_w-1;i>=0;i--)
	{
		for(int j=img_h-1;j>=0;j--)
		{
			Point p(i,j);
			Point q;
			q=p+Point(1,1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.414+dt[q]);
			q=p+Point(0,1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.0+dt[q]);
			q=p+Point(-1,1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.414+dt[q]);
			q=p+Point(1,0);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.0+dt[q]);
		}
	}

	// second time forward and backward

	//forward pass
	for(int j=0;j<img_h;j++)
	{
		for(int i=0;i<img_w;i++)
		{
			Point p(i,j);
			Point q;
			q=p+Point(-1,-1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.414+dt[q]);
			q=p+Point(0,-1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.0+dt[q]);
			q=p+Point(1,-1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.414+dt[q]);
			q=p+Point(-1,0);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.0+dt[q]);
		}
	}
	//backward pass
	for(int j=img_h-1;j>=0;j--)
	{
		for(int i=img_w-1;i>=0;i--)
		{
			Point p(i,j);
			Point q;
			q=p+Point(1,1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.414+dt[q]);
			q=p+Point(0,1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.0+dt[q]);
			q=p+Point(-1,1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.414+dt[q]);
			q=p+Point(1,0);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.0+dt[q]);
		}
	}

	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			if(labeling[i][j]==OBJ)
				dt[i][j] = dt[i][j]+0.5;
			else if(labeling[i][j]==BKG)
				dt[i][j] = -(dt[i][j]-0.5);
		}
	}
	return dt;
}

vector<int> getrandomvector(int n)
{
	vector<int> v(n,0);
	for (int i=0;i<n;i++)
		v[i]=i;
	int p;
	int tmp;

	for (int i=n-1;i>0;i--)
	{
		p=rand()%i;
		tmp=v[p];
		v[p]=v[i];
		v[i]=tmp;
	}
	return v;
}

vector<Point> getrandomvector2dim(int n)
{
	vector<Point> v;
	for(int i=0;i<n-1;i++)
		for(int j=i+1;j<n;j++)
			v.push_back(Point(i,j));
	int size_v = v.size();
	//outv(size_v);

	for (int i=size_v-1;i>0;i--)
	{
		int p=rand()%i;
		Point tmp=v[p];
		v[p]=v[i];
		v[i]=tmp;
	}
	return v;
}

Table2D<double> getoneDprob(const Table2D<RGB> & rgbimg,int temp_gap)
{
	int w = rgbimg.getWidth();
	int h = rgbimg.getHeight();
	Table2D<double> colorcounts(256,3,0);
	for(int i=0;i<w;i++){
		for(int j=0;j<h;j++){
			RGB c = rgbimg[i][j];
			c.r = c.r - c.r%temp_gap;
			c.g = c.g - c.g%temp_gap;
			c.b = c.b - c.b%temp_gap;
			colorcounts[c.r][0] += 1;
			colorcounts[c.g][1] += 1;
			colorcounts[c.b][2] += 1;
		}
	}
	Table2D<double> oneDprob(256,3,0);
	for(int c=0;c<3;c++){
		int total = 0;
		oneDprob[0][c] = (colorcounts[0][c]+colorcounts[1][c]) / 2;
		total += oneDprob[0][c];
		for(int i=1;i<255;i++){
			oneDprob[i][c] = (colorcounts[i-1][c]+colorcounts[i][c]+colorcounts[i+1][c]) / 3;
			total += oneDprob[i][c];
		}
		oneDprob[255][c] = (colorcounts[254][c]+colorcounts[255][c]) / 2;
		total += oneDprob[255][c];
		for(int i=0;i<256;i++){
			oneDprob[i][c] = oneDprob[i][c] / total;
		}
	}
	return oneDprob;
}

template <class T>
bool readbinfile(Table2D<T> & table, const char * filename, int w, int h){
	table.resize(w,h);
	//printf("read bin file %s\n",filename);
	FILE * pFile = fopen ( filename , "rb" );
    if (pFile==NULL) { fputs ("File error",stderr); exit (1);}
	T * buffer = new T[h];
	for(int i=0;i<w;i++){
		fread(buffer,sizeof(T),h,pFile);
		for(int j=0;j<h;j++){
			table[i][j] = buffer[j];
		}
	}
    fclose (pFile);
	delete [] buffer;
	//printf("bin closed\n");
	return true;
}

template <class T>
bool writebinfile(const Table2D<T> & table, const char * filename){
	int w = table.getWidth();
	int h = table.getHeight();
	printf("write bin file %s\n",filename);
	FILE * pFile = fopen ( filename , "wb" );
    if (pFile==NULL) { fputs ("File error",stderr); exit (1);}
	T * buffer = new T[h];
	for(int i=0;i<w;i++){
		for(int j=0;j<h;j++){
			buffer[j] = table[i][j]; 
		}
		fwrite(buffer,sizeof(T),h,pFile);
	}
    fclose (pFile);
	delete [] buffer;
	printf("bin closed\n");
	return true;
}

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

#endif
