#pragma once
#include "ezi/Image2D.h"
#include "ezi/Table2D.h"
#include "SparseMatrix.h"
#include <assert.h>


// get average sigma_square for smoothness term
double getsmoothnessterm(const Table2D<RGB> &img,vector<PointPair> & pointpairs,int connecttype);
void rgb2indeximg(Table2D<int> & indeximg,const Table2D<RGB> & img,double channelstep);
int getcompactlabel(Table2D<int> & colorlabel,double channelstep,vector<int> & compacthist);
// edge-constrast sensitive smoothness penalty
inline double fn(const double dI, double lambda,double sigma_square);

class Image{
public:
	Image();
	Image(Table2D<RGB> img_, double channelstep, const char * imgname_ = "", int connecttype_ = 16);
	Image(const char * imgpath, const char * imgname_, double channelstep, int connecttype_);
	void computesmoothnesscost();
	void print();

    void addboxsmooth(Table2D<bool> box);
	Table2D<RGB> img;
	const char * imgname;
	int img_w;
	int img_h;
	int img_size;
	double sigma_square;
	vector<PointPair> pointpairs;
	vector<double> smoothnesscosts;
	int colorbinnum;
	Table2D<int> colorlabel;
	vector<int> compacthist; // color bin histogram (whole image)
	int connecttype; // can be 4 or 8 or 16

	vector<Trituple<double> > pair_arcs;
};

Image::Image()
{

}

void Image::computesmoothnesscost()
{
	// number of neighboring pairs of pixels
	int numNeighbor = pointpairs.size();
	smoothnesscosts = vector<double>(numNeighbor);
	for(int i=0;i<numNeighbor;i++)
	{
		PointPair pp = pointpairs[i];
		smoothnesscosts[i] = fn(dI(img[pp.p1],img[pp.p2]),1.0,sigma_square)/(pp.p1-pp.p2).norm();
	}

	pair_arcs = vector< Trituple<double> >();
	for(int i=0;i<numNeighbor;i++)
	{
		Point p1 = pointpairs[i].p1;
		Point p2 = pointpairs[i].p2;
		pair_arcs.push_back(Trituple<double>(p1.x+p1.y*img_w,p2.x+p2.y*img_w,smoothnesscosts[i]));
		pair_arcs.push_back(Trituple<double>(p2.x+p2.y*img_w,p1.x+p1.y*img_w,smoothnesscosts[i]));
	}

}

void Image::addboxsmooth(Table2D<bool> box)
{
    int numNeighbor = pointpairs.size();
    for(int i=0;i<numNeighbor;i++)
	{
		PointPair pp = pointpairs[i];
		Point p1 = pointpairs[i].p1;
		Point p2 = pointpairs[i].p2;
		if( (box[p1]==false && box[p2] == true) || (box[p1]==true && box[p2] == false) )
		smoothnesscosts[i] = 1000000.0 / (pp.p1-pp.p2).norm();
	}
}

Image::Image(Table2D<RGB> img_, double channelstep, const char * imgname_, int connecttype_)
{
	Assert((connecttype_==4)||(connecttype_==8)||(connecttype_==16),"wrong connect type!");
	imgname = imgname_;
	img = Table2D<RGB>(img_);
	img_w = img.getWidth();
	img_h = img.getHeight();
	img_size = img_w * img_h;

	connecttype = connecttype_;
	sigma_square = getsmoothnessterm(img,pointpairs,connecttype);
	//sigma_square = sigma_square*2;
	
	colorlabel= Table2D<int>(img_w,img_h);
	rgb2indeximg(colorlabel,img,channelstep);
	colorbinnum = getcompactlabel(colorlabel,channelstep,compacthist);

	computesmoothnesscost();

}
Image::Image(const char * imgpath, const char * imgname_, double channelstep, int connecttype_)
{
	Table2D<RGB> img_ = loadImage<RGB>(imgpath);
	new (this)Image(img_, channelstep, imgname_, connecttype_);
}
void Image::print()
{
	cout<<"sigma_square: "<<sigma_square<<endl;
	cout<<"image width: "<<img_w<<endl;
	cout<<"image height: "<<img_h<<endl;
	cout<<"# of pairs of neighborhoods"<<pointpairs.size()<<endl;
	cout<<"# of color bins"<<colorbinnum<<endl;
	cout<<"size of compact hist: "<<compacthist.size()<<endl;
	cout<<"sigma_square = "<<sigma_square<<endl;
}

double getsmoothnessterm(const Table2D<RGB> &img,vector<PointPair> & pointpairs, int connecttype)
{
	int node_id = 0;
	int img_w = img.getWidth();
	int img_h = img.getHeight();
	double sigma_sum = 0;
	double sigma_square_count = 0;
	Point kernelshifts [] = {Point(1,0),Point(0,1),Point(1,1),Point(1,-1),
		Point(2,-1),Point(2,1),Point(1,2),Point(-1,2),};
	for (int y=0; y<img_h; y++) // adding edges (n-links)
	{
		for (int x=0; x<img_w; x++) 
		{ 
			Point p(x,y);
			for(int i=0;i<connecttype/2;i++)
			{
				Point q = p + kernelshifts[i];
				if(img.pointIn(q))
				{
					sigma_sum += dI(img[p],img[q]);
					pointpairs.push_back(PointPair(p,q));
					sigma_square_count ++;
				}
			}
		}
	}
	return sigma_sum/sigma_square_count;
}

void rgb2indeximg(Table2D<int> & indeximg,const Table2D<RGB> & img,double channelstep)
{
	RGB rgb_v;
	int r_idx =0, g_idx = 0, b_idx = 0, idx =0;
	int channelbin = (int)ceil(256.0/channelstep);
	for(unsigned int j=0;j<img.getHeight();j++)
	{
		for(unsigned int i=0;i<img.getWidth();i++)
		{
			rgb_v = img[i][j];
			r_idx = (int)(rgb_v.r/channelstep);
			g_idx = (int)(rgb_v.g/channelstep);
			b_idx = (int)(rgb_v.b/channelstep);
			idx = r_idx + g_idx*channelbin+b_idx*channelbin*channelbin;
			indeximg[i][j] = idx;
		}
	}
}

int getcompactlabel(Table2D<int> & colorlabel,double channelstep,vector<int> & compacthist)
{
	int channelbin = ceil(256/channelstep);
	vector<int> colorhist(channelbin*channelbin*channelbin,0);
	for(unsigned int j=0;j<colorlabel.getHeight();j++)
	{
		for(unsigned int i=0;i<colorlabel.getWidth();i++)
		{
			colorhist[colorlabel[i][j]] = colorhist[colorlabel[i][j]]+1;
		}
	}
	
	vector<int> correspondence(colorhist.size(),-1);
	int compactcount = 0;
	for(int i=0;i<colorhist.size();i++)
	{
		if(colorhist[i]!=0)
		{
			compacthist.push_back(colorhist[i]);
			correspondence[i] = compactcount;
			compactcount++;
		}
	}
	for(int j=0;j<colorlabel.getHeight();j++)
	{
		for(int i=0;i<colorlabel.getWidth();i++)
		{
			colorlabel[i][j] = correspondence[colorlabel[i][j]];
		}
	}
	return compacthist.size();
}

inline double fn(const double dI, double lambda,double sigma_square) 
{
	return lambda*(exp(-dI/2/sigma_square));
}
