#ifndef _PPBCGRABCUT_H_
#define _PPBCGRABCUT_H_
// PPBC for GrabCut (using histogram)
#include "PPBCBase.h"

#define GRABCUT_RGB   1
#define GRABCUT_RGBXY 2


Table2D<double> gethistdist(const Table2D<RGB> & rgbimg,int binsize, Table2D<bool> ROI);
Table2D<double> getxyhistdist(const Table2D<RGB> & rgbimg,int binsize, int xybinsize, Table2D<bool> ROI);

class PPBCgrabcut: public PPBCBase{
public:
	PPBCgrabcut(const Image & image_, double binsize_, int xybinsize_, double w_entropy_, double w_smooth_, const Table2D<Label> & hardconstraints_, double alpha_ = 1e-10)
	:image(image_),binsize(binsize_),w_entropy(w_entropy_),w_smooth(w_smooth_),alpha(alpha_){
		img_w = image.img_w;
		img_h = image.img_h;
		hardconstraints = hardconstraints_;
		xybinsize = xybinsize_;
		if(xybinsize > 0) GrabCutMode = GRABCUT_RGBXY;
		else GrabCutMode = GRABCUT_RGB;
	};
	virtual double computeenergy(const Table2D<Label> & labeling);
	virtual void updatemodel(); // collect statistics about current labeling
	virtual void * parabasegraph(double para_, UnknownRegion * unknownregion_p,
		double * flowoffset_p, vector<Point> * node_corr_p);
	virtual double getnewpara(ParaInterval paraInterval);
	
	int GrabCutMode;
private:
	Image image;
	double binsize;
	int xybinsize;
	double w_entropy;
	double w_smooth;
	Table2D<Label> hardconstraints;
	// unary auxiliary function
	Table2D<double> capsource;
	Table2D<double> capsink;
	double alpha;
};

double PPBCgrabcut::computeenergy(const Table2D<Label> & labeling)
{
	Table2D<double> obj_prob, bkg_prob;
	if(GRABCUT_RGB==GrabCutMode){
	    obj_prob = gethistdist(image.img,binsize, getROI(labeling,OBJ));
	    bkg_prob = gethistdist(image.img,binsize, getROI(labeling,BKG));
	}
	else if(GRABCUT_RGBXY==GrabCutMode){
	    obj_prob = getxyhistdist(image.img,binsize, xybinsize,getROI(labeling,OBJ));
	    bkg_prob = getxyhistdist(image.img,binsize, xybinsize,getROI(labeling,BKG));
	}
	double entropy_e =0;
	for(int i=0;i<img_w;i++){
	    for(int j=0;j<img_h;j++){
	        if(labeling[i][j]==OBJ)
	            entropy_e += -log(obj_prob[i][j]);
	        else if(labeling[i][j]==BKG)
	            entropy_e += -log(bkg_prob[i][j]);
	    }
	}

	if(w_smooth <1e-10)
		return entropy_e*w_entropy;
	double smoothenergy = 0;
	// number of neighboring pairs of pixels
	int numNeighbor = image.pointpairs.size();
	// n-link - smoothness term
	int node_id1 =0, node_id2 =0;
	for(int i=0;i<numNeighbor;i++)
	{
		PointPair pp = image.pointpairs[i];
		node_id1 = pp.p1.x+pp.p1.y*img_w;
		node_id2 = pp.p2.x+pp.p2.y*img_w;
		if((labeling[pp.p1]!=NONE)&&(labeling[pp.p2]!=NONE))
		{
			if(labeling[pp.p1]!=labeling[pp.p2])
			{
				double v = image.smoothnesscosts[i];
				smoothenergy += v;
			}
		}
	}
	//printf("grabcut energy: %.2f + %.2f = %.2f\n",entropy_e*w_entropy,smoothenergy * w_smooth,entropy_e*w_entropy + smoothenergy * w_smooth);
	return entropy_e*w_entropy + smoothenergy * w_smooth;
}

void PPBCgrabcut::updatemodel()
{
	capsource = Table2D<double>(img_w,img_h,0), capsink = Table2D<double>(img_w,img_h,0);
	Table2D<double> obj_prob, bkg_prob;
	if(GRABCUT_RGB==GrabCutMode){
	    obj_prob = gethistdist(image.img,binsize, getROI(initlabeling,OBJ));
	    bkg_prob = gethistdist(image.img,binsize, getROI(initlabeling,BKG));
	}
	else if(GRABCUT_RGBXY==GrabCutMode){
	    obj_prob = getxyhistdist(image.img,binsize, xybinsize,getROI(initlabeling,OBJ));
	    bkg_prob = getxyhistdist(image.img,binsize, xybinsize,getROI(initlabeling,BKG));
	}
	else{exit(-1);}
	
	for(int i=0;i<img_w;i++){
		for(int j=0;j<img_h;j++){
			capsink[i][j]   = -log(obj_prob[i][j]*(1-alpha)+alpha*(1.0/(256*256*256/binsize/binsize/binsize)));
			capsource[i][j] = -log(bkg_prob[i][j]*(1-alpha)+alpha*(1.0/(256*256*256/binsize/binsize/binsize)));
		}
	}

	for(int j=0;j<img_h;j++)
	{
		for(int i=0;i<img_w;i++)
		{
			if(BKG==hardconstraints[i][j]) // hard constraints to background
				capsink[i][j]=1e20;
			else if(OBJ==hardconstraints[i][j]) // hard constraints to foreground
				capsource[i][j]=1e20;
		}
	}
}

void * PPBCgrabcut::parabasegraph(double para_, UnknownRegion * unknownregion_p,
	double * flowoffset_p, vector<Point> * node_corr_p)
{
	if(unknownregion_p==NULL){
	GraphType * g = new GraphType(img_w*img_h, 4*img_w*img_h);
	g->add_node(img_w*img_h);    // adding nodes
	// add smoothness term
	if(w_smooth > 1e-10)
	    addsmoothnessterm(g, image, w_smooth, ROI,false);
	// add unary term
	for(int j=0;j<img_h;j++)
	{
		for(int i=0;i<img_w;i++)
		{
			if(ROI[i][j])
			{
				int n = i+j*img_w;
				g->add_tweights(n,capsource[i][j],capsink[i][j]+para_);
			}
		}
	}
	return g;}
}


double PPBCgrabcut::getnewpara(ParaInterval interval)
{
	// ballooning cost for lower parameter solution
	int ssize_low = interval.bplow.ssize;
	double ballooncost_low = ssize_low*interval.bplow.para;
	double entropy_low = interval.bplow.flow-ballooncost_low;
	// ballooning cost for upper parameter solution
	int ssize_up = interval.bpup.ssize;
	double ballooncost_up = ssize_up*interval.bpup.para;
	double entropy_up = interval.bpup.flow-ballooncost_up;
	//out(ssize_up - ssize_low);
	double newpara = (entropy_low - entropy_up) / (double)(ssize_up - ssize_low);
	return newpara;
}

Table2D<double> gethistdist(const Table2D<RGB> & rgbimg,int binsize, Table2D<bool> ROI){
	int img_w = rgbimg.getWidth();
	int img_h = rgbimg.getHeight();
	int r_idx =0, g_idx = 0, b_idx = 0, idx =0;
	int binnum = (int)ceil(256.0/binsize);
	Table2D<vector<double> > hist(binnum,binnum,vector<double>(binnum,0));
	for(unsigned int j=0;j<img_h;j++)
	{
		for(unsigned int i=0;i<img_w;i++)
		{
			if(ROI[i][j]){
				RGB c = rgbimg[i][j];
				r_idx = (int)(c.r/binsize);
				g_idx = (int)(c.g/binsize);
				b_idx = (int)(c.b/binsize);
				hist[r_idx][g_idx][b_idx] = hist[r_idx][g_idx][b_idx]+1;
			}
		}
	}
	double hist_sum = countintable(ROI,true);
	// normalize the histogram
	for(r_idx=0;r_idx<binnum;r_idx++){
		for(g_idx=0;g_idx<binnum;g_idx++){
			for(b_idx=0;b_idx<binnum;b_idx++){
				hist[r_idx][g_idx][b_idx] = hist[r_idx][g_idx][b_idx] / hist_sum;
			}
		}
	}

	// return probablity map
	Table2D<double> returnv(img_w,img_h,0);
	for(unsigned int j=0;j<img_h;j++)
	{
		for(unsigned int i=0;i<img_w;i++)
		{
			RGB c = rgbimg[i][j];
			r_idx = (int)(c.r/binsize);
			g_idx = (int)(c.g/binsize);
			b_idx = (int)(c.b/binsize);
			returnv[i][j] = hist[r_idx][g_idx][b_idx];
		}
	}
	return returnv;
}

Table2D<double> getxyhistdist(const Table2D<RGB> & rgbimg,int binsize, int xybinsize, Table2D<bool> ROI){
	int img_w = rgbimg.getWidth();
	int img_h = rgbimg.getHeight();
	int r_idx =0, g_idx = 0, b_idx = 0, x_idx=0, y_idx=0, idx =0;
	int binnum = (int)ceil(256.0/binsize);
	int x_binnum = (int)ceil((double)img_w/xybinsize);
	int y_binnum = (int)ceil((double)img_h/xybinsize);
	Table2D<vector<Table2D<double> > > hist(binnum,binnum,vector<Table2D<double> >(binnum,Table2D<double>(x_binnum,y_binnum,0)));
	for(unsigned int j=0;j<img_h;j++)
	{
		for(unsigned int i=0;i<img_w;i++)
		{
			if(ROI[i][j]){
				RGB c = rgbimg[i][j];
				r_idx = (int)(c.r/binsize);
				g_idx = (int)(c.g/binsize);
				b_idx = (int)(c.b/binsize);
				x_idx = (int)(i/xybinsize);
				y_idx = (int)(j/xybinsize);
				hist[r_idx][g_idx][b_idx][x_idx][y_idx] = hist[r_idx][g_idx][b_idx][x_idx][y_idx]+1;
			}
		}
	}
	double hist_sum = countintable(ROI,true);
	// normalize the histogram
	for(r_idx=0;r_idx<binnum;r_idx++){
		for(g_idx=0;g_idx<binnum;g_idx++){
			for(b_idx=0;b_idx<binnum;b_idx++){
				for(x_idx=0;x_idx<x_binnum;x_idx++){
					for(y_idx=0;y_idx<y_binnum;y_idx++){
						hist[r_idx][g_idx][b_idx][x_idx][y_idx] = hist[r_idx][g_idx][b_idx][x_idx][y_idx] / hist_sum;
					}
				}
			}
		}
	}

	// return probablity map
	Table2D<double> returnv(img_w,img_h,0);
	for(unsigned int j=0;j<img_h;j++)
	{
		for(unsigned int i=0;i<img_w;i++)
		{
			RGB c = rgbimg[i][j];
			r_idx = (int)(c.r/binsize);
			g_idx = (int)(c.g/binsize);
			b_idx = (int)(c.b/binsize);
			x_idx = (int)(i/xybinsize);
			y_idx = (int)(j/xybinsize);
			returnv[i][j] = hist[r_idx][g_idx][b_idx][x_idx][y_idx];
		}
	}
	return returnv;
}
#endif
