#ifndef _FTRKERNEL_H_
#define _FTRKERNEL_H_
// PPBC for kernel segmentation (Average association on KNN graph)
#include "PPBCBase.h"
#include "filterxy.h"

class PPBCncutg: public PPBCBase{
public:
	PPBCncutg(const Image & image_, Table2D<vector<float> > * vector_img_, double w_ncut_, double w_smooth_, double * sigmas_,const Table2D<Label> & hardconstraints_)
	:image(image_),w_ncut(w_ncut_),w_smooth(w_smooth_),hardconstraints(hardconstraints_){
		sigmas[0]=sigmas_[0];sigmas[1]=sigmas_[1];sigmas[2]=sigmas_[2];
		img_w = image.img_w;
		img_h = image.img_h;
		vector_img = vector_img_;
		weights = filterrgb(*vector_img, Table2D<bool>(img_w,img_h,true), sigmas);
		offset = -2;
		for(int i=0;i<img_w;i++){
			for(int j=0;j<img_h;j++){
				offset += 1 / weights[i][j];
			}
		}
		offset = offset * w_ncut;
	};
	virtual double computeenergy(const Table2D<Label> & labeling);
	virtual void updatemodel(); // collect statistics about current labeling
	virtual void * parabasegraph(double para_, UnknownRegion * unknownregion_p,
		double * flowoffset_p, vector<Point> * node_corr_p);
	virtual double getnewpara(ParaInterval paraInterval);
private:
	Image image;
	double w_ncut;
	double w_smooth;
	double sigmas[3];
	Table2D<Label> hardconstraints;
	// unary auxiliary function
	Table2D<double> capsource;
	Table2D<double> capsink;
    Table2D<vector<float> > * vector_img;
	Table2D<double> weights; // normailized cut is weighted kernel k-means
	double offset; // weighted kernel kmeans energy - normalized cut energy (scaled by w_ncut)
};

double PPBCncutg::computeenergy(const Table2D<Label> & labeling)
{
	double obj_weights = weights.sum(getROI(labeling,OBJ));
	double bkg_weights = weights.sum(getROI(labeling,BKG));

	Table2D<double> w_data_OBJ = filterrgb(*vector_img, getROI(labeling, OBJ), sigmas);
	double cuts = w_data_OBJ.sum(getROI(labeling, BKG));
	double ncut_e = (double)cuts / (double)(obj_weights+1e-10) + (double)cuts / (double)(bkg_weights+1e-10);
	if(w_smooth <1e-10)
		return ncut_e*w_ncut;
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
	printf("ncut energy: %.2f + %.2f = %.2f\n",ncut_e*w_ncut,smoothenergy * w_smooth,ncut_e*w_ncut + smoothenergy * w_smooth);
	return ncut_e*w_ncut + smoothenergy * w_smooth;
}

void PPBCncutg::updatemodel()
{
	capsource = Table2D<double>(img_w,img_h,0), capsink = Table2D<double>(img_w,img_h,0);

	Table2D<double> w_data_OBJ;
	w_data_OBJ = filterrgb(*vector_img, getROI(initlabeling, OBJ), sigmas);
	double obj_size  = weights.sum(getROI(initlabeling,OBJ));
	double obj_sum = w_data_OBJ.sum(getROI(initlabeling, OBJ));
	for(int i=0;i<img_w;i++){
		for(int j=0;j<img_h;j++){
			double wp = weights[i][j];
			capsink[i][j] = (1.0 / wp - 2*w_data_OBJ[i][j]/obj_size + obj_sum / obj_size / obj_size * wp)*w_ncut ;
		}
	}

	Table2D<double> w_data_BKG;
	w_data_BKG = filterrgb(*vector_img, getROI(initlabeling, BKG), sigmas);
	double bkg_size = weights.sum(getROI(initlabeling,BKG));
	double bkg_sum = w_data_BKG.sum(getROI(initlabeling, BKG));
	for(int i=0;i<img_w;i++){
		for(int j=0;j<img_h;j++){
			double wp = weights[i][j];
			capsource[i][j] = (1.0 / wp - 2*w_data_BKG[i][j]/bkg_size + bkg_sum / bkg_size / bkg_size * wp)*w_ncut;
		}
	}

	for(int j=0;j<img_h;j++)
	{
		for(int i=0;i<img_w;i++)
		{
			if(hardconstraints[i][j]==BKG) // hard constraints to background
				capsink[i][j]=INFTY;
			else if(hardconstraints[i][j]==OBJ) // hard constraints to foreground
				capsource[i][j]=INFTY;
		}
	}
}

void * PPBCncutg::parabasegraph(double para_, UnknownRegion * unknownregion_p,
	double * flowoffset_p, vector<Point> * node_corr_p)
{
	if(unknownregion_p==NULL){
	GraphType * g = new GraphType(img_w*img_h, 4*img_w*img_h); 
	g->add_node(img_w*img_h);    // adding nodes
	// add smoothness term
	addsmoothnessterm(g, image, w_smooth, ROI);
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

	// monotonic speedup
	// for monotonic mode
	double flowoffset = 0; 
	Table2D<bool> incompactgraph;
	GraphType * compact_g;
	Table2D<int> img_corr;
	vector<PointPair> compactpointpairs;
	vector<double> compactsmoothnesscosts;
	int compactsize = getcompactgraph(image,unknownregion_p->knownlabeling, incompactgraph,
		*node_corr_p,img_corr, compactpointpairs, compactsmoothnesscosts);
	compact_g = new GraphType(compactsize,10*compactsize);
	compact_g->add_node(compactsize);

	// number of neighboring pairs of pixels
	int compactnumNeighbor = compactpointpairs.size();
	// n-link - smoothness term
	for(int i=0;i<compactnumNeighbor;i++)
	{
		PointPair pp = compactpointpairs[i];
		if(incompactgraph[pp.p1]&&incompactgraph[pp.p2])
		{
			double v = w_smooth*compactsmoothnesscosts[i];
			compact_g->add_edge(img_corr[pp.p1],img_corr[pp.p2],v,v);
		}
	}

	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			int n = img_corr[i][j];
			if(incompactgraph[i][j])
			{
				if(unknownregion_p->knownlabeling[i][j]==OBJ)
					compact_g->add_tweights(n,INFTY,capsink[i][j]);
				else if(unknownregion_p->knownlabeling[i][j]==BKG)
					compact_g->add_tweights(n,capsource[i][j],INFTY);
				else if(unknownregion_p->knownlabeling[i][j]==UNKNOWN)
					compact_g->add_tweights(n,capsource[i][j],capsink[i][j]);
			}
			else
			{
				if(unknownregion_p->knownlabeling[i][j]==OBJ)
					flowoffset += capsink[i][j];
				else if(unknownregion_p->knownlabeling[i][j]==BKG)
					flowoffset += capsource[i][j];
			}
		}
	}
	if(flowoffset_p) *flowoffset_p = flowoffset;
	return compact_g;
}


double PPBCncutg::getnewpara(ParaInterval interval)
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
#endif
