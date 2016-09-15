#ifndef _PPBCAA_H_
#define _PPBCAA_H_
// PPBC for kernel segmentation (Average association on KNN graph)
#include "PPBCBase.h"

class PPBCaa: public PPBCBase{
public:
	PPBCaa(const Image & image_, double w_AA_, double w_smooth_, Table2D<int> * knntable_, int KNN_K_,const Table2D<Label> & hardconstraints_)
	:image(image_),w_AA(w_AA_),w_smooth(w_smooth_),hardconstraints(hardconstraints_){
		knntable = knntable_;
		KNN_K = KNN_K_;
		img_w = image.img_w;
		img_h = image.img_h;
	};
	virtual double computeenergy(const Table2D<Label> & labeling);
	virtual void updatemodel(); // collect statistics about current labeling
	virtual void * parabasegraph(double para_, UnknownRegion * unknownregion_p,
		double * flowoffset_p, vector<Point> * node_corr_p);
	virtual double getnewpara(ParaInterval paraInterval);
private:
	Image image;
	double w_AA;
	double w_smooth;
	Table2D<int> * knntable;
	int KNN_K;
	Table2D<Label> hardconstraints;
	// unary auxiliary function
	Table2D<double> capsource;
	Table2D<double> capsink;
};

double PPBCaa::computeenergy(const Table2D<Label> & labeling)
{
	int obj_size = countintable(labeling,OBJ);
	int bkg_size = countintable(labeling,BKG);
	int obj_quad_sum = 0, bkg_quad_sum = 0;
	Label l_p,l_q;
	for(int i=0;i<img_w;i++){
	    for(int j=0;j<img_h;j++){
			int p_idx = j+i*img_h;
			l_p = labeling[i][j];
		    for(int k=0;k<KNN_K;k++){
			    int q_idx = (*knntable)[p_idx][k];
				if(labeling[q_idx/img_h][q_idx%img_h]==l_p){
				    if(l_p==OBJ) obj_quad_sum ++;
					else if(l_p == BKG) bkg_quad_sum ++;
				}
			}
		}
	}
	double AA = (double)obj_quad_sum*2 / (double)(obj_size+1e-10) + (double)bkg_quad_sum*2 / (double)(bkg_size+1e-10);
	if(w_smooth <1e-10)
		return -AA*w_AA;
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
	//printf("AA energy: %.2f + %.2f = %.2f\n",-AA*w_AA,smoothenergy * w_smooth,-AA*w_AA + smoothenergy * w_smooth);
	return -AA*w_AA + smoothenergy * w_smooth;
}

void PPBCaa::updatemodel()
{
	capsource = Table2D<double>(img_w,img_h,0), capsink = Table2D<double>(img_w,img_h,0);

	Table2D<double> w_data_OBJ = zeroonekernel(*knntable, getROI(initlabeling,OBJ),KNN_K);
	double obj_size  = countintable(initlabeling, OBJ);
	double obj_sum = w_data_OBJ.sum(getROI(initlabeling, OBJ));
	for(int i=0;i<img_w;i++){
		for(int j=0;j<img_h;j++){
			capsink[i][j] = (- 2*w_data_OBJ[i][j]/obj_size + obj_sum / obj_size / obj_size)*w_AA ;
		}
	}

	Table2D<double> w_data_BKG = zeroonekernel(*knntable, getROI(initlabeling,BKG),KNN_K);
	double bkg_size = countintable(initlabeling, BKG);
	double bkg_sum = w_data_BKG.sum(getROI(initlabeling, BKG));
	for(int i=0;i<img_w;i++){
		for(int j=0;j<img_h;j++){
			capsource[i][j] = (- 2*w_data_BKG[i][j]/bkg_size + bkg_sum / bkg_size / bkg_size)*w_AA;
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

void * PPBCaa::parabasegraph(double para_, UnknownRegion * unknownregion_p,
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


double PPBCaa::getnewpara(ParaInterval interval)
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
