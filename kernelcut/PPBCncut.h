#ifndef _PPBCNCUT_H_
#define _PPBCNCUT_H_
// PPBC for kernel segmentation (Average association on KNN graph)
#include "PPBCBase.h"

class PPBCncut: public PPBCBase{
public:
	PPBCncut(const Image & image_, double w_ncut_, double w_smooth_, Table2D<int> * knntable_, int KNN_K_,const Table2D<Label> & hardconstraints_)
	:image(image_),w_ncut(w_ncut_),w_smooth(w_smooth_){
		knntable = knntable_;
		KNN_K = KNN_K_;
		img_w = image.img_w;
		img_h = image.img_h;
		hardconstraints = hardconstraints_;
		weights =  zeroonekernel(*knntable, Table2D<bool>(img_w,img_h,true), KNN_K);
		offset = -2;
		for(int i=0;i<img_w;i++){
			for(int j=0;j<img_h;j++){
				offset += 1 / weights[i][j];
			}
		}
		offset = offset * w_ncut;
		showflag = false;
	};
	virtual double computeenergy(const Table2D<Label> & labeling);
	double computeenergymulti(Table2D<int> multilabeling);
	virtual void updatemodel(); // collect statistics about current labeling
	virtual void * parabasegraph(double para_, UnknownRegion * unknownregion_p,
		double * flowoffset_p, vector<Point> * node_corr_p);
	virtual double getnewpara(ParaInterval paraInterval);

	bool showflag;
private:
	Image image;
	double w_ncut;
	double w_smooth;
	Table2D<int> * knntable;
	int KNN_K;
	Table2D<Label> hardconstraints;
	// unary auxiliary function
	Table2D<double> capsource;
	Table2D<double> capsink;

	Table2D<double> weights; // normailized cut is weighted kernel k-means
	double offset; // weighted kernel kmeans energy - normalized cut energy (scaled by w_ncut)
};

double PPBCncut::computeenergy(const Table2D<Label> & labeling)
{
	double obj_weights = weights.sum(getROI(labeling,OBJ));
	double bkg_weights = weights.sum(getROI(labeling,BKG));

	int cuts[3]={0,0,0}; // for NONE=0, OBJ=1, BKG=2
	Label l_p,l_q;
	for(int i=0;i<img_w;i++){
	    for(int j=0;j<img_h;j++){
			int p_idx = j+i*img_h;
			l_p = labeling[i][j];
		    for(int k=0;k<KNN_K;k++){
			    int q_idx = (*knntable)[p_idx][k];
				if(labeling[q_idx/img_h][q_idx%img_h]!=l_p){
				    cuts[l_p]++;
				    cuts[labeling[q_idx/img_h][q_idx%img_h]]++;
				}
			}
		}
	}
	double ncut_e = (double)cuts[1] / (double)(obj_weights+1e-10) + (double)cuts[2] / (double)(bkg_weights+1e-10);
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
	if(showflag) printf("ncut energy: %.2f + %.2f = %.2f\n",ncut_e*w_ncut,smoothenergy * w_smooth,ncut_e*w_ncut + smoothenergy * w_smooth);
	return ncut_e*w_ncut + smoothenergy * w_smooth;
}

double PPBCncut::computeenergymulti(Table2D<int> multilabeling)
{
    int numCluster = multilabeling.getMax()+1;
    vector<double> segment_weights(numCluster);
    for(int i=0;i<numCluster;i++)
        segment_weights[i] = weights.sum(getROI(multilabeling,i));

	vector<int>cuts(numCluster,0);
	int l_p,l_q;
	for(int i=0;i<img_w;i++){
	    for(int j=0;j<img_h;j++){
			int p_idx = j+i*img_h;
			l_p = multilabeling[i][j];
		    for(int k=0;k<KNN_K;k++){
			    int q_idx = (*knntable)[p_idx][k];
				if((multilabeling[q_idx/img_h][q_idx%img_h]!=l_p)){
				    cuts[l_p]++;
				    cuts[multilabeling[q_idx/img_h][q_idx%img_h]]++;
				}
			}
		}
	}
	double ncut_e = 0;
	for(int i=0;i<numCluster;i++)
	    ncut_e += (double)cuts[i] / (double)(segment_weights[i]+1e-10);
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
	    if(multilabeling[pp.p1]!=multilabeling[pp.p2])
		{
			double v = image.smoothnesscosts[i];
			smoothenergy += v;
		}
	}
	if(showflag) printf("ncut energy: %.2f + %.2f = %.2f\n",ncut_e*w_ncut,smoothenergy * w_smooth,ncut_e*w_ncut + smoothenergy * w_smooth);
	return ncut_e*w_ncut + smoothenergy * w_smooth;
}

void PPBCncut::updatemodel()
{
	capsource = Table2D<double>(img_w,img_h,0), capsink = Table2D<double>(img_w,img_h,0);

	Table2D<double> w_data_OBJ;
	w_data_OBJ = zeroonekernel(*knntable, getROI(initlabeling,OBJ),KNN_K);
	double obj_size  = weights.sum(getROI(initlabeling,OBJ));
	double obj_sum = w_data_OBJ.sum(getROI(initlabeling, OBJ));
	for(int i=0;i<img_w;i++){
		for(int j=0;j<img_h;j++){
			double wp = weights[i][j];
			capsink[i][j] = (1.0 / wp - 2*w_data_OBJ[i][j]/obj_size + obj_sum / obj_size / obj_size * wp)*w_ncut ;
		}
	}

	Table2D<double> w_data_BKG;
		w_data_BKG = zeroonekernel(*knntable, getROI(initlabeling,BKG),KNN_K);
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

void * PPBCncut::parabasegraph(double para_, UnknownRegion * unknownregion_p,
	double * flowoffset_p, vector<Point> * node_corr_p)
{
	if(unknownregion_p==NULL){
	GraphType * g = new GraphType(img_w*img_h, 4*img_w*img_h); 
	g->add_node(img_w*img_h);    // adding nodes
	// add smoothness term
	if(w_smooth> 1e-10)
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


double PPBCncut::getnewpara(ParaInterval interval)
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
