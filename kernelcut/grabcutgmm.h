#ifndef _GRABCUTGMM_H_
#define _GRABCUTGMM_H_
#include "gmm.h"

Table2D<Label> grabcutgmmsegmentation(const Image & image, Table2D<Vect3D> floatimg, int num_comps, double w_smooth, Table2D<Label> initlabeling, Table2D<Label> hardconstraints){
    int img_w = image.img_w;
    int img_h = image.img_h;
    Table2D<Label> solution = initlabeling;

    GMM objgmm(image.img,getROI(solution,OBJ),num_comps);
	GMM bkggmm(image.img,getROI(solution,BKG),num_comps);
    Table2D<double> capsource(img_w,img_h,0),capsink(img_w,img_h,0);
    double current_e = 1e+20;
    int itr_num = 0;
    while(1){
        // update segmentation
        GraphType * g = new GraphType(img_w*img_h, 4*img_w*img_h); 
        g->add_node(img_w*img_h);    // adding nodes
	    if(w_smooth>1e-10)
            addsmoothnessterm(g,image,w_smooth,Table2D<bool>(img_w,img_h,true),false);
        objgmm.feeddata(image.img,getROI(solution,OBJ));		
        objgmm.optimize(20);
        capsink = objgmm.getProbability(image.img, Table2D<bool>(img_w,img_h,true));

        bkggmm.feeddata(image.img,getROI(solution,BKG));
        bkggmm.optimize(20);
        capsource = bkggmm.getProbability(image.img, Table2D<bool>(img_w,img_h,true));
        for(int i=0;i<img_w;i++){
            for(int j=0;j<img_h;j++){
                if(hardconstraints[i][j]==OBJ) capsource[i][j] = 1e+20;
                if(hardconstraints[i][j]==BKG) capsink[i][j] = 1e+20;
                g->add_tweights(i+j*img_w,capsource[i][j],capsink[i][j]);
            }
        }
        double flow = g -> maxflow();
		//outv(flow);
		
		Table2D<Label> m_labeling(img_w,img_h);
		if(!getlabeling(g,m_labeling))
		{
			cout<<"trivial solution!"<<endl;
			delete g;
			break; // trivial solution
		}
		delete g;

		if(solution==m_labeling)
		{
			cout<<"labeling converged!"<<endl;
			break;
		}
		double energy = flow;
		if((current_e - energy < 1.0) || (itr_num > 20))
		    break;
		else{
		    current_e = energy;
		    solution = m_labeling;
		    outv(current_e);
		}
    }
    return solution;
}
#endif
