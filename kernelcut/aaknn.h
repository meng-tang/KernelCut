#ifndef _AAKNN_H_
#define _AAKNN_H_
#include "knn.h"
#include "PPBCaa.h"
Table2D<Label> aaknnsegmentation(const Image & image, Table2D<int> & knntable, double w_smooth, Table2D<Label> initlabeling, Table2D<Label> hardconstraints, bool PBO = false);

Table2D<Label> aaknnsegmentation(const Image & image, Table2D<int> & knntable, double w_smooth, Table2D<Label> initlabeling, Table2D<Label> hardconstraints, bool PBO)
{
    int KNN_K = knntable.getHeight();
    int img_w = image.img_w;
    int img_h = image.img_h;
    perturbconnection(knntable,img_w,img_h);
    
    PPBCaa ppbc = PPBCaa(image, 5.0, w_smooth, &knntable, KNN_K,hardconstraints);
    
    ppbc.setpara(-0.1,0.1,0.002,20,false);

    int bo_itr = 0,pbo_itr=0;
    double ppbc_energy = 1e+10;
    
    Table2D<Label> solution = initlabeling;
    
    while(++bo_itr){
        ppbc.setinitlabeling(solution);
        BreakPoint best_bp = ppbc.explorepara(0);
        //best_bp.print();
        best_bp.original_e = ppbc.computeenergy(best_bp.solution);
        if((best_bp.original_e<ppbc_energy-0.1)&&(best_bp.ssize!=0)&&(best_bp.ssize!=img_w*img_h)){
            solution = best_bp.solution;
            ppbc_energy = best_bp.original_e;
            //outv(ppbc_energy);
			//savebinarylabeling(rgbimg,init_labeling,(destdir+string(imgname)+"_ppbc_"+pch::to_string(ppbc_itr)+string(".bmp")).c_str());
        }
        else
            break;
        if(bo_itr>20) break;
    }
	
	if(PBO == false)
	    return solution;
	
    while(++pbo_itr){
        ppbc.setinitlabeling(solution);
        ppbc.explore();
        BreakPoint best_bp = ppbc.SelectBestBP();
        //best_bp.print();
        if((best_bp.original_e<ppbc_energy-0.1)&&(best_bp.ssize!=0)&&(best_bp.ssize!=img_w*img_h)){
            solution = best_bp.solution;
            ppbc_energy = best_bp.original_e;
            outv(ppbc_energy);
		    //savebinarylabeling(rgbimg,init_labeling,(destdir+string(imgname)+"_ppbc_"+pch::to_string(ppbc_itr)+string(".bmp")).c_str());
        }
        else
            break;
        if(pbo_itr>8) break;
    }
	return solution;
}
#endif
