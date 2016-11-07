#ifndef _NCUTKNN_H_
#define _NCUTKNN_H_
#include "knn.h"
#include "PPBCncut.h"
Table2D<Label> ncutknnbinarysegmentation(const Image & image, Table2D<int> & knntable, double w_smooth, Table2D<Label> initlabeling, Table2D<Label> hardconstraints, bool PBO = false);

Table2D<Label> ncutknnbinarysegmentation(const Image & image, Table2D<int> & knntable, double w_smooth, Table2D<Label> initlabeling, Table2D<Label> hardconstraints, bool PBO)
{
    int KNN_K = knntable.getHeight();
    int img_w = image.img_w;
    int img_h = image.img_h;
    
    //perturbconnection(knntable,img_w,img_h);
    
    Table2D<double> w_data_OBJ = zeroonekernel(knntable, getROI(initlabeling,OBJ),KNN_K);
    Table2D<double> w_data_BKG = zeroonekernel(knntable, getROI(initlabeling,BKG),KNN_K);
    for(int i=0;i<img_w;i++){
        for(int j=0;j<img_h;j++){
            if(initlabeling[i][j]!=NONE)
                continue;
            if(w_data_OBJ[i][j]>w_data_BKG[i][j]) 
                initlabeling[i][j]=OBJ;
			else
			    initlabeling[i][j]=BKG;
			if(hardconstraints[i][j] == OBJ) initlabeling[i][j]=OBJ;
			if(hardconstraints[i][j] == BKG) initlabeling[i][j]=BKG;
        }
    }
		
    PPBCncut ppbc = PPBCncut(image, 10000.0, w_smooth, & knntable, KNN_K,hardconstraints);
	ppbc.setpara(-0.5,0.5,0.008,20,true);
	int bo_itr = 0,pbo_itr = 0;
	double ppbc_energy = 1e+10;
	
	Table2D<Label> solution = initlabeling;
	//return solution;

	while(++bo_itr){
		ppbc.setinitlabeling(solution);
		BreakPoint best_bp = ppbc.explorepara(0);
		best_bp.original_e = ppbc.computeenergy(best_bp.solution);
		//best_bp.print();
		//return;
		if((best_bp.original_e<ppbc_energy-0.1)&&(best_bp.ssize!=0)&&(best_bp.ssize!=ppbc.getROISize()))
		{
			solution = best_bp.solution;
			ppbc_energy = best_bp.original_e;
			//outv(ppbc_energy);
			//outv(best_bp.ssize);
			//savebinarylabeling(rgbimg,init_labeling,(destdir+string(imgname)+"_ppbc_"+pch::to_string(ppbc_itr)+string(".bmp")).c_str());
		}
		else
			break;
		if(bo_itr>25)
			break;
	}
    
    if(PBO==false)
        return solution;
        
	while(++pbo_itr){
		break;
		if(pbo_itr>8)
			break;
		ppbc.setinitlabeling(solution);
		ppbc.explore();
		//ppbc.savesolutions(image,string(destdir));
		BreakPoint best_bp = ppbc.SelectBestBP();
		//best_bp.print();
		if((best_bp.original_e<ppbc_energy-0.1)&&(best_bp.ssize!=0)&&(best_bp.ssize!=img_w*img_h))
		{
			solution = best_bp.solution;
			ppbc_energy = best_bp.original_e;
			//outv(ppbc_energy);
		}
		else
			break;
	}
	return solution;
}
#endif
