#ifndef _GRABCUT_H_
#define _GRABCUT_H_
#include "ppbcgrabcut.h"
Table2D<Label> grabcutsegmentation(const Image & image, double binsize, int xybinsize, double w_smooth, Table2D<Label> initlabeling, Table2D<Label> hardconstraints,bool PBO = false);

Table2D<Label> grabcutsegmentation(const Image & image, double binsize, int xybinsize, double w_smooth, Table2D<Label> initlabeling, Table2D<Label> hardconstraints, bool PBO)
{
    PPBCgrabcut ppbc = PPBCgrabcut(image, binsize, xybinsize, 1.0, w_smooth, hardconstraints);
	ppbc.setpara(-4,4,0.01,20,false);
	
	Table2D<Label> solution = initlabeling;
	
	int bo_itr = 0,pbo_itr = 0;
	double ppbc_energy = 1e+20;

	while(++bo_itr){
		ppbc.setinitlabeling(solution);
		BreakPoint best_bp = ppbc.explorepara(0);
		best_bp.original_e = ppbc.computeenergy(best_bp.solution);
		//best_bp.print();
		//return;
		if((best_bp.original_e<ppbc_energy-1)&&(best_bp.ssize!=0)&&(best_bp.ssize!=image.img_w*image.img_h))
		{
			solution = best_bp.solution;
			ppbc_energy = best_bp.original_e;
			//outv(ppbc_energy);
			//savebinarylabeling(rgbimg,init_labeling,(destdir+string(imgname)+"_ppbc_"+pch::to_string(ppbc_itr)+string(".bmp")).c_str());
		}
		else
			break;
		if(bo_itr>20)
			break;
	}
	if(false == PBO) return solution;
	
	while(++pbo_itr){
		//break;
		if(pbo_itr>8)
			break;
		ppbc.setinitlabeling(solution);
		ppbc.explore();
		//ppbc.savesolutions(image,string(destdir));
		BreakPoint best_bp = ppbc.SelectBestBP();
		//best_bp.print();
		if((best_bp.original_e<ppbc_energy-0.1)&&(best_bp.ssize!=0)&&(best_bp.ssize!=image.img_w*image.img_h))
		{
			solution = best_bp.solution;
			ppbc_energy = best_bp.original_e;
			outv(ppbc_energy);
			//savebinarylabeling(rgbimg,init_labeling,(destdir+string(imgname)+"_pbo_"+pch::to_string(pbo_itr)+string(".bmp")).c_str());
		}
		else
			break;
	}
	return solution;
}
#endif
