#ifndef _NCUTKNNMULTI_H_
#define _NCUTKNNMULTI_H_
#include "knn.h"
Table2D<Label> getswaplabeling(Table2D<int> labeling, int alpha, int beta);
void applyswaplabeling(Table2D<int> & labeling, int alpha, int beta, Table2D<Label> swaplabeling);
bool ncutswap(const Image & image, Table2D<int> & solution,Table2D<int> & knntable,int KNN_K, double w_smooth, Table2D<int> hardconstraints, int alpha, int beta);
Table2D<int> ncutknnmultisegmentation(const Image & image, Table2D<int> & knntable, double w_smooth, Table2D<int> initlabeling, int numColor){
    int KNN_K = knntable.getHeight();
    int img_w = image.img_w;
    int img_h = image.img_h;
    if(KNN_K==0){
        printf("KNN table empty!\n");
        exit(-1);
    }
    Table2D<int> solution = initlabeling;
    Table2D<int> hardconstraints = initlabeling;
    hardconstraints.reset(-1);
    
    // initial solution is not complete
    if(countintable(solution,-1)>0){
    vector<Table2D<double> > w_data(numColor);
    for(int c=0;c<numColor;c++)
        w_data[c] = zeroonekernel(knntable, getROI(initlabeling,c),KNN_K);
    for(int i=0;i<img_w;i++){
        for(int j=0;j<img_h;j++){
            double maxsofar = -1e+6;
            for(int c=0;c<numColor;c++){
                if(w_data[c][i][j]>maxsofar){
                    maxsofar = w_data[c][i][j];
                    solution[i][j] = c;
                }
            }
			if(hardconstraints[i][j] >=0) solution[i][j]=initlabeling[i][j];
        }
    }
    }
    
    while(1){
        bool swapsuccess = false;
        for(int alpha=0;alpha<numColor;alpha++){
            for(int beta=alpha+1;beta<numColor;beta++){
                if(ncutswap(image, solution,knntable,KNN_K,w_smooth, hardconstraints, alpha,beta))
                    swapsuccess = true;
            }
        }
        if(swapsuccess==false)
            break;
    }
    return solution;
}

Table2D<Label> getswaplabeling(Table2D<int> labeling, int alpha, int beta){
    int img_w = labeling.getWidth();
    int img_h = labeling.getHeight();
    Table2D<Label> swaplabeling(img_w,img_h,NONE);
    for(int i=0;i<img_w;i++){
        for(int j=0;j<img_h;j++){
            if(labeling[i][j]==alpha)
                swaplabeling[i][j] = OBJ;
            else if(labeling[i][j]==beta)
                swaplabeling[i][j] = BKG;
        }
    }
    return swaplabeling;
}

void applyswaplabeling(Table2D<int> & labeling, int alpha, int beta, Table2D<Label> swaplabeling){
    int img_w = labeling.getWidth();
    int img_h = labeling.getHeight();
    for(int i=0;i<img_w;i++){
        for(int j=0;j<img_h;j++){
            if(swaplabeling[i][j]==OBJ)
                labeling[i][j] = alpha;
            else if(swaplabeling[i][j]==BKG)
                labeling[i][j] = beta;
        }
    }
}

bool ncutswap(const Image & image, Table2D<int> & solution,Table2D<int> & knntable,int KNN_K, double w_smooth, Table2D<int> hardconstraints, int alpha, int beta){
    bool success = false;
    int img_w = image.img_w;
    int img_h = image.img_h;
    Table2D<Label> swaplabeling = getswaplabeling(solution,alpha,beta);
    Table2D<Label> hardconstraints_binary = getswaplabeling(hardconstraints,alpha,beta);
    PPBCncut ppbc = PPBCncut(image, 10000.0, w_smooth, & knntable, KNN_K,hardconstraints_binary);
	ppbc.setpara(-0.5,0.5,0.008,20,true);
	printf("Swap label %d and %d\n",alpha,beta);
	
	double ppbc_energy = ppbc.computeenergy(swaplabeling);
	double multilabel_energy = ppbc.computeenergymulti(solution);
	outv(multilabel_energy);
	//printf("Init binary and multilabel energy: %.2f %.2f, diff %.2f\n",ppbc_energy,multilabel_energy,multilabel_energy-ppbc_energy);
	
	
	int bo_itr = 0;
	while(++bo_itr){
		ppbc.setinitlabeling(swaplabeling);
		BreakPoint best_bp = ppbc.explorepara(0);
		best_bp.original_e = ppbc.computeenergy(best_bp.solution);
		//best_bp.print();
		//return;
		if((best_bp.original_e<ppbc_energy-0.1)&&(best_bp.ssize!=0)&&(best_bp.ssize!=ppbc.getROISize()))
		{
			swaplabeling = best_bp.solution;
			ppbc_energy = best_bp.original_e;
			//outv(ppbc_energy);
			success = true;
		}
		else
			break;
		if(bo_itr>25)
			break;
	}
	if(success){
	     applyswaplabeling(solution, alpha, beta, swaplabeling);
	     multilabel_energy = ppbc.computeenergymulti(solution);
	     outv(multilabel_energy);
	     //printf("After Swap binary and multilabel energy: %.2f %.2f, diff %.2f\n",ppbc_energy,multilabel_energy,multilabel_energy-ppbc_energy);
	}
	return success;
}
#endif
