#ifndef _KMEANS_H_
#define _KMEANS_H_

Vect3D getmean(const Table2D<Vect3D> & floatimg, Table2D<bool> ROI);
double getstd(const Table2D<Vect3D> & floatimg, Table2D<bool> ROI);
Table2D<Label> kmeanssegmentation(const Image & image, Table2D<Vect3D> floatimg, double w_smooth, Table2D<Label> initlabeling, Table2D<Label> hardconstraints, bool elliptical = false);

Vect3D getmean(const Table2D<Vect3D> & floatimg, Table2D<bool> ROI)
{
    Vect3D mean_v(0,0,0);
    int img_w = floatimg.getWidth();
    int img_h = floatimg.getHeight();
    for(int i=0;i<img_w;i++){
        for(int j=0;j<img_h;j++){
            if(ROI[i][j])
                mean_v = mean_v + floatimg[i][j];
        }
    }
    return mean_v / (countintable(ROI,true));
}

double getstd(const Table2D<Vect3D> & floatimg, Table2D<bool> ROI)
{
    Vect3D mean_v = getmean(floatimg, ROI);
    int img_w = floatimg.getWidth();
    int img_h = floatimg.getHeight();
    double mystd = 0;
    for(int i=0;i<img_w;i++){
        for(int j=0;j<img_h;j++){
            if(ROI[i][j])
                mystd += (mean_v - floatimg[i][j]).norm();
        }
    }
    return mystd / (countintable(ROI,true));
}

Table2D<Label> kmeanssegmentation(const Image & image, Table2D<Vect3D> floatimg, double w_smooth, Table2D<Label> initlabeling, Table2D<Label> hardconstraints, bool elliptical)
{
    int img_w = image.img_w;
    int img_h = image.img_h;
    Vect3D obj_mean, bkg_mean;
    double obj_std =1, bkg_std = 1;
    Table2D<Label> solution = initlabeling;
    while(1){
        // update mean
        obj_mean = getmean(floatimg, getROI(solution,OBJ));
        bkg_mean = getmean(floatimg, getROI(solution,BKG));
        // elliptical kmeans?
        if(elliptical){
            obj_std = getstd(floatimg, getROI(solution,OBJ));
            bkg_std = getstd(floatimg, getROI(solution,BKG));
        }
        outv(obj_std);
        outv(bkg_std);
        // update labeling
        Table2D<Label> newsolution(img_w,img_h);
        for(int i=0;i<img_w;i++){
            for(int j=0;j<img_h;j++){
                double obj_logprob = 0.5* pow( (floatimg[i][j]-obj_mean).norm() / obj_std, 2.0 ) + 3 * log(obj_std);
                double bkg_logprob = 0.5* pow( (floatimg[i][j]-bkg_mean).norm() / bkg_std, 2.0 ) + 3 * log(bkg_std);
                if( obj_logprob < bkg_logprob )
                    newsolution[i][j] = OBJ;
                else
                    newsolution[i][j] = BKG;
                if(hardconstraints[i][j]==OBJ)
                    newsolution[i][j] = OBJ;
                else if(hardconstraints[i][j]==BKG)
                    newsolution[i][j] = BKG;
            }
        }
        // converged?
        if(newsolution == solution)
            return solution;
        solution = newsolution;
    }
}
#endif
