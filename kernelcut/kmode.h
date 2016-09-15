#ifndef _KMODE_H_
#define _KMODE_H_
Vect3D meanshiftinit(const vector<Vect3D> & data, double * kernelwidth, int num_init);
inline Vect3D samplepoint(const vector<Vect3D> & data){
    int N = data.size();
    return data[rand()%N];
}

inline double gaussiankernel(Vect3D v1, Vect3D v2, double * sigmas){
    Vect3D v3 = v1 - v2;
    double diff = pow(v3.x/sigmas[0],2.0) + pow(v3.y/sigmas[1],2.0) + pow(v3.z/sigmas[2],2.0);
    if(diff/2.0>300) {return 0;}
    return exp(-diff/2.0);
}

double meanshiftenergy(const vector<Vect3D> & data, Vect3D mu, double * kernelwidth){
    double kernelsum = 0;
    int N = data.size();
    for(int n=0;n<N;n++){
        double k = gaussiankernel(data[n],mu,kernelwidth);
        kernelsum += k;
    }
    return kernelsum;
}

Vect3D meanshiftonce(const vector<Vect3D> & data, Vect3D mu, double * kernelwidth){
    double kernelsum = 0;
    Vect3D wkernelsum;
    int N = data.size();
    for(int n=0;n<N;n++){
        double k = gaussiankernel(data[n],mu,kernelwidth);
        kernelsum += k;
        wkernelsum.x += k*(data[n].x);
        wkernelsum.y += k*(data[n].y);
        wkernelsum.z += k*(data[n].z);
    }
    return Vect3D(wkernelsum / (kernelsum+1e-10));
}

Vect3D meanshiftconverge(const vector<Vect3D> & data, Vect3D mu, double * kernelwidth, double delta){
    while(1){
        Vect3D new_mu = meanshiftonce(data,mu,kernelwidth);
        if((new_mu-mu).norm()>delta) // converged?
            mu = new_mu;
        else
            break;
        //outv(meanshiftenergy(data,mu,kernelwidth));
    }
    return mu;
}

Vect3D meanshiftinit(const vector<Vect3D> & data, double * kernelwidth, int num_init){
    double best_e = -1e+20;
    Vect3D mu;
    for(int i=0;i<num_init;i++){
        Vect3D new_mu = samplepoint(data);
        double e = meanshiftenergy(data,new_mu,kernelwidth);
        if(e > best_e) {best_e = e; mu = new_mu;}
    }
    return mu;
}

double computkmodeeenergy(const Image & image, const Table2D<Vect3D> & floatimg, const Table2D<Label> & labeling, Vect3D obj_mu, Vect3D bkg_mu, double * kernelwidth, double w_smooth)
{
	vector<Vect3D> obj_data = getvectordata(floatimg, getROI(labeling,OBJ));
    vector<Vect3D> bkg_data = getvectordata(floatimg, getROI(labeling,BKG));
    double obj_e = meanshiftenergy(obj_data, obj_mu, kernelwidth);
    double bkg_e = meanshiftenergy(bkg_data, bkg_mu, kernelwidth);
    double kmode_e = 2*image.img_w*image.img_h - obj_e - bkg_e;
	if(w_smooth <1e-10)
		return kmode_e;
	double smoothenergy = 0;
	// number of neighboring pairs of pixels
	int numNeighbor = image.pointpairs.size();
	// n-link - smoothness term
	int node_id1 =0, node_id2 =0;
	for(int i=0;i<numNeighbor;i++)
	{
		PointPair pp = image.pointpairs[i];
		node_id1 = pp.p1.x+pp.p1.y*image.img_w;
		node_id2 = pp.p2.x+pp.p2.y*image.img_w;
		if((labeling[pp.p1]!=NONE)&&(labeling[pp.p2]!=NONE))
		{
			if(labeling[pp.p1]!=labeling[pp.p2])
			{
				double v = image.smoothnesscosts[i];
				smoothenergy += v;
			}
		}
	}
	return kmode_e + smoothenergy * w_smooth;
}

Table2D<Label> kmodesegmentation(const Image & image, Table2D<Vect3D> floatimg, double * kernelwidth, double w_smooth, Table2D<Label> initlabeling, Table2D<Label> hardconstraints)
{
    int img_w = image.img_w;
    int img_h = image.img_h;
    Table2D<Label> solution = initlabeling;
    // initialize mu
    vector<Vect3D> obj_data = getvectordata(floatimg, getROI(solution,OBJ));
    vector<Vect3D> bkg_data = getvectordata(floatimg, getROI(solution,BKG));
    Vect3D obj_mu = meanshiftinit(obj_data,kernelwidth,100);
    Vect3D bkg_mu = meanshiftinit(bkg_data,kernelwidth,100);
    obj_mu = meanshiftconverge(obj_data,obj_mu,kernelwidth,0.1);
    bkg_mu = meanshiftconverge(bkg_data,bkg_mu,kernelwidth,0.1);
    
    Table2D<double> capsource(img_w,img_h,0),capsink(img_w,img_h,0);
    double current_e = 1e+20;
    int itr_num = 0;
    while(1){
        // update segmentation
        GraphType * g = new GraphType(img_w*img_h, 4*img_w*img_h); 
        g->add_node(img_w*img_h);    // adding nodes
	    if(w_smooth>1e-10)
            addsmoothnessterm(g,image,w_smooth,Table2D<bool>(img_w,img_h,true),false);
        for(int i=0;i<img_w;i++){
            for(int j=0;j<img_h;j++){
                capsource[i][j] = 2 - 2 * gaussiankernel(floatimg[i][j],bkg_mu,kernelwidth);
                capsink[i][j] = 2- 2 * gaussiankernel(floatimg[i][j],obj_mu,kernelwidth);
                //capsource[i][j] = (floatimg[i][j]-bkg_mu).norm();
                //capsink[i][j] = (floatimg[i][j]-obj_mu).norm();
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
		double energy = computkmodeeenergy(image, floatimg, m_labeling, obj_mu, bkg_mu, kernelwidth, w_smooth);
		if((current_e - energy < 0.1) || (itr_num > 20))
		    break;
		else{
		    current_e = energy;
		    solution = m_labeling;
		    outv(current_e);
		}
        // update mu
        vector<Vect3D> obj_data = getvectordata(floatimg, getROI(solution,OBJ));
        vector<Vect3D> bkg_data = getvectordata(floatimg, getROI(solution,BKG));
        obj_mu = meanshiftconverge(obj_data,obj_mu,kernelwidth,0.1);
        bkg_mu = meanshiftconverge(bkg_data,bkg_mu,kernelwidth,0.1);
        // sampled mu
        Vect3D obj_mu2 = meanshiftinit(obj_data,kernelwidth,50);
        Vect3D bkg_mu2 = meanshiftinit(bkg_data,kernelwidth,50);
        obj_mu2 = meanshiftconverge(obj_data,obj_mu2,kernelwidth,0.1);
        bkg_mu2 = meanshiftconverge(bkg_data,bkg_mu2,kernelwidth,0.1);
        if(meanshiftenergy(obj_data,obj_mu2,kernelwidth)<meanshiftenergy(obj_data,obj_mu,kernelwidth))
        {
            printf("sampled mu\n");
            obj_mu = obj_mu2;
        }
        if(meanshiftenergy(bkg_data,bkg_mu2,kernelwidth)<meanshiftenergy(bkg_data,bkg_mu,kernelwidth))
        {
            printf("sampled mu\n");
            bkg_mu = bkg_mu2;
        }
    }
    return solution;
}

#endif
