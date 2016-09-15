#ifndef _HDGRABCUT_H_
#define _HDGRABCUT_H_
#include <vector>

class HDPixel{
public:
    int x, y;
    int ndim;
    HDPixel(){ndim = 0;}
    HDPixel(int ndim_)
    {
        ndim = ndim_;
        features = vector<int>(ndim,0);
    }
    HDPixel(const HDPixel & p){
        x = p.x;
        y = p.y;
        ndim = p.ndim;
        if(ndim!=0){
            features = vector<int>(ndim,0);
            for(int d=0;d<ndim;d++)
                features[d] = p.features[d];
        }
    }
    ~HDPixel()
    {
        //if(features!=NULL)
            //delete [] features;
    }
    int & operator [] (int n) {return features[n];}
    bool operator == (HDPixel & p){
        for(int d=0;d<ndim;d++)
        {
            if(features[d]!=p[d])
                return false;
        }
        return true;
    }
    void print(){
        //cout<<"x and y: "<<x<<' '<<y<<endl;
        cout<<"features: ";
        for(int d=0;d<ndim;d++)
            cout<<features[d]<<' ';
        cout<<endl;
    }
private:
    vector<int> features;      
};

Table2D<int> sparsehist(Table2D<RGB> image, int binsize, int patchsize)
{
    assert(patchsize == 1 || patchsize == 3 || patchsize == 5 || patchsize == 9);
    int ndim = patchsize * 3;
    int img_w = image.getWidth();
    int img_h = image.getHeight();
    vector<HDPixel> pixels(img_w*img_h);
    // padded image
    Table2D<RGB> p_image(img_w+2,img_h+2);
    for(int i=0;i<img_w;i++){
        for(int j=0;j<img_h;j++){
            p_image[i+1][j+1] = image[i][j];
        }
    }
    for(int i=0;i<img_w;i++) p_image[i][0] = image[i][1];
    for(int i=0;i<img_w;i++) p_image[i][img_h] = image[i][img_h-2];
    for(int j=0;j<img_h;j++) p_image[0][j] = image[1][j];
    for(int j=0;j<img_h;j++) p_image[img_w][j] = image[img_w-2][j];
    p_image[0][0] = image[1][1];
    p_image[0][img_h+1] = image[1][img_h-2];
    p_image[img_w+1][0] = image[img_w-2][1];
    p_image[img_w+1][img_h+1] = image[img_w-2][img_h-2];
    
    int scale = 10;
    for(int i=0;i<img_w;i++){
        for(int j=0;j<img_h;j++){
            HDPixel p(ndim);
            p.x = i;
            p.y = j;
            if(patchsize==1){
                p[0] = (int)(image[i][j].r)/binsize;
                p[1] = (int)(image[i][j].g)/binsize;
                p[2] = (int)(image[i][j].b)/binsize;
            }
            if(patchsize==5){
                int idx = 0;
                for(int x_shift = -1; x_shift <=1; x_shift++){
                    for(int y_shift = -1; y_shift <=1; y_shift++){
                        if(abs(x_shift) + abs(y_shift) < 2)
                        {
                            int a;
                            (x_shift==0 && y_shift == 0)?(a=1):(a=scale);
                            p[idx*3+0] = (int)(p_image[i+x_shift+1][j+y_shift+1].r)/(binsize*a);
                            p[idx*3+1] = (int)(p_image[i+x_shift+1][j+y_shift+1].g)/(binsize*a);
                            p[idx*3+2] = (int)(p_image[i+x_shift+1][j+y_shift+1].b)/(binsize*a);
                            idx++;
                        }
                    }
                }
            }
            if(patchsize==9){
                int idx = 0;
                for(int x_shift = -1; x_shift <=1; x_shift++){
                    for(int y_shift = -1; y_shift <=1; y_shift++){
                        int a;
                        (x_shift==0 && y_shift == 0)?(a=1):(a=scale);
                        p[idx*3+0] = (int)(p_image[i+x_shift+1][j+y_shift+1].r)/(binsize*a);
                        p[idx*3+1] = (int)(p_image[i+x_shift+1][j+y_shift+1].g)/(binsize*a);
                        p[idx*3+2] = (int)(p_image[i+x_shift+1][j+y_shift+1].b)/(binsize*a);
                        idx++;
                    }
                }
            }
            if(patchsize==3){
                int idx = 0;
                for(int x_shift = -1; x_shift <=1; x_shift++){
                    for(int y_shift = 0; y_shift <=0; y_shift++){
                        p[idx*3+0] = (int)(p_image[i+x_shift+1][j+y_shift+1].r)/binsize;
                        p[idx*3+1] = (int)(p_image[i+x_shift+1][j+y_shift+1].g)/binsize;
                        p[idx*3+2] = (int)(p_image[i+x_shift+1][j+y_shift+1].b)/binsize;
                        idx++;
                    }
                }
            }
            pixels[i+j*img_w] = p;
        }
    }
    
    
    // radix sort of data
    for(int d=0;d<ndim;d++){
        // we need buckets
        int num_buckets = ceil( 256.0 / binsize);
        vector<vector<HDPixel> > buckets(num_buckets);
        for(int i=0;i<img_w*img_h;i++){
            buckets[pixels[i][d]].push_back(pixels[i]);
        }
        pixels.clear();
        for(int i=0;i<num_buckets;i++)
            pixels.insert(pixels.end(),buckets[i].begin(),buckets[i].end());
    }
    Table2D<int> labels(img_w,img_h,0);
    
    int index = 0;
    labels[pixels[0].x][pixels[0].y] = index;
    for(int i=1;i<img_w*img_h;i++){
        if(pixels[i] == pixels[i-1])
             labels[pixels[i].x][pixels[i].y] = index;
        else
            labels[pixels[i].x][pixels[i].y] = (++index);
    }
    return labels;
}

Table2D<double> gethistdist(Table2D<int> & compactlabels, int num_bins, Table2D<bool> ROI);

Table2D<Label> patchgrabcutsegmentation(const Image & image, double binsize, int patchsize, double w_smooth, Table2D<Label> initlabeling, Table2D<Label> hardconstraints)
{
    int img_w = image.img_w;
    int img_h = image.img_h;
    Table2D<int> compactlabels = sparsehist(image.img, (int)binsize, patchsize);
    
    //outv(compactlabels.getMax());
    int num_bins = compactlabels.getMax() + 1;
    
    Table2D<Label> solution = initlabeling;
    double current_flow = 1e+10;
    int imgid=0;
    
    Table2D<double> capsource(img_w,img_h), capsink(img_w,img_h);
	while(1)
	{
		GraphType * g;
		g = new GraphType(/*estimated # of nodes*/ img_w*img_h, /*estimated # of edges*/ 4*img_w*img_h); 
	    g->add_node(img_w*img_h);    // adding nodes
	    
	    // construct the graph
	    if(w_smooth > 1e-10)
	        addsmoothnessterm(g, image, w_smooth, Table2D<bool>(img_w,img_h,true),false);
	    
	    // add unary term
	    Table2D<double> obj_prob = gethistdist(compactlabels, num_bins, getROI(solution,OBJ));
	    Table2D<double> bkg_prob = gethistdist(compactlabels, num_bins, getROI(solution,BKG));
	    for(int j=0;j<img_h;j++)
	    {
		    for(int i=0;i<img_w;i++)
		    {
				int n = i+j*img_w;
				capsink[i][j]   = -log(obj_prob[i][j]*(1-1e-10)+1e-10*(1.0/(num_bins)));
			    capsource[i][j] = -log(bkg_prob[i][j]*(1-1e-10)+1e-10*(1.0/(num_bins)));
			    if(BKG==hardconstraints[i][j]) // hard constraints to background
				    capsink[i][j]=1e20;
			    else if(OBJ==hardconstraints[i][j]) // hard constraints to foreground
				    capsource[i][j]=1e20;
				g->add_tweights(n,capsource[i][j],capsink[i][j]);
			}
	    }
	    
	    double flow = g -> maxflow();
		//outv(flow);
		
		Table2D<Label> m_labeling(img_w,img_h);
		if(!getlabeling(g,m_labeling))
		{
			//cout<<"trivial solution!"<<endl;
			delete g;
			break; // trivial solution
		}
		delete g;

		if(solution==m_labeling)
		{
			//cout<<"labeling converged!"<<endl;
			break;
		}
		if(current_flow-flow<0.1)
		{
			//cout<<"energy increase!"<<endl;
			break;
		}
		solution = m_labeling;
		current_flow = flow;
		imgid++;
		
		if(imgid>20){
			//printf("max iterations.\n");
			break;
		}
		
	}
    return solution;
}

Table2D<double> gethistdist(Table2D<int> & compactlabels, int num_bins, Table2D<bool> ROI){
	int img_w = compactlabels.getWidth();
	int img_h = compactlabels.getHeight();
	
	vector<double> hist(num_bins,0);
	
	for(unsigned int j=0;j<img_h;j++)
	{
		for(unsigned int i=0;i<img_w;i++)
		{
			if(ROI[i][j]){
				hist[compactlabels[i][j]] = hist[compactlabels[i][j]] + 1.0;
			}
		}
	}
	
	double hist_sum = countintable(ROI,true);
	// normalize the histogram
	for(int i=0;i<num_bins;i++)
        hist[i] = hist[i] / hist_sum;

	// return probablity map
	Table2D<double> returnv(img_w,img_h,0);
	for(unsigned int j=0;j<img_h;j++)
	{
		for(unsigned int i=0;i<img_w;i++)
		{
			returnv[i][j] = hist[compactlabels[i][j]];
		}
	}
	return returnv;
}
#endif
