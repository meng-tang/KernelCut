#ifndef _FILTERING_H_
#define _FILTERING_H_
#include <math.h>
#include <vector>
#include "basicutil.h"

double unique_sum(Table2D<double> & table, const Table2D<RGB> & rgbimg, Table2D<bool> ROI);
Table2D<double> getcolorweights(const Table2D<RGB> & rgbimg, bool smoothed = false);
Table2D<vector<double> > getcolordensities(const Table2D<RGB> & rgbimg);

class Filtering{
private:
	double sigmas[3];
	double sample_rate;
	//double bin_size;
	double bin_sizes[3];
	int bin_nums[3]; // for each dimension
	Table2D<vector<double> > bins;
	Table2D<vector<double> > conv_w; // convolution results
	int img_w;
	int img_h;
public:
	Filtering()
	{}
	Filtering(double sigma_, double sample_rate_ = 2.0)
	{
		sigmas[2] = sigmas[1] = sigmas[0] = sigma_;
		sample_rate = sample_rate_;
		bin_sizes[2] = bin_sizes[1] = bin_sizes[0] = sigma_ / sample_rate;
	}
	Filtering(double sigma1, double sigma2, double sigma3, double sample_rate_ = 2.0)
	{
		sigmas[0] = sigma1;sigmas[1] = sigma2;sigmas[2] = sigma3;
		sample_rate = sample_rate_;
		bin_sizes[0] = sigmas[0] / sample_rate;
		bin_sizes[1] = sigmas[1] / sample_rate;
		bin_sizes[2] = sigmas[2] / sample_rate;
	}
	void splatter(const Table2D<RGB> & rgbimg, Table2D<bool> ROI,Table2D<double> weights = Table2D<double>(1,1,1))
	{
		img_w = rgbimg.getWidth();
		img_h = rgbimg.getHeight();

		bin_nums[0] = floor(255/bin_sizes[0]+0.5)+1;
		bin_nums[1] = floor(255/bin_sizes[1]+0.5)+1;
		bin_nums[2] = floor(255/bin_sizes[2]+0.5)+1;
		bool weighted = true;
		if(1==weights.getWidth()) weighted = false;
		bins = Table2D<vector<double> >(bin_nums[0],bin_nums[1],vector<double>(bin_nums[2],0));
		for(int i=0;i<img_w;i++){
			for(int j=0;j<img_h;j++){
				if(ROI[i][j]){
					RGB c = rgbimg[i][j];
					int index_r = floor(double(c.r) / bin_sizes[0]+0.5);
					int index_g = floor(double(c.g) / bin_sizes[1]+0.5);
					int index_b = floor(double(c.b) / bin_sizes[2]+0.5);
					if(weighted)
						bins[index_r][index_g][index_b] = bins[index_r][index_g][index_b] + weights[i][j];
					else
						bins[index_r][index_g][index_b] += 1.0;
				}
			}
		}
	}
	void uniform_splatter(const Table2D<RGB> & rgbimg, Table2D<bool> ROI)
	{
		img_w = rgbimg.getWidth();
		img_h = rgbimg.getHeight();
		
		Table2D<vector<bool> > present = Table2D<vector<bool> >(256,256,vector<bool>(256,false));

		bin_nums[0] = floor(255/bin_sizes[0]+0.5)+1;
		bin_nums[1] = floor(255/bin_sizes[1]+0.5)+1;
		bin_nums[2] = floor(255/bin_sizes[2]+0.5)+1;

		bins = Table2D<vector<double> >(bin_nums[0],bin_nums[1],vector<double>(bin_nums[2],0));
		for(int i=0;i<img_w;i++){
			for(int j=0;j<img_h;j++){
				if(ROI[i][j]){
					RGB c = rgbimg[i][j];
					if(present[c.r][c.g][c.b]) continue;
					int index_r = floor(double(c.r) / bin_sizes[0]+0.5);
					int index_g = floor(double(c.g) / bin_sizes[1]+0.5);
					int index_b = floor(double(c.b) / bin_sizes[2]+0.5);
					bins[index_r][index_g][index_b] +=1.0;
					present[c.r][c.g][c.b] = true;
				}
			}
		}
	}
	int count_unique(const Table2D<RGB> & rgbimg, Table2D<bool> ROI)
	{
		img_w = rgbimg.getWidth();
		img_h = rgbimg.getHeight();
		
		Table2D<vector<bool> > present = Table2D<vector<bool> >(256,256,vector<bool>(256,false));
		int unique = 0;
		for(int i=0;i<img_w;i++){
			for(int j=0;j<img_h;j++){
				if(ROI[i][j]){
					RGB c = rgbimg[i][j];
					if(false==present[c.r][c.g][c.b]){
						present[c.r][c.g][c.b] = true;
						unique ++;
					}
				}
			}
		}
		return unique;
	}
	void convolve()
	{
		int half_window = sample_rate * 3;
		double * kernel_w = new double[2*half_window+1];
		for(int shift=-half_window;shift<=half_window;shift++){
				kernel_w[shift+half_window] =exp(-0.5*(pow((double)shift,2))/pow(sample_rate,2));
		}
		
		Table2D<vector<double> > conv_w1 = Table2D<vector<double> >(bin_nums[0],bin_nums[1],vector<double>(bin_nums[2],0));
		Table2D<vector<double> > conv_w2 = Table2D<vector<double> >(bin_nums[0],bin_nums[1],vector<double>(bin_nums[2],0));
		conv_w = Table2D<vector<double> >(bin_nums[0],bin_nums[1],vector<double>(bin_nums[2],0));

		for(int i=0;i<bin_nums[0];i++){
			for(int j=0;j<bin_nums[1];j++){
				for(int k=0;k<bin_nums[2];k++){
					for(int shift=-half_window;shift<=half_window;shift++){
						if(((i+shift)>=0)&&((i+shift)<bin_nums[0]))
						{
							double weight = kernel_w[shift+half_window];
							conv_w1[i][j][k] += bins[i+shift][j][k]*weight;
						}
					}
				}
			}
		}
		
		for(int i=0;i<bin_nums[0];i++){
			for(int j=0;j<bin_nums[1];j++){
				for(int k=0;k<bin_nums[2];k++){
					for(int shift=-half_window;shift<=half_window;shift++){
						if(((j+shift)>=0)&&((j+shift)<bin_nums[1]))
						{
							double weight = kernel_w[shift+half_window];
							conv_w2[i][j][k] += conv_w1[i][j+shift][k]*weight;
						}
					}
				}
			}
		}

		for(int i=0;i<bin_nums[0];i++){
			for(int j=0;j<bin_nums[1];j++){
				for(int k=0;k<bin_nums[2];k++){
					for(int shift=-half_window;shift<=half_window;shift++){
						if(((k+shift)>=0)&&((k+shift)<bin_nums[2]))
						{
							double weight = kernel_w[shift+half_window];
							conv_w[i][j][k] += conv_w2[i][j][k+shift]*weight;
						}
					}
				}
			}
		}

		delete kernel_w;
		kernel_w = NULL;
	}
	inline double slice(RGB rgb_p)
	{
		int r_0 = min(floor((double)(rgb_p.r) / bin_sizes[0]),(double)bin_nums[0]-2);
		int g_0 = min(floor((double)(rgb_p.g) / bin_sizes[1]),(double)bin_nums[1]-2);
		int b_0 = min(floor((double)(rgb_p.b) / bin_sizes[2]),(double)bin_nums[2]-2);
		int r_1 = r_0+1;
		int g_1 = g_0+1;
		int b_1 = b_0+1;
		double ratio_r = ((double)(rgb_p.r) - r_0*bin_sizes[0])/bin_sizes[0];
		double ratio_g = ((double)(rgb_p.g) - g_0*bin_sizes[1])/bin_sizes[1];
		double ratio_b = ((double)(rgb_p.b) - b_0*bin_sizes[2])/bin_sizes[2];
		double c00 = conv_w[r_0][g_0][b_0]*(1-ratio_r)+conv_w[r_1][g_0][b_0]*ratio_r;
		double c10 = conv_w[r_0][g_1][b_0]*(1-ratio_r)+conv_w[r_1][g_1][b_0]*ratio_r;
		double c01 = conv_w[r_0][g_0][b_1]*(1-ratio_r)+conv_w[r_1][g_0][b_1]*ratio_r;
		double c11 = conv_w[r_0][g_1][b_1]*(1-ratio_r)+conv_w[r_1][g_1][b_1]*ratio_r;

		double c0 = c00*(1-ratio_g) + c10*ratio_g;
		double c1 = c01*(1-ratio_g) + c11*ratio_g;
		double c = c0*(1-ratio_b) + c1*ratio_b;
		return c;
	}
	inline Table2D<double> slice(const Table2D<RGB> & rgbimg,Table2D<double> weights = Table2D<double>(1,1,1))
	{
		Table2D<double> w_data(img_w,img_h,0);
		bool weighted = true;
		if(1==weights.getWidth()) weighted = false;
		for(int i=0;i<img_w;i++){
			for(int j=0;j<img_h;j++){
				RGB rgb_p = rgbimg[i][j];
				w_data[i][j] = this->slice(rgb_p);
				if(weighted) w_data[i][j] = w_data[i][j]*weights[i][j];
			}
		}
		return w_data;
	}
	inline double slicesum(const Table2D<RGB> & rgbimg,Table2D<bool> ROI)
	{
		double w_data_sum = 0;
		for(int i=0;i<img_w;i++){
			for(int j=0;j<img_h;j++){
				if(ROI[i][j]){
					RGB rgb_p = rgbimg[i][j];
					w_data_sum += this->slice(rgb_p);
				}
			}
		}
		return w_data_sum;
	}
	// for validation
	static Table2D<double> bruteforce(const Table2D<RGB> & rgbimg,Table2D<bool> ROI, double sigma_)
	{
		int w = rgbimg.getWidth();
		int h = rgbimg.getHeight();
		Table2D<double> w_data(w,h,0);
		for(int i=0;i<w;i++){
			for(int j=0;j<h;j++){
				double w_data_sum = 0;
				RGB c1 = rgbimg[i][j];
				for(int x=0;x<w;x++){
					for(int y=0;y<h;y++){
						if(ROI[x][y]){
							RGB c2 = rgbimg[x][y];
							double d_r = (double)(c1.r) - (double)(c2.r);
							double d_g = (double)(c1.g) - (double)(c2.g);
							double d_b = (double)(c1.b) - (double)(c2.b);
							w_data_sum += exp(-(d_r*d_r+d_g*d_g+d_b*d_b)/2.0/sigma_/sigma_);
						}
					}
				}
				w_data[i][j]=w_data_sum;
			}
		}
		return w_data;
	}

	static Table2D<double> bruteforce_adaptive(const Table2D<RGB> & rgbimg,Table2D<bool> ROI, double sigmafactor, Table2D<double> sigmas, double * imgstds)
	{
		int w = rgbimg.getWidth();
		int h = rgbimg.getHeight();
		Table2D<double> w_data(w,h,0);
		for(int i=0;i<w;i++){
			for(int j=0;j<h;j++){
				double w_data_sum = 0;
				RGB c1 = rgbimg[i][j];
				for(int x=0;x<w;x++){
					for(int y=0;y<h;y++){
						if(ROI[x][y]){
							RGB c2 = rgbimg[x][y];
							double d_r = ((double)(c1.r) - (double)(c2.r)) / imgstds[0];
							double d_g = ((double)(c1.g) - (double)(c2.g)) / imgstds[1];
							double d_b = ((double)(c1.b) - (double)(c2.b)) / imgstds[2];
							double sigma_p = sigmafactor * sigmas[i][j];
							double sigma_q = sigmafactor * sigmas[x][y];
							w_data_sum += exp(-(d_r*d_r+d_g*d_g+d_b*d_b)/sigma_p/sigma_p)+exp(-(d_r*d_r+d_g*d_g+d_b*d_b)/sigma_q/sigma_q);
						}
					}
				}
				w_data[i][j]=w_data_sum;
			}
		}
		return w_data;
	}

	static Table2D<double> bruteforce_adaptive2(const Table2D<RGB> & rgbimg,Table2D<bool> ROI, double sigmafactor, Table2D<double> sigmas, double * imgstds)
	{
		int w = rgbimg.getWidth();
		int h = rgbimg.getHeight();
		static Table2D<vector<int> > neighbors(w,h,vector<int>());
		static bool neighborcomputed = false;
		if(!neighborcomputed){
			printf("computing neighbors\n");
			for(int i=0;i<w;i++){
				for(int j=0;j<h;j++){
					RGB c1 = rgbimg[i][j];
					double sigma_p = sigmafactor * sigmas[i][j];
					for(int x=0;x<w;x++){
						for(int y=0;y<h;y++){
							RGB c2 = rgbimg[x][y];
							double d_r = ((double)(c1.r) - (double)(c2.r)) / imgstds[0];
							double d_g = ((double)(c1.g) - (double)(c2.g)) / imgstds[1];
							double d_b = ((double)(c1.b) - (double)(c2.b)) / imgstds[2];
							double dist = sqrt(d_r*d_r+d_g*d_g+d_b*d_b);
							if(dist < 2 * sigma_p){
								neighbors[i][j].push_back(x+j*w);
							}
						}
					}
					outv(neighbors[i][j].size());
				}
			}
			neighborcomputed = true;
		}
		Table2D<double> w_data(w,h,0);
		exit(-1);
		return w_data;
	}

	// accelerated bruteforce
	static Table2D<double> bruteforce_A(const Table2D<RGB> & rgbimg,Table2D<bool> ROI, double sigma_,double sigma2_=-1,double sigma3_=-1,Table2D<double> weights = Table2D<double>(1,1,1))
	{
		bool weighted = true;
		if(1==weights.getWidth()) weighted = false;
		//cout<<"Accelerated bruteforce\n";
		int w = rgbimg.getWidth();
		int h = rgbimg.getHeight();
		static bool firsttime =true;
		static Table2D<vector<vector<int> > > colorcounts= Table2D<vector<vector<int> > >(256,256,vector<vector<int> >(256));;
		static vector<RGB> compactcolors;
		static vector<double> compactdensities;
		static Table2D<double> oneDprob;
		static vector<int> hists;
		static Table2D<int>  compactlabel(w,h); // starts from zero
		static vector<RGB> sigmas;
		if(firsttime){
			int temp_gap = 1.0;
			while(1){
				Table2D<vector<bool> > colorflag = Table2D<vector<bool> > (256,256,vector<bool>(256,false));
				for(int i=0;i<w;i++){
					for(int j=0;j<h;j++){
						RGB c = rgbimg[i][j];
						c.r = c.r - c.r%temp_gap;
						c.g = c.g - c.g%temp_gap;
						c.b = c.b - c.b%temp_gap;
						colorflag[c.r][c.g][c.b] = true;
					}
				}
				int color_count = 0;
				for(int r =0;r<256;r++){
					for(int g =0;g<256;g++){
						for(int b =0;b<256;b++){
							if(colorflag[r][g][b]) color_count++;
						}
					}
				}
				if(color_count<20000) break;
				else temp_gap ++;
			}
			for(int i=0;i<w;i++){
				for(int j=0;j<h;j++){
					int pixel_idx = i+j*w;
					RGB c = rgbimg[i][j];
					c.r = c.r - c.r%temp_gap;
					c.g = c.g - c.g%temp_gap;
					c.b = c.b - c.b%temp_gap;
					colorcounts[c.r][c.g][c.b].push_back(pixel_idx);
				}
			}
			oneDprob = getoneDprob(rgbimg,temp_gap);
			// normalzie oneDprob so that max is 1.0
			double maxprob = oneDprob.getMax();
			for(int i=0;i<256;i++){
				for(int c=0;c<3;c++){
					oneDprob[i][c] = oneDprob[i][c] / maxprob;
				}
			}
			outv(oneDprob.getMax());
			outv(oneDprob.getMin());
			cout<<"densiteis in color space"<<endl;
			Table2D<vector<double> > densities = getcolordensities(rgbimg);
			cout<<"linear scan of the color space\n";
			for(int r =0;r<256;r++){
				for(int g =0;g<256;g++){
					for(int b =0;b<256;b++){
						if(!colorcounts[r][g][b].empty()){
							hists.push_back(colorcounts[r][g][b].size());
							compactcolors.push_back(RGB(r,g,b));
							//compactdensities.push_back((double)(colorcounts[r][g][b].size())/w/h);
							compactdensities.push_back(densities[r][g][b]);
							int compact_l = hists.size()-1;
							for(int i=0;i<hists[compact_l];i++)
							{
								int pixel_idx = colorcounts[r][g][b][i];
								compactlabel[pixel_idx%w][pixel_idx/w] = compact_l;
							}
						
						}
					}
				}
			}
			outv(hists.size());
			outv(temp_gap);
			firsttime = false;
			double maxdensity = -1e20;
			double mindensity = 1e20;
			for(int i=0;i<hists.size();i++){
				if(compactdensities[i]>maxdensity) maxdensity = compactdensities[i];
				if(compactdensities[i]<mindensity) mindensity = compactdensities[i];
			}
			outv(maxdensity);
			outv(mindensity);
			for(int i=0;i<hists.size();i++){ // densities in the range [0 1]
				compactdensities[i] = compactdensities[i] / maxdensity;
			}
			// for computing adaptive sigmas
			for(int i=0;i<hists.size();i++){
				RGB c = compactcolors[i];
				double f_r = oneDprob[c.r][0];
				double f_g = oneDprob[c.g][1];
				double f_b = oneDprob[c.b][2];
				double maxratio = 2.0;
				double minratio = 0.5;
				double sigma_r = min(sigma_  / (f_r*(1.0/minratio - 1.0/maxratio)+1.0/maxratio),255.0);
				double sigma_g = min(sigma2_ / (f_g*(1.0/minratio - 1.0/maxratio)+1.0/maxratio),255.0);
				double sigma_b = min(sigma3_ / (f_b*(1.0/minratio - 1.0/maxratio)+1.0/maxratio),255.0);
				//double sigma_r = sigma_;
				//double sigma_g = sigma2_;
				//double sigma_b = sigma3_;
				sigmas.push_back(RGB(sigma_r,sigma_g,sigma_b));
			}
		}
		
		int num_color = hists.size();
		// linear scan of ROI to collect hist
		vector<double> ROI_hist(num_color,0);
		for(int i=0;i<w;i++){
			for(int j=0;j<h;j++){
				if(ROI[i][j]&&weighted)
					ROI_hist[compactlabel[i][j]] += weights[i][j];
				else if(ROI[i][j]&&(!weighted))
					ROI_hist[compactlabel[i][j]] += 1;
			}
		}
		// linear scan of compact colors
		vector<double> linear_sum(num_color,0);
		for(int i=0;i<num_color;i++){
			RGB c1 = compactcolors[i];
			for(int j=0;j<num_color;j++){
				if(ROI_hist[j]){
					RGB c2 = compactcolors[j];
					double d_r = (double)(c1.r) - (double)(c2.r);
					double d_g = (double)(c1.g) - (double)(c2.g);
					double d_b = (double)(c1.b) - (double)(c2.b);
					if(sigma2_<0)
						linear_sum[i] += exp(-(d_r*d_r+d_g*d_g+d_b*d_b)/2.0/sigma_/sigma_)*ROI_hist[j];
					else{
						//double sigmap = 6.0/(sqrt(compactdensities[i])+0.08);
						//double sigmaq = 6.0/(sqrt(compactdensities[j])+0.08);
						//linear_sum[i] += exp(-(d_r*d_r+d_g*d_g+d_b*d_b)/2.0/sigmap/sigmaq)*ROI_hist[j];
						linear_sum[i] += exp(-d_r*d_r/2.0/sigma_/sigma_-d_g*d_g/2.0/sigma2_/sigma2_-d_b*d_b/2.0/sigma3_/sigma3_)*ROI_hist[j];
						//if(d_r*d_r+d_g*d_g+d_b*d_b>75)
						//linear_sum[i] += exp(-d_r*d_r/2.0/sigmas[i].r/sigmas[j].r-d_g*d_g/2.0/sigmas[i].g/sigmas[j].g-d_b*d_b/2.0/sigmas[i].b/sigmas[j].b)*ROI_hist[j];
					}
				}
			}
		}
		Table2D<double> w_data(w,h,0);
		for(int i=0;i<w;i++){
			for(int j=0;j<h;j++){
				w_data[i][j] = linear_sum[compactlabel[i][j]]*(weighted?weights[i][j]:1.0);
			}
		}
		//exit(-1);
		return w_data;
	}
};

double unique_sum(Table2D<double> & table, const Table2D<RGB> & rgbimg, Table2D<bool> ROI)
{
	int img_w = rgbimg.getWidth();
	int img_h = rgbimg.getHeight();
		
	Table2D<vector<bool> > present = Table2D<vector<bool> >(256,256,vector<bool>(256,false));
	double returnv = 0;
	for(int i=0;i<img_w;i++){
		for(int j=0;j<img_h;j++){
			if(ROI[i][j]){
				RGB c = rgbimg[i][j];
				if(false==present[c.r][c.g][c.b]){
					present[c.r][c.g][c.b] = true;
					returnv += table[i][j];
				}
			}
		}
	}
	return returnv;
}

Table2D<double> getcolorweights(const Table2D<RGB> & rgbimg, bool smoothed)
{
	int img_w = rgbimg.getWidth();
	int img_h = rgbimg.getHeight();
	Table2D<double> weights(img_w,img_h,0);
	Table2D<vector<double> > colorcount = Table2D<vector<double> >(256,256,vector<double>(256,0));

	for(int i=0;i<img_w;i++){
		for(int j=0;j<img_h;j++){
			RGB c = rgbimg[i][j];
			colorcount[c.r][c.g][c.b] += 1 ;
			if(smoothed){
				if(c.r>=1) colorcount[c.r-1][c.g][c.b] += 0.5 ;
				if(c.g>=1) colorcount[c.r][c.g-1][c.b] += 0.5 ;
				if(c.b>=1) colorcount[c.r][c.g][c.b-1] += 0.5 ;
				if(c.r<=254) colorcount[c.r+1][c.g][c.b] += 0.5 ;
				if(c.g<=254) colorcount[c.r][c.g+1][c.b] += 0.5 ;
				if(c.b<=254) colorcount[c.r][c.g][c.b+1] += 0.5 ;
			}
		}
	}
	for(int i=0;i<img_w;i++){
		for(int j=0;j<img_h;j++){
			RGB c = rgbimg[i][j];
			weights[i][j] = 1.0 / colorcount[c.r][c.g][c.b];
		}
	}
	outv(weights.getMax());
	outv(weights.getMin());
	return weights;
}

Table2D<vector<double> > getcolordensities(const Table2D<RGB> & rgbimg){
	Table2D<vector<double> > densities = Table2D<vector<double> > (256,256,vector<double>(256,0));
	Filtering density_fil(8,8,8);
	density_fil.splatter(rgbimg,Table2D<bool>(rgbimg.getWidth(),rgbimg.getHeight(),true));
	density_fil.convolve();
	for(int r =0;r<256;r++){
		for(int g =0;g<256;g++){
			for(int b =0;b<256;b++){
				densities[r][g][b] = density_fil.slice(RGB(r,g,b));
			}
		}
	}
	return densities;
}

#endif
