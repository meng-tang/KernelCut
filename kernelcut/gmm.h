#ifndef _GMM_H_
#define _GMM_H_
#include <vector>
#include "matrixcal.h"
#include "kmeansforgmm.h"

#define INFSMALL 1.0e-60
#define GMM_PI 3.1415926
struct RGB_STD
{
	double r_std;
	double g_std;
	double b_std;
};
inline double normdist1D(double x,double mu,double sigma);
inline double normdist3D(RGB c,RGB mu,RGB_STD sigma);
inline double normdist3D(RGB c,RGB mu,double * cov_inv, double cov_det);

inline double normdist1D(double x,double mu,double sigma)
{
	return exp(-(x-mu)*(x-mu)/(2*sigma*sigma))/sigma/sqrt(2*GMM_PI);
}

inline double normdist3D(RGB c,RGB mu,RGB_STD sigma)
{
	return normdist1D(c.r,mu.r,sigma.r_std)*normdist1D(c.g,mu.g,sigma.g_std)*normdist1D(c.b,mu.b,sigma.b_std)
		+INFSMALL;
}

inline double normdist3D(RGB c,RGB mu,double * cov_inv, double cov_det)
{
	double * xminusmu = new double[3]; // 1 row 3 col
	xminusmu[0] = (double)(c.r) - (double)(mu.r);
	xminusmu[1] = (double)(c.g) - (double)(mu.g);
	xminusmu[2] = (double)(c.b) - (double)(mu.b);
	double * xminusmu_t = matrixTranspose(xminusmu,3,1); // 3 row 1 col
	double * mul1 = matrixMultiply(xminusmu,cov_inv,1,3,3); // 1 row 3 col
	double * mul2 = matrixMultiply(mul1,xminusmu_t,1,3,1); // 1 row 1 col (schalar)
	double mul = mul2[0];
	delete [] xminusmu;
	delete [] xminusmu_t;
	delete [] mul1;
	delete [] mul2;
	return 1.0/pow(2*GMM_PI,3/2.0)/sqrt(cov_det)*exp(-0.5*mul)+INFSMALL;
}
class GMM{
private:
	vector<RGB> pixels;
	int num_pix;
	int num_com; // number of Gaussian components
	int num_dim; // number of dimensions
	vector<double> wei_com; // weight of Gaussian components
	vector<RGB> mu_com; // mean of Gaussian components
	vector<RGB_STD> std_com; // variance of Gaussian components
	int num_itr;
	vector<vector<double> > wpk;

	vector<double *> cov_com; // covariance matrix
	vector<double *> cov_inv_com; // inverse of covariance matrix
	vector<double> cov_det_com; // determinant of covariance matrix
public:
	GMM(){}
	/*GMM(const GMM & gmm_){
		pixels = gmm_.pixels;
		num_pix = gmm_.num_pix;
		num_com = gmm_.num_com;
		num_dim = gmm_.num_dim;
		wei_com = gmm_.wei_com; 
		mu_com = gmm_.mu_com;
		std_com = gmm_.std_com; 
		num_itr = gmm_.num_itr;
		wpk = gmm_.wpk;

		cov_com = vector<double *>(num_com);
		cov_inv_com = vector<double *>(num_com);
		for(int k=0;k<num_com;k++){
			cov_com[k] = new double[9];
			cov_inv_com[k] = new double[9];
			for(int i=0;i<9;i++){
				cov_com[k][i] = gmm_.cov_com[k][i];
				cov_inv_com[k][i] = gmm_.cov_inv_com[k][i];
			}
		}
		cov_det_com = gmm_.cov_det_com; 
	}*/
	GMM(const Table2D<RGB> & rgbimg, Table2D<bool> ROI,int _num_com=10){
		num_com = _num_com;
		num_dim = 3;
		num_itr = 10;
		feeddata(rgbimg, ROI);
		initbykmeans();
	}
	void naiveinit()
	{
		wei_com = vector<double>(num_com,1.0/num_com);
		RGB_STD a={10, 10, 10};
		std_com = vector<RGB_STD>(num_com,a);
		cov_com = vector<double *>(num_com);
		cov_inv_com = vector<double *>(num_com);
		cov_det_com = vector<double>(num_com);
		mu_com = vector<RGB>(num_com);
		for(int k=0;k<num_com;k++){
			mu_com[k] = pixels[rand()%num_pix];
			cov_com[k] = diamatrix(num_dim, 100);
			cov_inv_com[k] = matrixInverse(cov_com[k],num_dim);
			cov_det_com[k] = matrixDeterminant(cov_com[k],num_dim);
		}
	}
	void feeddata(const Table2D<RGB> & rgbimg, Table2D<bool> ROI){
		pixels = vector<RGB>(countintable(ROI,true));
		int img_w = rgbimg.getWidth();
		int img_h = rgbimg.getHeight();
		int pix_id =0;
		for(int x=0;x<img_w;x++){
			for(int y=0;y<img_h;y++){
				if(ROI[x][y])
					pixels[pix_id++] = rgbimg[x][y];
			}
		}
		num_pix = pixels.size();
		wpk = vector<vector<double> >(num_pix,vector<double>(num_com,0));
	}
	void initbykmeans(){
		vector<vector<double> > data(num_pix,vector<double>(3,0));
		for(int i=0;i<num_pix;i++){
			data[i][0] = pixels[i].r;
			data[i][1] = pixels[i].g;
			data[i][2] = pixels[i].b;
		}
		KMeans kmeans(data,num_com);
		kmeans.optimizemultiinit(5);
		//kmeans.print();
		kmeans.computegmm();
		this->wei_com = kmeans.wei_com;
		vector<vector<double> > init_mu = kmeans.mu_com;
		cov_com = vector<double *>(num_com);
		cov_inv_com = vector<double *>(num_com);
		cov_det_com = vector<double>(num_com);
		mu_com = vector<RGB>(num_com);
		for(int k=0;k<this->num_com;k++){
			this->mu_com[k].r = floor(init_mu[k][0]);
			this->mu_com[k].g = floor(init_mu[k][1]);
			this->mu_com[k].b = floor(init_mu[k][2]);
			cov_com[k] = kmeans.cov_com[k];
			cov_com[k][0] = max(cov_com[k][0],1.0);
			cov_com[k][4] = max(cov_com[k][4],1.0);
			cov_com[k][8] = max(cov_com[k][8],1.0);
			cov_inv_com[k] = matrixInverse(cov_com[k],num_dim);
			cov_det_com[k] = matrixDeterminant(cov_com[k],num_dim);
		}
	}
	~GMM()
	{
		for(int k=0;k<num_com;k++){
			delete [] cov_com[k];
			delete [] cov_inv_com[k];
		}
	}
	void optimize(int _num_itr=0)
	{
		if(_num_itr)
			num_itr = _num_itr;
		if(this->num_com == 1)
			num_itr = 1;
		double best_logprob = 1e+100;
		for(int i =0;i<num_itr;i++){
			eStep();
			mStep();
			//double logprob = sumlogprobability();
			//outv(logprob);
			//if(best_logprob-logprob<10)
				//break;
		}
	}
	inline void eStep(){
		for(int i=0;i<num_pix;i++)
		{
			RGB c = pixels[i];
			double probility_sum = getProbability(c);
			for(int k=0;k<num_com;k++)
			{
				wpk[i][k] = wei_com[k]*normdist3D(c,mu_com[k],cov_inv_com[k],cov_det_com[k])/probility_sum;
				//wpk[i][k] = wei_com[k]*normdist3D(c,mu_com[k],std_com[k])/probility_sum;
			}
		}
	}
	inline void mStep(){
		for(int k = 0;k<num_com;k++)
		{
			double sumwp = 0;
			double sumcpwpr = 0;
			double sumcpwpg = 0;
			double sumcpwpb = 0;
			for(int j=0;j<num_pix;j++)
			{
				sumwp+=wpk[j][k];
				sumcpwpr+=wpk[j][k]*(pixels[j].r);
				sumcpwpg+=wpk[j][k]*(pixels[j].g);
				sumcpwpb+=wpk[j][k]*(pixels[j].b);
			}
			wei_com[k]=sumwp/num_pix;
			mu_com[k].r=floor(sumcpwpr/(sumwp+INFSMALL));
			mu_com[k].g=floor(sumcpwpg/(sumwp+INFSMALL));
			mu_com[k].b=floor(sumcpwpb/(sumwp+INFSMALL));
		}
		for(int k = 0;k<num_com;k++)
		{
			RGB mu = mu_com[k];
			double * new_cov = new double[3*3];
			for(int i=0;i<3*3;i++)
				new_cov[i] = 0;
			for(int j=0;j<num_pix;j++)
			{
				RGB c = pixels[j];
				double * xminusmu = new double[3]; // 1 row 3 col
				xminusmu[0] = (double)(c.r) - (double)(mu.r);
				xminusmu[1] = (double)(c.g) - (double)(mu.g);
				xminusmu[2] = (double)(c.b) - (double)(mu.b);
				double * xminusmu_t = matrixTranspose(xminusmu,3,1); // 3 row 1 col
				double * mul1 = matrixMultiply(xminusmu,xminusmu_t,3,1,3); // 3 row 3 col
				delete [] xminusmu;
				delete [] xminusmu_t;
				matrixMultiply(mul1, 3, 3, wpk[j][k]);
				matrixSum(mul1, new_cov, 3, 3);
				delete [] mul1;
			}
			matrixMultiply(new_cov,3,3,1.0/(num_pix*wei_com[k]));
			delete  cov_com[k];
			//new_cov[0] += 0.01;new_cov[4] += 0.01;new_cov[8] += 0.01;
			cov_com[k] = new_cov;
			cov_inv_com[k] = matrixInverse(cov_com[k],num_dim);
			if(cov_inv_com[k] == NULL){ // determinant is zero
				resetsingular(k);
			}
			else
				cov_det_com[k] = matrixDeterminant(cov_com[k],num_dim);
		}
	}
	void resetsingular(int k){
		cout<<"reset singular component "<<k<<endl;
		mu_com[k].r = rand()%255;mu_com[k].g = rand()%255;mu_com[k].b = rand()%255;
		cov_com[k] = diamatrix(3, 100);
		cov_inv_com[k] = matrixInverse(cov_com[k],num_dim);
		cov_det_com[k] = matrixDeterminant(cov_com[k],num_dim);
	}
	void print()
	{
		cout<<"GMM for "<<num_pix<<" pixels:"<<endl;
		for(int k=0;k<num_com;k++)
		{
			cout<<"weight mu "<<wei_com[k]<<" [";
			cout<<(int)mu_com[k].r<<' '<<(int)mu_com[k].g<<' '<<(int)mu_com[k].b<<"]\n";
			//cout<<std_com[k].r_std<<' '<<std_com[k].g_std<<' '<<std_com[k].b_std<<endl;
			cout<<"covariance matrix: (with determinant "<<cov_det_com[k]<<" )\n";
			matrixPrint(cov_com[k],num_dim);
		}
	}
	//mixture of guassian probility in R3
	inline double getProbability(RGB c)
	{
		double probility = 0;
		for(int k=0;k<num_com;k++)
		{
			//probility+=wei_com[k]*normdist3D(c,mu_com[k],std_com[k]);
			probility+=wei_com[k]*normdist3D(c,mu_com[k],cov_inv_com[k],cov_det_com[k]);
		}
		return probility;
	}

	double sumlogprobability()
	{
		double sumlog = 0;
		for(int i=0;i<num_pix;i++){
			sumlog = sumlog - log(getProbability(pixels[i]));
		}
		return sumlog;
	}
	Table2D<double> getProbability(const Table2D<RGB> & rgbimg, Table2D<bool> ROI)
	{
		int img_w = rgbimg.getWidth();
		int img_h = rgbimg.getHeight();
		Table2D<double> probility_t(img_w,img_h,0);
		for(int i=0;i<img_w;i++){
			for(int j=0;j<img_h;j++){
				if(ROI[i][j]){
					probility_t[i][j] = -log(getProbability(rgbimg[i][j]));
				}
			}
		}
		return probility_t;
	}
};
#endif
