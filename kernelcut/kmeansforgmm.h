#ifndef _KMEANSFORGMM_H_
#define _KMEANSFORGMM_H_
#include <iostream>
#include <vector>
using namespace std;

class KMeans{
public: // for GMM only
	vector<double> wei_com; // weight of Gaussian components
	vector<vector<double> > mu_com; // mean of Gaussian components
	vector<double *> cov_com;
private:
	vector<int> labeling_hist;
	vector<vector<double> > data;
	vector<vector<double> > centers;
	vector<int> labeling;
	int num_dims;
	int num_points;
	int K;
public:
	KMeans(){}
	KMeans(vector<vector<double> > data, int K)
	{
		this->data = data;
		this->num_points = data.size();
		this->num_dims = data[0].size();
		centers = vector<vector<double> >(K,vector<double>(num_dims,0));
		this->K = K;
		this->labeling = vector<int>(this->num_points,0);
	}
	double getSSD(){
		double ssd = 0;
		for(int i=0;i<num_points;i++){
			ssd = ssd + getdist(data[i],centers[labeling[i]]);
		}
		return ssd;
	}
	void optimizelabeling(){
		for(int i=0;i<num_points;i++){
			double min_dist = 1e+20;
			for(int k=0;k<K;k++){
				double dist_i_k = getdist(data[i],centers[k]);
				if(dist_i_k<min_dist){
					labeling[i] = k;
					min_dist = dist_i_k;
				}
			}
		}
	}
	void optimizecenters(){
		vector<vector<double> > centers_sum(K,vector<double>(num_dims,0));
		labeling_hist =vector<int>(K,0);
		for(int i=0;i<num_points;i++){
			vectorsum(data[i],centers_sum[labeling[i]]);
			labeling_hist[labeling[i]] = labeling_hist[labeling[i]]+1;
		}
		for(int k=0;k<K;k++){
			if(labeling_hist[k])
				vectormultiply(centers_sum[k],1.0/labeling_hist[k]);
		}
		centers = centers_sum;
	}
	// optimize from one set of initializations
	double optimizeoneinit(const vector<vector<double> > & initial_centers){
		centers = initial_centers;
		resetlabeling();
		double init_ssd = getSSD();
		//cout<<"init ssd: "<<init_ssd<<endl;
		//print();
		while(1){
			optimizelabeling();
			optimizecenters();
			double new_ssd = getSSD();
			//cout<<"new ssd: "<<new_ssd<<endl;
			//print();
			if(init_ssd - new_ssd > 10)
				init_ssd = new_ssd;
			else
				break;
		}
		return init_ssd;
	}
	void optimizemultiinit(int num_init){
		vector<vector<double> > best_initial_centers;
		double init_ssd = 1e+20;
		for(int i=0;i<num_init;i++){
			vector<vector<double> > initial_centers(K,vector<double>(num_dims,0));
			for(int k=0;k<K;k++)
				initial_centers[k] = data[rand()%num_points];
			double new_ssd = optimizeoneinit(initial_centers);
			//printf("kmeans ssd: %f\n",new_ssd);
			if(init_ssd - new_ssd > 1){
				init_ssd = new_ssd;
				best_initial_centers = initial_centers;
			}
		}
		// optimize with the best initialization
		//printf("best kmeans ssd: %f\n",optimizeoneinit(best_initial_centers));
	}
	void computegmm(){
		wei_com = vector<double>(K,0);
		for(int k=0;k<K;k++){
			wei_com[k] = labeling_hist[k]/(double)num_points;
		}
		mu_com = centers;
		// covariance matrix
		cov_com = vector<double *>(K,NULL);
		for(int k=0;k<K;k++){
			cov_com[k] = new double[num_dims*num_dims];
			for(int d_row=0;d_row<num_dims;d_row++){
				for(int d_col=0;d_col<num_dims;d_col++){
					cov_com[k][d_row*num_dims+d_col] = 0;
				}
			}
		}
		for(int i=0;i<num_points;i++){
			for(int d_row=0;d_row<num_dims;d_row++){
				for(int d_col=0;d_col<num_dims;d_col++){
					cov_com[labeling[i]][d_row*num_dims+d_col] += (data[i][d_row]-mu_com[labeling[i]][d_row])
						*(data[i][d_col]-mu_com[labeling[i]][d_col])/(double)(labeling_hist[labeling[i]]);
				}
			}
		}
	}
	inline void resetlabeling(){
		labeling = vector<int>(num_points,0);
	}
	inline double getdist(const vector<double>& point, const vector<double> & center){
		double dist = 0;
		for(int d=0;d<num_dims;d++)
			dist = dist + (point[d]-center[d])*(point[d]-center[d]);
		return dist;
	}
	inline void vectorsum(const vector<double>& source, vector<double> & dest){
		for(int d=0;d<num_dims;d++)
			dest[d] = dest[d] + source[d];
	}
	inline void vectormultiply(vector<double> & dest, double w){
		for(int d=0;d<num_dims;d++)
			dest[d] = dest[d]*w;
	}
	void print(){
		cout<<"Kmeans of "<<K<<" centers:\n";
		for(int k=0;k<K;k++){
			for(int d=0;d<num_dims;d++){
				cout<<centers[k][d]<<'\t';
			}
			cout<<endl;
		}
	}
};
#endif
