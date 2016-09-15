#ifndef _MATRIXCAL_H_
#define _MATRIXCAL_H_
#include <math.h>
#include <iostream>
using namespace std;

inline double * matrixMultiply(double a[], double b[],int m,int n,int k);
inline void matrixMultiply(double a[], int M, int N, double weight);
inline void matrixSum(double source[], double dest[], int M, int N);
inline double * matrixTranspose(double a[], int m, int n);
inline double matrixDeterminant(double a[], int n);
inline double * matrixInverse(double a[], int N);
inline void matrixPrint(double a[], int N);
inline double * diamatrix(int N, double v);

// matrix multiply, return the new matrix
inline double * matrixMultiply(double a[], double b[],int m,int n,int k)  
{  
    int i,j,l,u;
	double * c = new double[m*k];
    for(i=0;i<m;i++){  
		for(j=0;j<k;j++){  
            u=i*k+j;
			c[u]=0.0;  
            for(l=0;l<=n-1;l++)  
                c[u] += a[i*n+l]*b[l*k+j];
        }  
	}
	return c;
}

// matrix multiply, return the new matrix
inline void matrixMultiply(double a[], int M, int N, double weight)  
{  
    for(int m=0;m<M;m++){
		for(int n=0;n<N;n++){
			a[n+m*N] = a[n+m*N]*weight;
		}
	}
}

inline void matrixSum(double source[], double dest[], int M, int N)
{
	for(int m=0;m<M;m++){
		for(int n=0;n<N;n++){
			dest[n+m*N] = dest[n+m*N] + source[n+m*N];
		}
	}
}
// matrix transpose, return the new matrix
inline double * matrixTranspose(double a[], int m, int n)
{
	double * b = new double[m*n];
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			b[i+j*m]=a[j+i*n];
		}
	}
	return b;
}

inline double matrixDeterminant(double a[], int N)
{
	if(1==N)
		return a[0];
	else if(2==N)
		return a[0]*a[3]-a[1]*a[2];
	// for n>=3
	double det = 0;
	for(int i=0;i<N;i++){
		double * Mi0=new double[(N-1)*(N-1)];
		int idx=0;
		for(int j=0;j<N;j++){
			for(int k=0;k<N;k++){
				if(k==0)continue;
				if(j==i)continue;
				Mi0[idx++] = a[k+j*N];
			}
		}
		double det_Mi0 = matrixDeterminant(Mi0,N-1);
		delete [] Mi0;
		det = det + pow((double)-1,i+1+1)*det_Mi0*a[i*N+0];
	}
	return det;
}

inline void matrixPrint(double a[], int N)
{
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			cout<<a[j+i*N]<<'\t';
		}
		cout<<endl;
	}
}

inline double * matrixAdj(double a[], int N)
{
	double * b = new double[N*N];
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			double * Mij=new double[(N-1)*(N-1)];
			int idx=0;
			for(int m=0;m<N;m++){
				for(int n=0;n<N;n++){
					if(n==j)continue;
					if(m==i)continue;
					Mij[idx++] = a[n+m*N];
				}
			}
			b[j+i*N] = pow((double)-1,i+j+2)*matrixDeterminant(Mij, N-1);
			delete [] Mij;
		}
	}
	double * b_t = matrixTranspose(b,N,N);
	delete [] b;
	return b_t;
}

inline double * matrixInverse(double a[], int N)
{
	double * a_adj = matrixAdj(a, N);
	double a_det = matrixDeterminant(a,N);
	if(abs(a_det)<1e-3){
		//cout<<"determinant is zero!!!!!!!!!!!!\n";
		return NULL;
	}
	double * a_inv = new double[N*N];
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			a_inv[j+i*N] = a_adj[j+i*N]/a_det;
		}
	}
	delete [] a_adj;
	return a_inv;
}


inline double * diamatrix(int N, double v)
{
	double * a = new double[N*N];
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			if(i==j)
				a[j+i*N] = v;
			else
				a[j+i*N] = 0;
		}
	}
	return a;
}
/* client program
double a[9]={1,0,0,0,2,0,0,0,4};
matrixPrint(a,3);
double * a_t=matrixTranspose(a,3,3);
matrixPrint(a_t,3);
delete [] a_t;
double a_det =matrixDeterminant(a,3);
outv(a_det);
double * a_adj=matrixAdj(a,3);
matrixPrint(a_adj,3);
delete [] a_adj;
double * a_inv = matrixInverse(a, 3);
matrixPrint(a_inv,3);
delete [] a_inv;
return -1;
*/
#endif
