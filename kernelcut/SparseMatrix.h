#pragma once
#include<vector>

template<class T>
// Tri-tuple
struct Trituple
{
public:
	Trituple():row(0),col(0),val(0){}
	Trituple(int row_, int col_, T val_)
		:row(row_),col(col_),val(val_){}
	~Trituple(){};
	int row;
	int col;
	T val;
};

// sparse matrix
template<class T>
class SparseMatrix
{
public:
	SparseMatrix():width(0),height(0),nnz(0){}
	SparseMatrix(int width_,int height_)
	{
		width = width_;
		height = height_;
		nnz=0;
		data = vector<Trituple<T> >();
	}
	SparseMatrix(int width_,int height_,const vector<Trituple<T> > & data_)
		:width(width_),height(height_)
	{
		data = data_;
		nnz = data.size();
	}
	~SparseMatrix(){}
	void print()
	{
		for(int i=0;i<nnz;i++)
			cout<<data[i].row<<" "<<data[i].col<<" "<<data[i].val<<endl;
	}
	Trituple<T> operator[](int i)
	{
		return data[i];
	}
	void add(Trituple<T> entry);
	int getsize(){return nnz;}
	int width,height;
private:
	vector<Trituple<T> > data;
	int nnz;
};

template<class T>
void SparseMatrix<T>::add(Trituple<T> entry)
{
	data.push_back(entry);
	nnz++;
}
