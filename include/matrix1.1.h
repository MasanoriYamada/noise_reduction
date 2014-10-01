#ifndef INCLUDED_MATRIX_H
#define INCLUDED_MATRIX_H

#include <complex>
#define I       std::complex<double>(0.0,1.0)
typedef std::complex<double> COMPLEX;

class Matrix{
	
	int row;  //行
	int column;  //列
	
	COMPLEX** p_top; //配列の最初を指すポインタ
	
public:
	Matrix(int i=1, int j=1);//コンストラクタ
	Matrix(const Matrix &cp);//コピーコンストラクタ
	
	~Matrix();//デストラクタ
	
	int row_size(){ return(row); }
	int column_size(){ return(column); }
	COMPLEX rd(int,int);
	
	void change_size(int i, int j);//行列のサイズを変更する
	
	//演算子のオーバーロード
	COMPLEX* &operator[](int i){ return(p_top[i]); }
	Matrix operator=(const Matrix &a);
	Matrix operator+(const Matrix &a);
	Matrix operator-(const Matrix &a);
	Matrix operator*(const Matrix &a);
	bool operator==(const Matrix &a);
	bool operator!=(const Matrix &a);
	
	friend Matrix operator*(const Matrix &a, COMPLEX b);
	friend Matrix operator*(COMPLEX b, const Matrix &a);
	
	//行列の変換など
	void unit_matrix();//単位行列にする
	Matrix transposed();//転置行列をかえす
};

#endif