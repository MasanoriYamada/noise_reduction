/*
 *  noise_redection0.1.cpp
 *
 *
 *  Created by yamada on 2/4/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "../include/noise_redection.h"
#include "../include/matrix.h"
#include "../include/analys.h"



using namespace std;
Matrix ** MA=new Matrix*[24];
//C1 1
Matrix E(3,3);
//C2 9
Matrix C2_1(3,3);
Matrix C2_2(3,3);
Matrix C2_3(3,3);
Matrix C2_4(3,3);
Matrix C2_5(3,3);
Matrix C2_6(3,3);
Matrix C2_7(3,3);
Matrix C2_8(3,3);
Matrix C2_9(3,3);

//C3 8
Matrix C3_1(3,3);
Matrix C3_2(3,3);
Matrix C3_3(3,3);
Matrix C3_4(3,3);
Matrix C3_5(3,3);
Matrix C3_6(3,3);
Matrix C3_7(3,3);
Matrix C3_8(3,3);

//C4 6
Matrix C4_1(3,3);
Matrix C4_2(3,3);
Matrix C4_3(3,3);
Matrix C4_4(3,3);
Matrix C4_5(3,3);
Matrix C4_6(3,3);



int main(int argc , char** argv){
    //make out put dir
	dir_path=argv[1];
	cout <<"Directory path  ::"<<dir_path<<endl;
	in_dir_path = dir_path;
	out_dir_path = dir_path;
	root_mkdir(dir_path.c_str());
    out_dir_path=out_dir_path + "/impBSwave";
	root_mkdir(out_dir_path.c_str());
    out_dir_path = out_dir_path + "/noise_reduction";
	root_mkdir(out_dir_path.c_str());

    SinkSpin=argv[2];
    SrcSpin=argv[3];



	setC4();

	make_rep();
	MA[0]=&(E);
	MA[1]=&(C2_1);
	MA[2]=&(C2_2);
	MA[3]=&(C2_3);
	MA[4]=&(C2_4);
	MA[5]=&(C2_5);
	MA[6]=&(C2_6);
	MA[7]=&(C2_7);
	MA[8]=&(C2_8);
	MA[9]=&(C2_9);
	MA[10]=&(C3_1);
	MA[11]=&(C3_2);
	MA[12]=&(C3_3);
	MA[13]=&(C3_4);
	MA[14]=&(C3_5);
	MA[15]=&(C3_6);
	MA[16]=&(C3_7);
	MA[17]=&(C3_8);
	MA[18]=&(C4_1);
	MA[19]=&(C4_2);
	MA[20]=&(C4_3);
	MA[21]=&(C4_4);
	MA[22]=&(C4_5);
	MA[23]=&(C4_6);

    for (int tmp=0; tmp<Tshiftsize; tmp++) {
	for (int it=T_in; it<T_fi+1; it++) {
	for (int j=0; j<Confsize; j++) {
	//for (int id=0; id<24; id++) {
	//	cout << "Matrix"<<id<<endl;
	//	display(MA[id]);
	//	check(MA[id]);
	//}
	COMPLEX *Wave_in=new COMPLEX[XYZnodeSites];
	COMPLEX *Wave_out=new COMPLEX[XYZnodeSites];
	COMPLEX *Wave_ave=new COMPLEX[XYZnodeSites];
	call_data(it,tmp,j,Wave_in);
	cout << "map	start"<<endl;
	for (int id=0; id<24; id++) {
		map(MA[id],Wave_in,Wave_out);
	for (int ix=0; ix<XnodeSites; ix++) {
	for (int iy=0; iy<YnodeSites; iy++) {
	for (int iz=0; iz<ZnodeSites; iz++) {
		wave_ave(ix,iy,iz)=wave_out(ix,iy,iz)/24.0+wave_ave(ix,iy,iz);      //you sholud modify here (charactor part)
	}}}}

	cout << "map	end"<<endl;
		out_data(it,tmp,j,Wave_ave);
	delete []Wave_ave;
	delete []Wave_out;
	delete []Wave_in;
	}}}
	delete [] MA;
	cout << "finish"<<endl;
	return 0;
}



//----------------------------------------------
//set to matrix C4
//----------------------------------------------
int setC4(){
	E[1][1]=1;	E[1][2]=0;	E[1][3]=0;
	E[2][1]=0;	E[2][2]=1;	E[2][3]=0;
	E[3][1]=0;	E[3][2]=0;	E[3][3]=1;

	//4x
	C4_1[1][1]=1;	C4_1[1][2]=0;	C4_1[1][3]=0;
	C4_1[2][1]=0;	C4_1[2][2]=0;	C4_1[2][3]=-1;
	C4_1[3][1]=0;	C4_1[3][2]=1;	C4_1[3][3]=0;
	//4y
	C4_2[1][1]=0;	C4_2[1][2]=0;	C4_2[1][3]=1;
	C4_2[2][1]=0;	C4_2[2][2]=1;	C4_2[2][3]=0;
	C4_2[3][1]=-1;	C4_2[3][2]=0;	C4_2[3][3]=0;
	//4z
	C4_3[1][1]=0;	C4_3[1][2]=-1;	C4_3[1][3]=0;
	C4_3[2][1]=1;	C4_3[2][2]=0;	C4_3[2][3]=0;
	C4_3[3][1]=0;	C4_3[3][2]=0;	C4_3[3][3]=1;
	//4x^-
	C4_4[1][1]=1;	C4_4[1][2]=0;	C4_4[1][3]=0;
	C4_4[2][1]=0;	C4_4[2][2]=0;	C4_4[2][3]=1;
	C4_4[3][1]=0;	C4_4[3][2]=-1;	C4_4[3][3]=0;
	//4y^-
	C4_5[1][1]=0;	C4_5[1][2]=0;	C4_5[1][3]=-1;
	C4_5[2][1]=0;	C4_5[2][2]=1;	C4_5[2][3]=0;
	C4_5[3][1]=1;	C4_5[3][2]=0;	C4_5[3][3]=0;
	//4z^-
	C4_6[1][1]=0;	C4_6[1][2]=1;	C4_6[1][3]=0;
	C4_6[2][1]=-1;	C4_6[2][2]=0;	C4_6[2][3]=0;
	C4_6[3][1]=0;	C4_6[3][2]=0;	C4_6[3][3]=1;


	//debug

	assert(C4_1*C4_4 == E);
	assert(C4_2*C4_5 == E);
	assert(C4_3*C4_6 == E);

	return 0;
}




//----------------------------------------------
//construct all representation
//----------------------------------------------
void make_rep(){
cout << "make	rep"<<endl;
	C2_1=C4_1*C4_1;
	C2_2=C4_2*C4_2;
	C2_3=C4_3*C4_3;
	C2_4=C4_1*C4_2*C4_3;
	C2_5=C4_2*C4_3*C4_1;
	C2_6=C4_3*C4_1*C4_2;
	C2_7=C4_4*C4_5*C4_3;
	C2_8=C4_5*C4_6*C4_1;
	C2_9=C4_6*C4_4*C4_2;


	C3_1=C4_1*C4_2;
	C3_2=C4_1*C4_2*C4_1*C4_2;
	C3_3=C4_2*C4_1;
	C3_4=C4_3*C4_2;
	C3_5=C4_1*C4_3;
	C3_6=C4_2*C4_1*C4_2*C4_1;
	C3_7=C4_3*C4_2*C4_3*C4_2;
	C3_8=C4_1*C4_3*C4_1*C4_3;
	cout << "make	rep		end"<<endl;

}

//----------------------------------------------
//input file name
//----------------------------------------------
void call_data(int it,int tmp,int j,COMPLEX data[XYZnodeSites]){
            char fname[200]={0};
	    int mj = j+251;
	    if(mj< 430)
	      {
	      }
	    else if(429 < mj)
	      {
		mj = mj + 1;
	      }
	    sprintf(fname,"%s/Projwave/xyz/spin%sz+0.%sz+0/%s-b-%05d0/NBSwave.+%03d+%03d.000.000.000.%s-b-%05d0.DD_NR.DD_NR",in_dir_path.c_str(),SinkSpin.c_str(),SrcSpin.c_str(),base,mj,it,tshift[tmp],base,mj);

		cout << fname<<"reading now"<<endl;
            read_data(&(fname[0]),data);
    cout<<"it = "<<it<<" Tshift count = "<<tmp<<"  conf =  "<<j<<"reading"<<endl;
}
//----------------------------------------------
//read data
//----------------------------------------------
int read_data(char fname[300],COMPLEX data[XYZnodeSites]){
	fstream infile;

	infile.open(fname,ios::in|ios::binary);
	if (!infile.is_open()) {
		cout << "ERROR file can't open (no exist)"<<endl;
		return EXIT_FAILURE;
	}
	cout << "SUCCESS open reading file"<<endl;
	int id=0;
	while(!infile.eof()){
		infile.read( ( char * ) &data[id], sizeof( COMPLEX ) );
		id=id+1;
	}
	static int tmp=0;
	if (tmp==0) {
		cout <<"reading data size is    ;;"<<id-1<<endl;
		tmp=tmp+1;
	}
	endian_convert((double*)data,XYZnodeSites*2);
	for (int point=0; point<id; point++) {
		//data[point]=rand()%24;
		//cout << data[point]<<endl;
	}
	cout <<"check "<<data[0]<<endl;
	return 0;
}

//----------------------------------------------
//output file name
//----------------------------------------------
void out_data(int it,int tmp,int j,COMPLEX data[XYZnodeSites]){

    char fname[200]={0};
    sprintf(fname,"%s/OmgOmgwave_PH1.+%03d+%03d.%s-%06d",out_dir_path.c_str(),it,tshift[tmp],base,j);//hachi
	write_data(fname,data);
    cout<<"it = "<<it<<" Tshift count = "<<tmp<<"  conf =  "<<j<<"writing"<<endl;

}
//----------------------------------------------
//write data
//----------------------------------------------
int write_data(char fname[300],COMPLEX Wave_ave[XYZnodeSites]){
	fstream outfile;

	outfile.open(fname,ios::out|ios::binary|ios::trunc);
	if (!outfile) {
        cout << "ERR	output file can't open"<<endl;
        return 1;
    }
	for (int ix=0; ix<XnodeSites; ix++) {
	for (int iy=0; iy<YnodeSites; iy++) {
	for (int iz=0; iz<ZnodeSites; iz++) {
			outfile.write((const char * ) &wave_ave(ix,iy,iz),sizeof( COMPLEX ) );
			//cout << ix<<iy<<iz<<"	"<<wave_ave(ix,iy,iz)<<endl;
	}}}
		return 0;
}
//----------------------------------------------
//display matrix for debug
//----------------------------------------------
int display(Matrix &a){
cout << a[1][1]<<"	"<<a[1][2]<<"	"<<a[1][3]<<endl;
cout << a[2][1]<<"	"<<a[2][2]<<"	"<<a[2][3]<<endl;
cout << a[3][1]<<"	"<<a[3][2]<<"	"<<a[3][3]<<endl;
cout << endl;
	return 0;
}
int display(Matrix* &b){
	Matrix a=b[0];
	cout << a[1][1]<<"	"<<a[1][2]<<"	"<<a[1][3]<<endl;
	cout << a[2][1]<<"	"<<a[2][2]<<"	"<<a[2][3]<<endl;
	cout << a[3][1]<<"	"<<a[3][2]<<"	"<<a[3][3]<<endl;
	cout << endl;
	return 0;
}
//----------------------------------------------
//map xyz==>xPyPzP
//----------------------------------------------
void map(Matrix* &a,COMPLEX Wave_in[],COMPLEX Wave_out[]){
	Matrix ALL=a[0];
	int ixP=0;
	int iyP=0;
	int izP=0;

	for (int iz=0; iz<ZnodeSites; iz++) {
	for (int iy=0; iy<YnodeSites; iy++) {
	for (int ix=0; ix<XnodeSites; ix++) {

		ixP=ALL[1][1]*ix+ALL[1][2]*iy+ALL[1][3]*iz;
		iyP=ALL[2][1]*ix+ALL[2][2]*iy+ALL[2][3]*iz;
		izP=ALL[3][1]*ix+ALL[3][2]*iy+ALL[3][3]*iz;

		wave_out(ixP,iyP,izP)=wave_in(ix,iy,iz);
	//denug
		//cout <<"in "<<ix<<"	"<<iy<<"	"<<iz<<"	" <<wave_in(ix,iy,iz)<<endl;
		//cout <<"out"<<ixP<<"	"<<iyP<<"	"<<izP<<"	" <<wave_out(ixP,iyP,izP)<<endl;
		assert(wave_out(ixP,iyP,izP).real()!=0 && wave_out(ixP,iyP,izP).imag()!=0);
	}}}
	//debug end
}
//--------------------------------------------------------------------------
// brief check endian
//--------------------------------------------------------------------------
static bool machine_is_little_endian()
{
  int int0=0;
  *((char*) &int0)=1;
  if (int0==1)
    return true;
  else
    return false;
}

//----------------------------------------------
//endian convert
//----------------------------------------------
static void endian_convert(double* buf,int len)
{
 	for(int i=0; i<len; i++){
		char tmp[8];
		((double*)tmp)[0] = *buf;
		for(int j=0; j<8; j++){
			((char*)buf)[j] = ((char*)tmp)[7-j];
		}
		buf++;
	}
}
//----------------------------------------------
//matrix check
//----------------------------------------------
int check(Matrix* &a){
	Matrix ALL=a[0];
	static int countpxx=0;
	static int countmxx=0;
	static int countpxy=0;
	static int countmxy=0;
	static int countpxz=0;
	static int countmxz=0;

	static int countpyx=0;
	static int countmyx=0;
	static int countpyy=0;
	static int countmyy=0;
	static int countpyz=0;
	static int countmyz=0;

	static int countpzx=0;
	static int countmzx=0;
	static int countpzy=0;
	static int countmzy=0;
	static int countpzz=0;
	static int countmzz=0;

	if(ALL[1][1]==1) countpxx++;
	cout <<"x==>+x	::"<<countpxx<<endl;
	if(ALL[1][1]==-1) countmxx++;
	cout <<"x==>-x	::"<<countmxx<<endl;
	if(ALL[2][1]==1) countpxy++;
	cout <<"x==>+y	::"<<countpxy<<endl;
	if(ALL[2][1]==-1) countmxy++;
	cout <<"x==>-y	::"<<countmxy<<endl;
	if(ALL[3][1]==1) countpxz++;
	cout <<"x==>+z	::"<<countpxz<<endl;
	if(ALL[3][1]==-1) countmxz++;
	cout <<"x==>-z	::"<<countmxz<<endl;

	if(ALL[1][2]==1) countpyx++;
	cout <<"y==>+x	::"<<countpyx<<endl;
	if(ALL[1][2]==-1) countmyx++;
	cout <<"y==>-x	::"<<countmyx<<endl;
	if(ALL[2][2]==1) countpyy++;
	cout <<"y==>+y	::"<<countpyy<<endl;
	if(ALL[2][2]==-1) countmyy++;
	cout <<"y==>-y	::"<<countmyy<<endl;
	if(ALL[3][2]==1) countpyz++;
	cout <<"y==>+z	::"<<countpyz<<endl;
	if(ALL[3][2]==-1) countmyz++;
	cout <<"y==>-z	::"<<countmyz<<endl;

	if(ALL[1][3]==1) countpzx++;
	cout <<"z==>+x	::"<<countpzx<<endl;
	if(ALL[1][3]==-1) countmzx++;
	cout <<"z==>-x	::"<<countmzx<<endl;
	if(ALL[2][3]==1) countpzy++;
	cout <<"z==>+y	::"<<countpzy<<endl;
	if(ALL[2][3]==-1) countmzy++;
	cout <<"z==>-y	::"<<countmzy<<endl;
	if(ALL[3][3]==1) countpzz++;
	cout <<"z==>+z	::"<<countpzz<<endl;
	if(ALL[3][3]==-1) countmzz++;
	cout <<"z==>-z	::"<<countmzz<<endl;

	//display(ALL);
	return 0;

}
