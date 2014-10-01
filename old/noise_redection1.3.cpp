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



int main(){
	
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
	
	
	for (int it=T_i; it<T_f+1; it++) {
	for (int j=0; j<Confsize; j++) {
	//for (int id=0; id<24; id++) {
	//	cout << "Matrix"<<id<<endl;
	//	display(MA[id]);
	//	check(MA[id]);
	//}
	COMPLEX *Wave_in=new COMPLEX[XYZsites];
	COMPLEX *Wave_out=new COMPLEX[XYZsites];
	COMPLEX *Wave_ave=new COMPLEX[XYZsites];
	call_data(it,j,Wave_in);
	cout << "map	start"<<endl;
	for (int id=0; id<24; id++) {
		map(MA[id],Wave_in,Wave_out);
	for (int ix=0; ix<Xsites; ix++) {
	for (int iy=0; iy<Ysites; iy++) {
	for (int iz=0; iz<Zsites; iz++) {
		wave_ave(ix,iy,iz)=wave_out(ix,iy,iz)/24.0+wave_ave(ix,iy,iz);
	}}}}
	
	cout << "map	end"<<endl;
		out_data(it,j,Wave_ave);
	delete []Wave_ave;
	delete []Wave_out;
	delete []Wave_in;
	}}
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
void call_data(int it,int b,COMPLEX data[XYZsites]){
	
	char xyz_fname[300]={0};
		//sprintf(xyz_fname,"./results/xyz/OmgOmgwave_PH1.+%03d+%03d.RC16x32_B1830Kud013760Ks013710C1761-1-%05d0",it,tshift,j);//lab 
		//sprintf(xyz_fname,"/Users/SINYAMADA/lab/src/OmgOmg4point/results/Proj-OmgOmg_bibary/xyz/s%dz0s%dz0.RC16x32_B1830Kud013760Ks013710C1761-1-%05d0.kap_013710.000_000_000_000_it%02d",0,0,j,it);//test 4^3
		//sprintf(xyz_fname,"/Users/SINYAMADA/lab/src/noise_reduction/results/read/OmgOmgRwave.004200.+006.RC16x32_B1830Kud013760Ks013710C1761-000000");//test 16^3
	sprintf(xyz_fname,"/work/data2/yamada/OmgOmg/results/Rcor/binR/xyz/OmgOmgRwave.%06d.+%03d.RC16x32_B1830Kud013760Ks013710C1761-%05d0",Confsize,it,b);

	
	read_data(xyz_fname,data);
}
//----------------------------------------------
//read data
//----------------------------------------------
int read_data(char fname[300],COMPLEX data[XYZsites]){
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
	//endian_convert((double*)data,XYZsites*2);
	for (int point=0; point<id; point++) {
		//data[point]=rand()%24;
		//cout << data[point]<<endl;                                                                                                                                                          
	}
	return 0;
}

//----------------------------------------------
//output file name
//----------------------------------------------
void out_data(int it,int b,COMPLEX data[XYZsites]){
	
	char xyz_fname[300]={0};
	//sprintf(xyz_fname,"/work/data2/yamada/OmgOmg/results/xyz/OmgOmgwave_PH1.+%03d+%03d.RC16x32_B1830Kud013760Ks013710C1761-1-%05d0",it,tshift,j); 
	//sprintf(xyz_fname,"/work/data2/yamada/OmgOmg/results/xyz/OmgOmgwave_PH1.+%03d+%03d.RC16x32_B1830Kud013760Ks013710C1761-1-%05d0",it,tshift,j); 
	//sprintf(xyz_fname,"../results/xyz/impOmgOmgwave_PH1.+%03d+%03d.RC16x32_B1830Kud013760Ks013710C1761-1-%05d0",it,tshift,j); 
	sprintf(xyz_fname,"/work/data2/yamada/OmgOmg/results/Rcor/redbin/xyz/OmgOmgRwave.%06d.+%03d.RC16x32_B1830Kud013760Ks013710C1761-%05d0",Confsize,it,b);
	write_data(xyz_fname,data);
}
//----------------------------------------------
//write data
//----------------------------------------------
int write_data(char fname[300],COMPLEX Wave_ave[XYZsites]){
	fstream outfile;
	
	outfile.open(fname,ios::out|ios::binary|ios::trunc);
	if (!outfile) {
        cout << "ERR	output file can't open"<<endl;
        return 1;
    }
	for (int ix=0; ix<Xsites; ix++) {
	for (int iy=0; iy<Ysites; iy++) {
	for (int iz=0; iz<Xsites; iz++) {		
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

	for (int iz=0; iz<Xsites; iz++) {
	for (int iy=0; iy<Ysites; iy++) {
	for (int ix=0; ix<Xsites; ix++) {	
		
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
