/*
 *  noise_redection0.1.h
 *  
 *
 *  Created by yamada on 2/4/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "./matrix.h"
#include <complex>

#ifndef INCLUDED_NOISEREDUCTION_H
#define INCLUDED_NOISEREDUCTION_H
#include <cassert>

using namespace std;
#define I       std::complex<double>(0.0,1.0)
typedef std::complex<double> COMPLEX;


//static const int Xsites=16;
//static const int Ysites=16;
//static const int Zsites=16;
static const int Matrix_size=3;
//static const int T_i=6;
//static const int T_f=10;

void call_data(int,int,int,COMPLEX[]);
int read_data(char[],COMPLEX[]);
void out_data(int,int,int,COMPLEX[]);
int write_data(char[],COMPLEX[]);
static void endian_convert(double* ,int);
static bool machine_is_little_endian();

int setC4();
void make_rep();
int display(Matrix&);
int display(Matrix* &);
void sum();
void map(Matrix *&,COMPLEX[],COMPLEX[]);
int check(Matrix *&);

#define wave_in(x,y,z) Wave_in[(x) +XnodeSites*((y) + YnodeSites*((z)))]
#define wave_out(x,y,z) Wave_out[((x+XnodeSites)%XnodeSites) +XnodeSites*(((y+YnodeSites)%YnodeSites) + YnodeSites*(((z+ZnodeSites)%ZnodeSites)))]
#define wave_ave(x,y,z) Wave_ave[(x) +XnodeSites*((y) + YnodeSites*((z)))]

#endif
