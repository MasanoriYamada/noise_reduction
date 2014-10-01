/*
 *  Rep.h
 *  
 *
 *  Created by yamada on 2/4/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef INCLUDED_REP_H
#define INCLUDED_REP_H

#include "./Matrix.h"
#include "./noise_redection.h"

class Rep : public Matrix{
public:
	
//C1 1
	Matrix EMatrix_size,Matrix_size;
//C2 9	
	Matrix C2_1Matrix_size,Matrix_size;
	Matrix C2_2Matrix_size,Matrix_size;
	Matrix C2_3Matrix_size,Matrix_size;
	Matrix C2_4Matrix_size,Matrix_size;
	Matrix C2_5Matrix_size,Matrix_size;
	Matrix C2_6Matrix_size,Matrix_size;
	Matrix C2_7Matrix_size,Matrix_size;
	Matrix C2_8Matrix_size,Matrix_size;
	Matrix C2_9Matrix_size,Matrix_size;

//C3 8	
	Matrix C3_1Matrix_size,Matrix_size;
	Matrix C3_2Matrix_size,Matrix_size;
	Matrix C3_3Matrix_size,Matrix_size;
	Matrix C3_4Matrix_size,Matrix_size;
	Matrix C3_5Matrix_size,Matrix_size;
	Matrix C3_6Matrix_size,Matrix_size;
	Matrix C3_7Matrix_size,Matrix_size;
	Matrix C3_8Matrix_size,Matrix_size;

//C4 6	
	Matrix C4_1Matrix_size,Matrix_size;
	Matrix C4_2Matrix_size,Matrix_size;
	Matrix C4_3Matrix_size,Matrix_size;
	Matrix C4_4Matrix_size,Matrix_size;
	Matrix C4_5Matrix_size,Matrix_size;
	Matrix C4_6Matrix_size,Matrix_size;

	Rep();	//constract
	int set();
	
};
#endif