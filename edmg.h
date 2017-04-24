/*Copyright (c) 2017. Xinzhe WU in Maison de la Simulation. All rights reserved */

#ifndef _EDMG_H
#define _EDMG_H

#include "petscmat.h"
#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h> 
#include <petscvec.h>
#include <petscksp.h>
#include "petsc.h"


typedef struct _MatrixInfo{
	int n;
	int m;
	int nnz;
} MatrixInfo;

#endif
