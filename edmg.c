/*Copyright (c) 2017. Xinzhe WU in Maison de la Simulation. All rights reserved */

#include "edmg.h"

//#define CLOCKS_PER_SEC ((clock_t)1000) 
#define PI 3.1415926
#define epsilon 1

#define max(a,b) (a>=b?a:b)
#define min(a,b) (a<=b?a:b)

static char help[] = "Dense Matrix Generator by select eigenvalues";

PetscErrorCode getFileSize(const char * name, long * size){
  FILE * fptr;
  *size = 0L;

#ifdef LINUX
  struct stat fs;

  if(stat(name,&fs)!=0){
    perror("Cannot state file\n");
  }
  *size=fs.st_size;

#else
fptr=fopen(name,"rb");
  if(fptr!=NULL){
    fseek(fptr,0L,SEEK_END);
    *size = ftell(fptr);
    fclose(fptr);
  }
#endif

  return 0;
}

PetscErrorCode readBinaryScalarArray(const char * name, int * nb, PetscScalar * array){
  int file_descriptor;
  PetscErrorCode ierr;
  long size;

  getFileSize(name,&size);

  if(*nb<=0) *nb=(int)size/((int)sizeof(PetscScalar));
  if(size/sizeof(PetscScalar)!=*nb) {
    return 1;
  }


  ierr=PetscBinaryOpen(name,FILE_MODE_READ,&file_descriptor);CHKERRQ(ierr);
  ierr=PetscBinarySynchronizedRead(PETSC_COMM_WORLD,file_descriptor,array,*nb,PETSC_SCALAR);CHKERRQ(ierr);
  ierr=PetscBinaryClose(file_descriptor);CHKERRQ(ierr);

  return ierr;
}


void random_selection(PetscScalar *ret, PetscInt nombre)
{

  PetscInt		i;
  PetscScalar 	flg;

  int 			my_seed,my_rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&my_rank);
  my_seed=time(NULL)+my_rank;

  srand(my_seed); 

  if(rand()/(double)RAND_MAX > 0.5){
  	flg = 1.0;
  }
  else{
  	flg = 0.0;
  }


  for(i = 0; i < nombre; i++){

    if(i < nombre / 2)
    	ret[i] = (cos(2*PI*i/nombre)+2) + PETSC_i*sin(2*PI*i/nombre);
    else ret[i] = 2 * (cos(2*PI*i/nombre)+2) + PETSC_i*sin(2*PI*i/nombre);
  }
}

void shuffer(PetscScalar *array, PetscInt n)
{
	int index, i;
	PetscScalar tmp;
	srand(time(NULL));

	for(i = n-1; i>0;i--){
		index = rand() % i;
		tmp = array[i];
		array[i]=array[index];
		array[index]=tmp;
	}
}

int *indexShuffer(PetscInt n)
{
	int index, i;
	int *a;
	a = (int *) malloc(n*sizeof(int));
	PetscInt tmp;
	srand(time(NULL));
	for (i = 0; i < n; i++)
		a[i] = i;

	for(i = n-1; i>0;i--){
		index = rand() % i;
		tmp = a[i];
		a[i]=a[index];
		a[index]=tmp;
	}
	return a;
}

void printarray(PetscInt n, PetscScalar *a) {

	PetscScalarView(n, a, PETSC_VIEWER_STDOUT_WORLD);
}

int main(int argc, char ** argv){

	PetscInitialize(&argc,&argv,(char *)0,help);

	char fileb[PETSC_MAX_PATH_LEN];
  
  	Mat            A;

  	PetscErrorCode ierr;

  	PetscInt       n, k, nzeros;

  	PetscInt    	size;

  	PetscBool      flagb, flagn, flagnzeros;

	MatInfo     	Ainfo;
	double        	gnnz;
	char           	matrixOutputFile[PETSC_MAX_PATH_LEN];
	PetscViewer    	output_viewer;

	int             world_size;

	MPI_Comm_size(PETSC_COMM_WORLD, &world_size);
        PetscPrintf(PETSC_COMM_WORLD, "--------------------------\n");
	PetscPrintf(PETSC_COMM_WORLD, "Using number of %d precessors for the generation...\n", world_size);

    ierr=PetscOptionsGetInt(PETSC_NULL,"-n",&n,&flagn);CHKERRQ(ierr);
    if (!flagn){
		ierr = PetscPrintf(PETSC_COMM_WORLD,"!!!Please set the dimension of matrix to be generated\n");CHKERRQ(ierr);
		return 0;
	}
    ierr=PetscOptionsGetInt(PETSC_NULL,"-nzeros",&nzeros,&flagnzeros);CHKERRQ(ierr);
    if (!flagnzeros){
		ierr = PetscPrintf(PETSC_COMM_WORLD,"!!!Please number of zeros per row of matrix to be generated\n");CHKERRQ(ierr);
		return 0;
	}
        PetscPrintf(PETSC_COMM_WORLD, "To generate matrix with dim = %d, number zeros = %d ... \n", n, nzeros);

  	int *Apermut;
  	Apermut = indexShuffer(n);

  	PetscReal *column_condensed;
	PetscReal *row_condensed;

	PetscMalloc1(n-nzeros, &column_condensed);
	PetscMalloc1(n-nzeros, &row_condensed);

	for (k = 0; k < (n-nzeros); k++) {
 		column_condensed[k]=(PetscReal)rand()/RAND_MAX;
 		row_condensed[k]=(PetscReal)rand()/RAND_MAX;
	}

  	int *RCpermut;
  	RCpermut = indexShuffer(n);

  	PetscScalar *Deigenvalues;
  	PetscMalloc1(n,&Deigenvalues);
        PetscPrintf(PETSC_COMM_WORLD, "--------------------------\n");
    ierr=PetscOptionsGetString(PETSC_NULL,"-vfile",fileb,PETSC_MAX_PATH_LEN-1,&flagb);CHKERRQ(ierr);
	
	if (!flagb){
		PetscPrintf(PETSC_COMM_WORLD, "Not providing the outside eigenvalues files, using the internal functions to generate them...\n");
		random_selection(Deigenvalues,n);
		shuffer(Deigenvalues,n);
	}
	else{
        	PetscPrintf(PETSC_COMM_WORLD, "Using the eigenvalues provides by outside files: %s ...\n", fileb);
		readBinaryScalarArray(fileb, &size, Deigenvalues);
		shuffer(Deigenvalues,size);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"read file size = %d\n", size);CHKERRQ(ierr);
		if(size != n){
			ierr = PetscPrintf(PETSC_COMM_WORLD,"!!!read file size and vec dimemson do not match and mat dim set to be equal to vec dim\n");CHKERRQ(ierr);
			return 0;
		}
		n = size;
	}

//	printarray(n,Deigenvalues);
        PetscPrintf(PETSC_COMM_WORLD, "--------------------------\n");
        PetscPrintf(PETSC_COMM_WORLD, "@>Generating ...\n");
	PetscReal RinvAC=0;
	PetscReal RinvADC=0;

	PetscReal access_column, access_row;
	PetscInt k1;
	PetscInt k2;

	for (k1 = 0; k1 < n; k1++) {
		if (RCpermut[k1]>=(n-nzeros))
			access_row=0;
		else
			access_row=row_condensed[RCpermut[k1]];

		if (RCpermut[Apermut[k1]]>=(n-nzeros))
			access_column=0;
		else
			access_column=column_condensed[RCpermut[Apermut[k1]]];

		RinvAC+=access_column*access_row;
		RinvADC+=Deigenvalues[Apermut[k1]]*access_column*access_row;
	}

	PetscReal inv_lambda;

	inv_lambda=-RinvAC+epsilon;
	clock_t start, finish;
	double  duration;	
	PetscInt Istart, Iend;
	start = clock();
	ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
	ierr = MatSetType(A,MATMPIAIJ);CHKERRQ(ierr);
	ierr = MatSetFromOptions(A);CHKERRQ(ierr); 
	ierr = MatSetUp(A);CHKERRQ(ierr);
	MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	//ierr = MatMPIAIJSetPreallocation(A,n-nzeros,NULL,n-nzeros,NULL);CHKERRQ(ierr);
	//ierr = MatSeqAIJSetPreallocation(A,n-nzeros,NULL);
	ierr = MatGetOwnershipRange(A,&Istart,&Iend); CHKERRQ(ierr);

	PetscScalar matrix_element=.0;

	int cnt = 0;
	int index[n];
	PetscScalar val[n];
  	 for (k1 = Istart; k1 < Iend; k1++) {
		for(k2=0; k2 < n; k2++){
			if(k2 == k1){
				matrix_element = Deigenvalues[Apermut[k1]];
                                MatSetValues(A,1,&k1,1,&k2,&matrix_element,INSERT_VALUES); 
			}
			
			if(RCpermut[k2]<(n-nzeros) && RCpermut[Apermut[k1]] < (n-nzeros)){
                                {
				access_row=row_condensed[RCpermut[k2]];
				access_column=column_condensed[RCpermut[Apermut[k1]]];
				index[cnt]=k2;
				val[cnt] = matrix_element + (1/(inv_lambda))*Deigenvalues[Apermut[k1]]*access_column*access_row 
					   -(1/epsilon)*Deigenvalues[Apermut[k2]]*access_column*access_row 
					   -(RinvADC/(inv_lambda*epsilon))*access_column*access_row;
///			 	printf("k1 = %d ; k2 = %d, cnt = %d, index[%d] = %d \n", k1, k2, cnt, cnt, index[cnt]);		
				cnt++;
				}
		                MatSetValues(A,1,&k1,cnt,index,val,INSERT_VALUES);	
		} 
  	}
	cnt = 0;
}

	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	MatGetInfo(A,MAT_GLOBAL_SUM,&Ainfo);
	gnnz = Ainfo.nz_used;

	sprintf(matrixOutputFile,"Eigen_known_matrix_nb_%d_%dx%d_%g_nnz.gz",n, n,n,gnnz);
		
	PetscPrintf(PETSC_COMM_WORLD,"\n@>Dumping matrix to PETSc binary %s\n",matrixOutputFile);
			
	PetscViewerBinaryOpen(PETSC_COMM_WORLD,matrixOutputFile,FILE_MODE_WRITE,&output_viewer);
	PetscViewerSetFormat(output_viewer,PETSC_VIEWER_ASCII_INFO_DETAIL);
///	MatView(A,PETSC_VIEWER_STDOUT_WORLD);
	MatView(A,output_viewer);
	PetscViewerDestroy(&output_viewer);
		
	PetscPrintf(PETSC_COMM_WORLD,"\n@>Matrix %s Dumped\n\n",matrixOutputFile);
	PetscPrintf(PETSC_COMM_WORLD,"\n>>>>>>Please use the command 'gzip -d **' to unzip the file to binary file\n\n");

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
        PetscPrintf(PETSC_COMM_WORLD, "--------------------------\n");
        PetscPrintf(PETSC_COMM_WORLD,"\nElapsed time is %f seconds\n\n", duration);
        PetscPrintf(PETSC_COMM_WORLD, "--------------------------\n");
	PetscFinalize();

	return 0;
}
