/*Copyright (c) 2011—2017. Xinzhe WU in Maison de la Simulation. All rights reserved */

#include "gen.h"

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


void random_selection(PetscReal *ret, PetscInt nombre)
{

  PetscInt		i, indice;
  int 			my_seed,my_rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&my_rank);
  my_seed=time(NULL)+my_rank;

  srand(my_seed); 
  for(i = 0; i < nombre; i++){
    ret[i] = rand()/(double)RAND_MAX;
  }
}

void selection(PetscScalar *ret, PetscInt nombre, PetscInt min, PetscInt max)
{

  PetscScalar		*tab;
  PetscInt		i, indice, maxi = max - min;
 
  if(min >= max || nombre > maxi + 1 || nombre < 1)
    PetscPrintf(PETSC_COMM_WORLD,"Input values of tirage() are wrong.\n");
 
  PetscMalloc((maxi + 1)*sizeof(tab[0]),&tab);

  for(i = 0; i < maxi + 1; i++)
    tab[i] = i + min;
 
  for(i = 0; i < nombre; i++){
    indice = rand() % (maxi + 1);
    ret[i] = tab[indice];
    tab[indice] = tab[maxi];
    maxi--;
  }
  PetscFree(tab);
}

void shuffer(PetscReal *array, PetscInt n)
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

void printarray(PetscInt n, PetscReal *a) {
    int k = 0;

    for (k = 0; k < n; k++) {
		printf("%e   ", a[k]);
		if (k % 8 == 7)
	    	printf("\n");
    }
}


PetscScalar Random (PetscInt _iMin, PetscInt _iMax) 
{ 
	return (_iMin + (rand () % (_iMax-_iMin+1))); 
} 

PetscInt IRandom (PetscInt _iMin, PetscInt _iMax) 
{ 
	return (_iMin + (rand () % (_iMax-_iMin+1))); 
} 

int main(int argc, char ** argv){

  	PetscInt sizen;

	PetscInitialize(&argc,&argv,(char *)0,help);


	char buf[PETSC_MAX_PATH_LEN];
	char fileb[PETSC_MAX_PATH_LEN];

	Vec            eigenvalues, v;
  	PetscInt       vl_start,vl_end;

  	Mat            Mt, AM, MA, matAop;
  
  	Mat            A;

  	PetscErrorCode ierr;

  	PetscInt       n,i,j,k,degree, nzeros;
  	PetscScalar    rRandom1, rRandom2;
  	PetscInt       iRandom, d1, d2;

  	PetscInt    	size;


  	PetscScalar    *array, *tmp_array, *read_array;
  	PetscBool      flagb;

  	PetscMPIInt    my_rank;
  	PetscViewer fd;

	MatInfo     	Ainfo;
	double        	gnnz;
	char           	matrixOutputFile[PETSC_MAX_PATH_LEN];
	PetscViewer    	output_viewer, viewer;


  	n = 274;
  	nzeros =270; 
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

    ierr=PetscOptionsGetString(PETSC_NULL,"-vfile",fileb,PETSC_MAX_PATH_LEN-1,&flagb);CHKERRQ(ierr);
	
	if (!flagb){
		random_selection(Deigenvalues,n);
	}
	else{
		readBinaryScalarArray(fileb, &size, Deigenvalues);
		shuffer(Deigenvalues,size);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"read file size = %d\n", size);CHKERRQ(ierr);
		if(size != n){
			ierr = PetscPrintf(PETSC_COMM_WORLD,"read file size and vec dimemson do not match\n");CHKERRQ(ierr);
			return 0;
		}
	}

	printarray(n,Deigenvalues);

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

	ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
	ierr = MatSetType(A,MATMPIAIJ);CHKERRQ(ierr);
	ierr = MatSetFromOptions(A);CHKERRQ(ierr); 
	ierr = MatSetUp(A);CHKERRQ(ierr);

	PetscScalar matrix_element;

	for (k1 = 0; k1 < n; k1++) {
		for (k2 = 0; k2 < n; k2++) {
			matrix_element=0;
			if (k1==k2)
				matrix_element+=Deigenvalues[Apermut[k1]];
			if (RCpermut[k2]>=(n-nzeros))
				access_row=0;
			else
				access_row=row_condensed[RCpermut[k2]];

			if (RCpermut[Apermut[k1]]>=(n-nzeros))
				access_column=0;
			else
				access_column=column_condensed[RCpermut[Apermut[k1]]];

			matrix_element+=(1/(inv_lambda))*Deigenvalues[Apermut[k1]]*access_column*access_row;
			matrix_element+=-(1/epsilon)*Deigenvalues[Apermut[k2]]*access_column*access_row;
			matrix_element+=-(RinvADC/(inv_lambda*epsilon))*access_column*access_row;

			if(matrix_element != 0){
				ierr = MatSetValues(A,1,&k1,1,&k2,&matrix_element,INSERT_VALUES);CHKERRQ(ierr);
			}
		}
	}

	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	MatGetInfo(A,MAT_GLOBAL_SUM,&Ainfo);
	gnnz = Ainfo.nz_used;

	sprintf(matrixOutputFile,"Eigen_known_matrix_nb_%d_%dx%d_%g_nnz.gz",n, n,n,gnnz);
		
	PetscPrintf(PETSC_COMM_WORLD,"\n@>Dumping matrix to PETSc binary %s\n",matrixOutputFile);
			
	PetscViewerBinaryOpen(PETSC_COMM_WORLD,matrixOutputFile,FILE_MODE_WRITE,&output_viewer);
	PetscViewerSetFormat(output_viewer,PETSC_VIEWER_ASCII_INFO_DETAIL);
	MatView(A,PETSC_VIEWER_STDOUT_WORLD);
	MatView(A,output_viewer);
	PetscViewerDestroy(&output_viewer);
		
	PetscPrintf(PETSC_COMM_WORLD,"\n@>Matrix %s Dumped\n\n",matrixOutputFile);
	PetscPrintf(PETSC_COMM_WORLD,"\n>>>>>>Please use the command 'gzip -d **' to unzip the file to binary file\n\n");
	PetscFinalize();

	return 0;
}
