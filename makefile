ALL: blib exec 

#compilation and various flags
DIRS    = 
EXEC    = generateur
CFLAGS	= 
FFLAGS	= 
CPPFLAGS	= 
FPPFLAGS	=
CLEANFILES  = ${EXEC}
OFILES= ${wildcard ./*.o}
NBLOCK = 2
DEBUG_VALGRIND = valgrind

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

blib :
	-@echo "BEGINNING TO COMPILE LIBRARIES "
	-@echo "========================================="
	-@${OMAKE}  PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} ACTION=libfast tree
	-@echo "Completed building libraries"
	-@echo "========================================="

distclean :
	-@echo "Cleaning application and libraries"
	-@echo "========================================="
	-@${OMAKE} PETSC_ARCH=${PETSC_ARCH}  PETSC_DIR=${PETSC_DIR} clean
	-${RM} ${OFILES}
	-@echo "Finised cleaning application and libraries"
	-@echo "========================================="	

exec: edmg.o
	-@echo "BEGINNING TO COMPILE APPLICATION "
	-@echo "========================================="
	-@${CLINKER} -o ${EXEC} edmg.o ${PETSC_LIB}
	-@echo "Completed building application"
	-@echo "========================================="

run1:
	-@${MPIEXEC} -np 1 ./generateur -vfile lsqr.bin -n 274 -nzeros 100 ${DEBUG_VALGRIND}
run2:
<<<<<<< HEAD
	-@${MPIEXEC} -np 20 ./generateur -n 1000000 -nzeros 999999
=======
	-@${MPIEXEC} -np 2 ./generateur -n 300 -nzeros 200
>>>>>>> 8354f04d89e19bf96d92ea0c1328fb1e3cb3df5c
run3:
	-@${MPIEXEC} -np 2 ./generateur -n 300 -nzeros 200 -realMat
run4:
	-@${MPIEXEC} -np 1 ./generateur -vfile lsqr.bin -n 274 -nzeros 200 -realMat
<<<<<<< HEAD

run5:
	./generateur -n 330001 -nzeros 330000
=======
>>>>>>> 8354f04d89e19bf96d92ea0c1328fb1e3cb3df5c
