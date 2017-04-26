CC=gcc
UMFPACK_FLAGS=-I/user/include/suitesparse -lumfpack
CFLAGS=-Wall

TESTS: solver_test

solver_test: solver_test.c local_solvers.c
	${CC} ${CFLAGS} ${UMFPACK_FLAGS} solver_test.c -o solver_test
