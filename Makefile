CC=gcc
LIB_FLAGS=-lumfpack -lcxsparse
CFLAGS=-Wall

ALL: TESTS

NOTEST:

TESTS: solver_test diag_test

clean:
	rm -f solver_test diag_test 

solver_test: solver_test.c local_solvers.c
	${CC} ${CFLAGS} ${LIB_FLAGS} solver_test.c -o solver_test; ./solver_test

diag_test: diag_test.c local_solvers.c
	${CC} ${CFLAGS} ${LIB_FLAGS} diag_test.c -o diag_test; ./diag_test
