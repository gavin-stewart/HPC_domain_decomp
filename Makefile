CC=gcc
MPICC=mpicc
LIB_FLAGS=-lumfpack -lcxsparse -lcholmod
CFLAGS=-Wall -Wno-unused-function

.PHONY: ALL NOTEST TEST clean

ALL: NOTEST TESTS

NOTEST: domain_decomp

TESTS: umfpack_test diag_test 

clean:
	rm -f umfpack_test diag_test cholmod_test domain_decomp

domain_decomp: domain_decomp.c local_solvers.c residue.c
	${MPICC} ${CFLAGS} ${LIB_FLAGS} domain_decomp.c -o domain_decomp

umfpack_test: umfpack_test.c local_solvers.c
	${CC} ${CFLAGS} ${LIB_FLAGS} umfpack_test.c -o umfpack_test; ./umfpack_test

diag_test: diag_test.c local_solvers.c
	${CC} ${CFLAGS} ${LIB_FLAGS} diag_test.c -o diag_test; ./diag_test

#cholmod_test: cholmod_test.c local_solvers.c
#	${CC} ${CFLAGS} ${LIB_FLAGS} cholmod_test.c -o cholmod_test; ./cholmod_test
