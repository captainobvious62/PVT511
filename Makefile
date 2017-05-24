# Makefile for PTE 511 Project

FFLAGS =
FPPFLAGS =
LDLIBS =

-include ${PETSC_DIR}/conf/variables
-include ${PETSC_DIR}/conf/rules
-include ${PETSC_DIR}/lib/petsc/conf/variables
-include ${PETSC_DIR}/lib/petsc/conf/rules

OBJS = main.o m_global.o m_roots.o

m_roots.o   : m_roots.F90
m_global.o  : m_global.F90 m_roots.o
main.o      : main.F90 m_global.o
all: ${OBJS}
	-${FLINKER} ${OBJS} -o ./characterization ${PETSC_LIB} ${LDLIBS}
