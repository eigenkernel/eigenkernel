-include ../../../Makefile.inc
-include ../Makefile.in

EIGENEXA_WITH_TIMER_FLAG = $(MACRO_PREFIX)EIGENEXA_WITH_TIMER=0
#EIGENEXA_WITH_TIMER_FLAG = $(MACRO_PREFIX)EIGENEXA_WITH_TIMER=1

# Remove period in a version string. Supported: 2014.06.001, 2015.02.001. 2015.02.002, 2015.05.001
ifndef ELPA_VERSION
	ELPA_VERSION = 201505001
endif
ELPA_VERSION_FLAG = $(MACRO_PREFIX)ELPA_VERSION=$(ELPA_VERSION)

FFLAGS_INTERNAL = $(FFLAGS)
LDFLAGS_INTERNAL = $(LDFLAGS)

TARGET_MAIN = eigbench

TARGET_LIB = libeigenkernel.a

.PHONY: all lib clean dep

all: $(TARGET_MAIN)

lib: $(TARGET_LIB)

OBJS_MAIN = main.o

OBJS_SOLVER = solver/solver_main.o solver/solver_lapack.o solver/solver_scalapack_all.o solver/solver_scalapack_select.o solver/verifier.o solver/generalized_to_standard.o

OBJS_EIGENEXA = solver/solver_eigenexa.o
ifeq ($(WITH_EIGENEXA), 1)
	SRCS_EIGENEXA = solver/solver_eigenexa.f90
	FFLAGS_EIGENEXA = $(FFLAGS_INTERNAL) $(CPPFLAG) $(EIGENEXA_WITH_TIMER_FLAG)
else
	SRCS_EIGENEXA = solver/solver_eigenexa_dummy.f90
	FFLAGS_EIGENEXA = $(FFLAGS_INTERNAL) -Wno-unused-dummy-argument
endif

OBJS_ELPA = solver/solver_elpa.o
ifeq ($(WITH_ELPA), 1)
	SRCS_ELPA = solver/solver_elpa.f90
	FFLAGS_ELPA = $(FFLAGS_INTERNAL) $(CPPFLAG) $(ELPA_VERSION_FLAG)
else
	SRCS_ELPA = solver/solver_elpa_dummy.f90
	FFLAGS_ELPA = $(FFLAGS_INTERNAL) -Wno-unused-dummy-argument
endif

OBJS_ELPA_EIGENEXA = solver/solver_elpa_eigenexa.o
ifeq ($(WITH_ELPA), 1)
	ifeq ($(WITH_EIGENEXA), 1)
		SRCS_ELPA_EIGENEXA = solver/solver_elpa_eigenexa.f90
		FFLAGS_ELPA_EIGENEXA = $(FFLAGS_INTERNAL) $(CPPFLAG) $(ELPA_VERSION_FLAG)
	else
		SRCS_ELPA_EIGENEXA = solver/solver_elpa_eigenexa_dummy.f90
		FFLAGS_ELPA_EIGENEXA = $(FFLAGS_INTERNAL) -Wno-unused-dummy-argument
	endif
else
	SRCS_ELPA_EIGENEXA = solver/solver_elpa_eigenexa_dummy.f90
	FFLAGS_ELPA_EIGENEXA = $(FFLAGS_INTERNAL) -Wno-unused-dummy-argument
endif

OBJS_UTIL_F77 = util/mmio.o

OBJS_UTIL = util/descriptor_parameters.o util/global_variables.o util/command_argument.o util/matrix_io.o util/distribute_matrix.o util/processes.o util/eigenpairs_types.o util/event_logger.o util/fson.o util/modules.o

OBJS_LIB = $(OBJS_SOLVER) $(OBJS_EIGENEXA) $(OBJS_ELPA) $(OBJS_ELPA_EIGENEXA) $(OBJS_UTIL) $(OBJS_UTIL_F77)

OBJS = $(OBJS_MAIN) $(OBJS_LIB)

$(OBJS_UTIL_F77): %.o: %.f
	$(FC) -c $(FFLAGS_INTERNAL) $(LIBS) $< -o $@

$(OBJS_MAIN) $(OBJS_SOLVER) $(OBJS_UTIL): %.o: %.f90
	$(FC) -c $(FFLAGS_INTERNAL) $(LIBS) $< -o $@

$(OBJS_EIGENEXA): $(SRCS_EIGENEXA)
	$(FC) -c $(FFLAGS_EIGENEXA) $(LIBS) $< -o $@

$(OBJS_ELPA): $(SRCS_ELPA)
	$(FC) -c $(FFLAGS_ELPA) $(LIBS) $< -o $@

$(OBJS_ELPA_EIGENEXA): $(SRCS_ELPA_EIGENEXA)
	$(FC) -c $(FFLAGS_ELPA_EIGENEXA) $(LIBS) $< -o $@

$(TARGET_MAIN): $(OBJS)
	$(FC) $(LDFLAGS_INTERNAL) $^ $(LIBS) -o $@
	@mkdir -p ../bin
	cp $(TARGET_MAIN) ../bin/$(TARGET_MAIN)

$(TARGET_LIB): $(OBJS_LIB)
	$(AR) r $(TARGET_LIB) $^

clean:
	@rm -f *.exe *.o *.mod *~ ../bin/* */*.o */*.mod */*~ $(TARGET_MAIN) $(TARGET_LIB)

dep:
	find . -name \*.f90 -and ! -name \*_dummy.f90 | xargs makedepf90 -nosrc > Makefile.dep

include Makefile.dep
