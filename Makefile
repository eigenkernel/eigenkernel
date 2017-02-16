-include ./Makefile.inc

EIGENEXA_WITH_TIMER_FLAG = $(MACRO_PREFIX)USE_EIGENEXA_WITH_TIMER=0
#EIGENEXA_WITH_TIMER_FLAG = $(MACRO_PREFIX)USE_EIGENEXA_WITH_TIMER=1

# Remove period in a version string. Supported: 2014.06.001, 2015.02.001. 2015.02.002, 2015.05.001
ifndef ELPA_VERSION
	ELPA_VERSION = 201505001
endif
ELPA_VERSION_FLAG = $(MACRO_PREFIX)ELPA_VERSION=$(ELPA_VERSION)

FFLAGS_INTERNAL = $(FFLAGS)
LDFLAGS_INTERNAL = $(LDFLAGS)

TARGET_MAIN = bin/eigbench

TARGET_LIB = libeigenkernel.a

.PHONY: all main lib clean dep

all: main lib

main: $(TARGET_MAIN)

lib: $(TARGET_LIB)

OBJS_MAIN = src/main.o

OBJS_SOLVER = \
	src/solver_main.o \
	src/solver_lapack.o \
	src/solver_scalapack_all.o \
	src/solver_scalapack_select.o \
	src/verifier.o \
	src/generalized_to_standard.o

OBJS_EIGENEXA = src/solver_eigenexa.o
ifeq ($(WITH_EIGENEXA), 1)
	SRCS_EIGENEXA = src/solver_eigenexa.f90
	FFLAGS_EIGENEXA = $(FFLAGS_INTERNAL) $(CPPFLAG) $(EIGENEXA_WITH_TIMER_FLAG)
else
	SRCS_EIGENEXA = src/solver_eigenexa_dummy.f90
	FFLAGS_EIGENEXA = $(FFLAGS_INTERNAL) $(WNO_UNUSED_DUMMY_ARGUMENT_FLAG)
endif

OBJS_ELPA = src/solver_elpa.o
ifeq ($(WITH_ELPA), 1)
	SRCS_ELPA = src/solver_elpa.f90
	FFLAGS_ELPA = $(FFLAGS_INTERNAL) $(CPPFLAG) $(ELPA_VERSION_FLAG)
else
	SRCS_ELPA = src/solver_elpa_dummy.f90
	FFLAGS_ELPA = $(FFLAGS_INTERNAL) $(WNO_UNUSED_DUMMY_ARGUMENT_FLAG)
endif

OBJS_ELPA_EIGENEXA = src/solver_elpa_eigenexa.o
ifeq ($(WITH_ELPA), 1)
	ifeq ($(WITH_EIGENEXA), 1)
		SRCS_ELPA_EIGENEXA = src/solver_elpa_eigenexa.f90
		FFLAGS_ELPA_EIGENEXA = $(FFLAGS_INTERNAL) $(CPPFLAG) $(ELPA_VERSION_FLAG)
	else
		SRCS_ELPA_EIGENEXA = src/solver_elpa_eigenexa_dummy.f90
		FFLAGS_ELPA_EIGENEXA = $(FFLAGS_INTERNAL) $(WNO_UNUSED_DUMMY_ARGUMENT_FLAG)
	endif
else
	SRCS_ELPA_EIGENEXA = src/solver_elpa_eigenexa_dummy.f90
	FFLAGS_ELPA_EIGENEXA = $(FFLAGS_INTERNAL) $(WNO_UNUSED_DUMMY_ARGUMENT_FLAG)
endif

OBJS_UTIL_F77 = src/mmio.o

OBJS_UTIL = \
	src/descriptor_parameters.o \
	src/global_variables.o \
	src/command_argument.o \
	src/matrix_io.o \
	src/distribute_matrix.o \
	src/processes.o \
	src/eigenpairs_types.o \
	src/event_logger.o \
	src/fson.o \
	src/modules.o

OBJS_LIB = $(OBJS_SOLVER) $(OBJS_EIGENEXA) $(OBJS_ELPA) $(OBJS_ELPA_EIGENEXA) $(OBJS_UTIL) $(OBJS_UTIL_F77)

OBJS = $(OBJS_MAIN) $(OBJS_LIB)

$(OBJS_UTIL_F77): %.o: %.f
	$(FC) -c $(FFLAGS_INTERNAL) $< -o $@

$(OBJS_MAIN) $(OBJS_SOLVER) $(OBJS_UTIL): %.o: %.f90
	$(FC) -c $(FFLAGS_INTERNAL) $< -o $@

$(OBJS_EIGENEXA): $(SRCS_EIGENEXA)
	$(FC) -c $(FFLAGS_EIGENEXA) $< -o $@

$(OBJS_ELPA): $(SRCS_ELPA)
	$(FC) -c $(FFLAGS_ELPA) $< -o $@

$(OBJS_ELPA_EIGENEXA): $(SRCS_ELPA_EIGENEXA)
	$(FC) -c $(FFLAGS_ELPA_EIGENEXA) $< -o $@

$(TARGET_MAIN): $(OBJS)
	@mkdir -p ./bin
	$(FC) $(FFLAGS_INTERNAL) $(LDFLAGS_INTERNAL) $^ $(LIBS) -o $@
#	cp $(TARGET_MAIN) ./bin/$(TARGET_MAIN)

$(TARGET_LIB): $(OBJS_LIB)
	$(AR) r $(TARGET_LIB) $^

clean:
	@rm -f *.exe bin/* src/*.o *.mod $(TARGET_MAIN) $(TARGET_LIB)

dep:
	find . -name \*.f90 -and ! -name \*_dummy.f90 | xargs makedepf90 -nosrc > src/Makefile.dep

include src/Makefile.dep
