##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of test for DiscreteScenarioSet                                 #
#                                                                            #
#                              Antonio Frangioni                             #
#                          Dipartimento di Informatica                       #
#                              Universita' di Pisa                           #
#                                                                            #
##############################################################################

# module name - must match the binary name and source file
NAME = test_discretescenarioset

# basic directory
DIR = .

# common flags
COMMON_SW = -std=c++20 -ferror-limit=1 -Wno-enum-compare

# debug switches
SW_DEBUG = -g3 -glldb -fno-inline $(COMMON_SW) -DCLANG_1200_0_32_27_PATCH
# debug switches with address sanitizer and extra pedantic warning
SW_DEBUG_ASAN = -g3 -glldb -fno-inline $(COMMON_SW) -fsanitize=undefined -fsanitize=address -fno-omit-frame-pointer -Wpedantic -Wextra -Wno-unused-parameter -Wno-enum-compare -Wno-deprecated-copy-with-user-provided-copy
# production switches
SW_RELEASE = -O3 $(COMMON_SW) -DNDEBUG -DCLANG_1200_0_32_27_PATCH

# compiler
CC = clang++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: release

# debug target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

debug: SW = $(SW_DEBUG)
debug: $(DIR)/$(NAME)

# debug_asan target - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

debug_asan: SW = $(SW_DEBUG_ASAN)
debug_asan: $(DIR)/$(NAME)

# release target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

release: SW = $(SW_RELEASE)
release: $(DIR)/$(NAME)

# clean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)/*.o $(DIR)/*~ $(NAME)

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object(s) / library
# *LIB is the external libraries + -L< libdirs >
# *H   is the list of all include files
# *INC is the -I< include directories >

# define input macros for SMS++ complete makefile, then include it
SMS++SDR = ../../SMS++
include $(SMS++SDR)/lib/makefile-inc

# StochasticBlock
StcBlkSDR = ../../StochasticBlock
include $(StcBlkSDR)/makefile

# CapacitatedFacilityLocationBlock
CFLBkSDR = ../../CapacitatedFacilityLocationBlock
include $(CFLBkSDR)/makefile-s

# BinaryKnapsackBlock (needed by CapacitatedFacilityLocationBlock)
BKPBkSDR = ../../BinaryKnapsackBlock
include $(BKPBkSDR)/makefile

# MCFBlock (needed by CapacitatedFacilityLocationBlock)  
MCFBkSDR = ../../MCFBlock
include $(MCFBkSDR)/makefile

# BundleSolver
BNDSLVSDR = ../../BundleSolver
include $(BNDSLVSDR)/makefile-s

# LagrangianDualSolver
LgDSLVSDR = ../../LagrangianDualSolver
include $(LgDSLVSDR)/makefile

# MILPSolver - we need it for ScenarioReductionSolver
MILPSSDR = ../../MILPSolver
include $(MILPSSDR)/makefile-c

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files - includes all necessary objects
MOBJ = $(SMS++OBJ) $(StcBlkOBJ) $(BNDSLVOBJ) $(LgDSLVOBJ) $(CFLBkOBJ) $(BKPBkOBJ) $(MCFLEOBJ)

# libraries - includes all necessary libraries
MLIB =  $(BNDSLVLIB) $(LgDSLVLIB) $(CFLBkLIB) $(BKPBkLIB) $(SMS++LIB) $(MCFLELIB) $(MILPSLIB) $(StcBlkLIB)

# Final compilation target
$(DIR)/$(NAME): $(MOBJ) $(DIR)/$(NAME).o
	$(CC) -o $(DIR)/$(NAME) $^ $(MLIB) $(SW)

# dependencies: every .o from its .cpp + every recursively included .h- - - -

# include directives - add MILPS_INC as an extra safety measure
MINC =  $(SMS++INC) $(StcBlkINC) $(BNDSLVINC) $(LgDSLVINC) $(CFLBkINC) $(BKPBkINC) $(MCFBkINC) $(MCFLEINC) $(MILPSINC) $(MILPS_INC) -I../include -I../build

# includes
MH =    $(SMS++H) $(StcBlkH) $(BNDSLVH) $(LgDSLVH) $(CFLBkH) $(BKPBkH) $(MCFBkH) $(MCFLEH) $(MILPSH)

# compile command
$(DIR)/$(NAME).o: $(DIR)/$(NAME).cpp $(MH)
	$(CC) -c $(DIR)/$(NAME).cpp -o $@ $(MINC) $(SW)

# distclean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

distclean: clean

# phony targets - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.PHONY: debug debug_asan release clean distclean

############################ End of makefile #################################