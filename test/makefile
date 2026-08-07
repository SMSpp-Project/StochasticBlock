##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of the StochasticBlock tests                                    #
#                                                                            #
#                              Antonio Frangioni                             #
#                          Dipartimento di Informatica                       #
#                              Universita' di Pisa                           #
#                                                                            #
##############################################################################

# the executables built here
NAMES = StochasticBlock_test test_discrete test_mssds test_imssg

# basic directory
DIR = .

# common flags
COMMON_SW = -std=c++20 -ferror-limit=1

# debug switches
SW_DEBUG = -g -glldb -fno-inline $(COMMON_SW)
# debug switches with address sanitizer and extra pedantic warning
SW_DEBUG_ASAN = -g3 -glldb -fno-inline $(COMMON_SW) -fsanitize=undefined -fsanitize=address -fno-omit-frame-pointer -Wpedantic -Wextra -Wno-unused-parameter
# production switches
SW_RELEASE = -O3 $(COMMON_SW) -DNDEBUG -DCLANG_1200_0_32_27_PATCH

# compiler
CC = clang++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: release

# debug target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

debug: SW = $(SW_DEBUG)
debug: $(addprefix $(DIR)/,$(NAMES))

# debug_asan target - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

debug_asan: SW = $(SW_DEBUG_ASAN)
debug_asan: $(addprefix $(DIR)/,$(NAMES))

# release target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

release: SW = $(SW_RELEASE)
release: $(addprefix $(DIR)/,$(NAMES))

# clean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)/*.o $(DIR)/*~ $(addprefix $(DIR)/,$(NAMES))

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# each module outputs some macros to be used here:
# *OBJ is the final object(s) / library
# *LIB is the external libraries + -L< libdirs >
# *H   is the list of all include files
# *INC is the -I< include directories >

# StochasticBlock complete + core SMS++
StcBlkSDR = ..
include $(StcBlkSDR)/makefile-c

# MILPSolver (without core SMS++), only needed by test_discrete
MILPSSDR = ../../MILPSolver
include $(MILPSSDR)/makefile-s

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# StochasticBlock_test needs only StochasticBlock + core
$(DIR)/StochasticBlock_test: $(StcBlkOBJ) $(DIR)/test.o
	$(CC) -o $@ $(DIR)/test.o $(StcBlkOBJ) $(StcBlkLIB) $(SW)

# test_discrete additionally needs MILPSolver
$(DIR)/test_discrete: $(StcBlkOBJ) $(MILPSOBJ) $(DIR)/test_discrete.o
	$(CC) -o $@ $(DIR)/test_discrete.o $(StcBlkOBJ) $(MILPSOBJ) \
	      $(StcBlkLIB) $(MILPSLIB) $(SW)

# test_mssds needs only StochasticBlock + core
$(DIR)/test_mssds: $(StcBlkOBJ) $(DIR)/test_mssds.o
	$(CC) -o $@ $(DIR)/test_mssds.o $(StcBlkOBJ) $(StcBlkLIB) $(SW)

# test_imssg needs only StochasticBlock + core
$(DIR)/test_imssg: $(StcBlkOBJ) $(DIR)/test_imssg.o
	$(CC) -o $@ $(DIR)/test_imssg.o $(StcBlkOBJ) $(StcBlkLIB) $(SW)

# dependencies: every .o from its .cpp + every recursively included .h- - - -

# include directives
MINC =  $(StcBlkINC) $(MILPSINC)

# includes
MH =    $(StcBlkH) $(MILPSH)

# compile command

$(DIR)/%.o: $(DIR)/%.cpp $(MH)
	$(CC) -c $< -o $@ $(MINC) $(SW)

# distclean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

distclean: clean

# phony targets - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.PHONY: debug debug_asan release clean distclean

############################ End of makefile #################################
