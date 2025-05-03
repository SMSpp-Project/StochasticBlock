##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of test (for StochasticBlock)                                   #
#                                                                            #
#                              Antonio Frangioni                             #
#                          Dipartimento di Informatica                       #
#                              Universita' di Pisa                           #
#                                                                            #
##############################################################################

# module names
NAME = StochasticBlock_test_discretescenarioset

# basic directory
DIR = .

# debug switches
#SW = -g -glldb -fno-inline -std=c++20 -ferror-limit=1
# debug switches with address sanitizer and extra pedantic warning
#SW = -g3 -glldb -fno-inline -std=c++20 -ferror-limit=1 -fsanitize=undefined -fsanitize=address -fno-omit-frame-pointer -Wpedantic -Wextra -Wno-unused-parameter
# production switches
SW = -O3 -std=c++20 -DNDEBUG -DCLANG_1200_0_32_27_PATCH

# debug build switches 
DSW = -g -O0 -std=c++20 -fno-inline

# compiler
CC = clang++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: $(DIR)/$(NAME)

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)/*.o $(DIR)/*~ $(DIR)/$(NAME) temp_*.nc

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
StcBlkSDR = ..
include $(StcBlkSDR)/makefile

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files
MOBJ =  $(SMS++OBJ) $(StcBlkOBJ)

# libraries
MLIB =  $(SMS++LIB) $(StcBlkLIB)

$(DIR)/$(NAME): $(DIR)/test_discretescenarioset.o $(StcBlkOBJ) $(SMS++OBJ)
	$(CC) -o $(DIR)/$(NAME) $^ $(MLIB) $(DSW)

# dependencies: every .o from its .C + every recursively included .h- - - - -

# include directives
MINC =  $(SMS++INC) $(StcBlkINC)

# includes
MH =    $(SMS++H) $(StcBlkH)

# compile command

$(DIR)/test_discretescenarioset.o: $(DIR)/test_discretescenarioset.cpp $(MH)
	$(CC) -c $*.cpp -o $@ $(MINC) -I$(StcBlkSDR)/include $(DSW)

############################ End of makefile #################################
