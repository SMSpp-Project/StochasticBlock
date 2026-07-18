##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of StochasticBlock                                              #
#                                                                            #
#   The makefile takes in input the -I directives for all the external       #
#   libraries needed by StochasticBlock, i.e., core SMS++.                   #
#                                                                            #
#   Note that, conversely, $(SMS++INC) is also assumed to include any        #
#   -I directive corresponding to external libraries needed by SMS++, at     #
#   least to the extent in which they are needed by the parts of SMS++       #
#   used by StochasticBlock.                                                 #
#                                                                            #
#   Input:  $(CC)          = compiler command                                #
#           $(SW)          = compiler options                                #
#           $(SMS++INC)    = the -I$( core SMS++ directory )                 #
#           $(SMS++OBJ)    = the core SMS++ library                          #
#           $(StcBlkSDR)   = the directory where the source is               #
#                                                                            #
#   Output: $(StcBlkOBJ)   = the final object(s) / library                   #
#           $(StcBlkH)     = the .h files to include                         #
#           $(StcBlkINC)   = the -I$( source directory )                     #
#                                                                            #
#                              Antonio Frangioni                             #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

StcBlkOBJ = $(StcBlkSDR)/obj/DiscreteScenarioSet.o \
	$(StcBlkSDR)/obj/IndependentMultiStageScenarioGenerator.o \
	$(StcBlkSDR)/obj/MultiStageDiscreteScenarioSet.o \
	$(StcBlkSDR)/obj/ScenarioReductionBlock.o \
	$(StcBlkSDR)/obj/StochasticBlock.o

StcBlkLIB = 

StcBlkINC = -I$(StcBlkSDR)/include

StcBlkH   = $(StcBlkSDR)/include/ScenarioGenerator.h \
	$(StcBlkSDR)/include/DiscreteScenarioSet.h \
	$(StcBlkSDR)/include/IndependentMultiStageScenarioGenerator.h \
	$(StcBlkSDR)/include/MultiStageDiscreteScenarioSet.h \
	$(StcBlkSDR)/include/ScenarioReductionBlock.h \
	$(StcBlkSDR)/include/StochasticBlock.h

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(StcBlkOBJ) $(StcBlkSDR)/*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(StcBlkSDR)/obj/StochasticBlock.o: $(StcBlkSDR)/src/StochasticBlock.cpp \
	$(StcBlkSDR)/include/StochasticBlock.h $(SMS++OBJ)
	$(CC) -c $(StcBlkSDR)/src/StochasticBlock.cpp -o $@ $(StcBlkINC) \
	$(SMS++INC) $(SW)

# Note: DiscreteScenarioSet requires CapacitatedFacilityLocationBlock for
# scenario reduction functionality. The $(CFLBkINC) dependency must be
# provided by the including makefile when DiscreteScenarioSet is used.
$(StcBlkSDR)/obj/DiscreteScenarioSet.o: \
	$(StcBlkSDR)/src/DiscreteScenarioSet.cpp \
	$(StcBlkSDR)/include/ScenarioGenerator.h \
	$(StcBlkSDR)/include/DiscreteScenarioSet.h $(SMS++OBJ)
	$(CC) -c $(StcBlkSDR)/src/DiscreteScenarioSet.cpp -o $@ \
	$(StcBlkINC) $(CFLBkINC) $(SMS++INC) $(SW)

# IndependentMultiStageScenarioGenerator only needs the (abstract)
# ScenarioGenerator and core SMS++ (netCDF, factory).
$(StcBlkSDR)/obj/IndependentMultiStageScenarioGenerator.o: \
	$(StcBlkSDR)/src/IndependentMultiStageScenarioGenerator.cpp \
	$(StcBlkSDR)/include/ScenarioGenerator.h \
	$(StcBlkSDR)/include/IndependentMultiStageScenarioGenerator.h $(SMS++OBJ)
	$(CC) -c $(StcBlkSDR)/src/IndependentMultiStageScenarioGenerator.cpp \
	-o $@ $(StcBlkINC) $(SMS++INC) $(SW)

# MultiStageDiscreteScenarioSet only needs the (abstract) ScenarioGenerator
# and core SMS++ (netCDF, factory); no CapacitatedFacilityLocationBlock.
$(StcBlkSDR)/obj/MultiStageDiscreteScenarioSet.o: \
	$(StcBlkSDR)/src/MultiStageDiscreteScenarioSet.cpp \
	$(StcBlkSDR)/include/ScenarioGenerator.h \
	$(StcBlkSDR)/include/MultiStageDiscreteScenarioSet.h $(SMS++OBJ)
	$(CC) -c $(StcBlkSDR)/src/MultiStageDiscreteScenarioSet.cpp -o $@ \
	$(StcBlkINC) $(SMS++INC) $(SW)

$(StcBlkSDR)/obj/ScenarioReductionBlock.o: \
	$(StcBlkSDR)/src/ScenarioReductionBlock.cpp \
	$(StcBlkSDR)/include/ScenarioReductionBlock.h \
	$(StcBlkSDR)/include/ScenarioGenerator.h $(SMS++OBJ)
	$(CC) -c $(StcBlkSDR)/src/ScenarioReductionBlock.cpp -o $@ \
	$(StcBlkINC) $(SMS++INC) $(SW)

########################## End of makefile ###################################
