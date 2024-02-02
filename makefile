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

StcBlkOBJ = $(StcBlkSDR)/obj/StochasticBlock.o 

StcBlkINC = -I$(StcBlkSDR)/include

StcBlkH   = $(StcBlkSDR)/include/StochasticBlock.h 

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(StcBlkOBJ) $(StcBlkSDR)/*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(StcBlkOBJ): $(StcBlkSDR)/src/StochasticBlock.cpp $(StcBlkH) $(SMS++OBJ)
	$(CC) -c $(StcBlkSDR)/src/StochasticBlock.cpp -o $@ $(StcBlkINC) \
	$(SMS++INC) $(SW)

########################## End of makefile ###################################
