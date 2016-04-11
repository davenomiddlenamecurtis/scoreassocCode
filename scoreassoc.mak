# Makefile for all programs relating to scoreassoc and pscoreassoc

C = gcc
CC = g++

MAX_LOCI_MAK = 12000
MAX_ALL_MAK = 40
MAX_SUB_MAK = 10000

OURFLAGS = $(CFLAGS) -DMAX_LOCI=$(MAX_LOCI_MAK) -DMAX_ALL=$(MAX_ALL_MAK) -DMAX_SUB=$(MAX_SUB_MAK) 

HEADERS = cdflib.h  dcerror.hpp  dcexpr.hpp  fisher.h  sagcutils.h  safilterfuncs.hpp  scoreassoc.hpp
# cheat and just assume all code dependent on all of these

ifdef INOBJ
all: scoreassoc pscoreassoc pathwayAssoc permPathwayAssoc 
else
all:
	if [ ! -e ../obj ] ; then mkdir ../obj ; fi ; \
	cd ../obj; \
	make -f ../scoreassocCode/scoreassoc.mak INOBJ=INOBJ ; \
	cp scoreassoc pscoreassoc pathwayAssoc permPathwayAssoc ${DCBIN} ; \
	echo copied executables to ${DCBIN} ; \
	cd ../scoreassocCode
endif
# unless you have a folder on path called $DCBIN you should delete that line to copy the executables to it

clean:
	rm ../obj/*.o

VPATH=../scoreassocCode
	
%.o: ../scoreassocCode/%.cpp $(HEADERS)
	$(CC) $(OURFLAGS) -c $< -o ../obj/$@
	
%.o: ../scoreassocCode/%.c $(HEADERS)
	$(C) $(OURFLAGS) -c $< -o ../obj/$@

scoreassoc: scoreassoc.o saglobals.o scoreassocfuncs.o sarecfuncs.o sahaprecfuncs.o satriofuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o saFilterFuncs.o 
	$(CC) -o scoreassoc scoreassoc.o saglobals.o scoreassocfuncs.o sarecfuncs.o sahaprecfuncs.o satriofuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o saFilterFuncs.o -lm

pscoreassoc: pscoreassoc.o saglobals.o scoreassocfuncs.o sarecfuncs.o sahaprecfuncs.o satriofuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o saFilterFuncs.o 
	$(CC) -o pscoreassoc pscoreassoc.o saglobals.o scoreassocfuncs.o sarecfuncs.o sahaprecfuncs.o satriofuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o saFilterFuncs.o -lm

pathwayAssoc: pathwayAssoc.o dcdflib.o ipmpar.o dcerror.o 
	$(CC) -o pathwayAssoc pathwayAssoc.o dcdflib.o ipmpar.o dcerror.o -lm

permPathwayAssoc: permPathwayAssoc.o dcdflib.o ipmpar.o dcerror.o 
	$(CC) -o permPathwayAssoc permPathwayAssoc.o dcdflib.o ipmpar.o dcerror.o -lm

testMRVSpread: testMRVSpread.o scoreassocfuncs.o sarecfuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o saFilterFuncs.o 
	$(CC) -o testMRVSpread testMRVSpread.o scoreassocfuncs.o sarecfuncs.o sagcutils.o dcdflib.o ipmpar.o dcerror.o dcexpr.o saFilterFuncs.o -lm

