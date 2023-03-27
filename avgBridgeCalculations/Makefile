all:
	#
	# MAKEFILE INSTRUCTIONS:
	# ~~~~~~~~~~~~~~~~~~~~~~
	#
	# {~} 'make compile' --> compile the code to generate executable
	# {~} 'make run' --> run the executable
	#
compile: countBridges.c inputFunctions.c helperFunctions.c computeBridgesBetweenBins.c computeBridgeYDistribution.c computeBridgeCenterDistribution.c structDefinitions.h
	export OMP_DYNAMIC=true
	gcc -c main.c -Wall -fstack-protector -g -fopenmp -lm
	gcc -c inputFunctions.c -Wall -fstack-protector -g -fopenmp -lm
	gcc -c helperFunctions.c -Wall -fstack-protector -g -fopenmp -lm
	gcc -c computeBridgesBetweenBins.c -Wall -fstack-protector -g -fopenmp -lm
	gcc -c computeBridgeYDistribution.c -Wall -fstack-protector -g -fopenmp -lm
	gcc -c computeBridgeCenterDistribution.c -Wall -fstack-protector -g -fopenmp -lm
	gcc main.o inputFunctions.o helperFunctions.o computeBridgesBetweenBins.o computeBridgeYDistribution.o computeBridgeCenterDistribution.o -o countBridges -lm -fopenmp -Wall -g -fstack-protector
run: countBridges
	./countBridges
