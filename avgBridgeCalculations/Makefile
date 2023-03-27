all:
	#
	# MAKEFILE COMMANDS:
	# ~~~~~~~~~~~~~~~~~
	#
	# {~} 'make compile' --> compile the code to generate executable
	# {~} 'make run'     --> run the executable
	# {~} 'make clean'   --> remove all *.o files
	#
	# NOTE: Input variables can be edited from headers/inputParameters.h
	#
compile: main.c definitions/inputFunctions.c definitions/helperFunctions.c definitions/computeBridgesBetweenBins.c definitions/computeBridgeYDistribution.c definitions/computeBridgeCenterDistribution.c headers/structDefinitions.h
	export OMP_DYNAMIC=true
	gcc -c main.c -Wall -fstack-protector -g -fopenmp -lm
	gcc -c definitions/inputFunctions.c -Wall -fstack-protector -g -fopenmp -lm
	gcc -c definitions/helperFunctions.c -Wall -fstack-protector -g -fopenmp -lm
	gcc -c definitions/computeBridgesBetweenBins.c -Wall -fstack-protector -g -fopenmp -lm
	gcc -c definitions/computeBridgeYDistribution.c -Wall -fstack-protector -g -fopenmp -lm
	gcc -c definitions/computeBridgeCenterDistribution.c -Wall -fstack-protector -g -fopenmp -lm
	gcc main.o inputFunctions.o helperFunctions.o computeBridgesBetweenBins.o computeBridgeYDistribution.o computeBridgeCenterDistribution.o -o countBridges -lm -fopenmp -Wall -g -fstack-protector
run: countBridges
	./countBridges
clean:
	rm -rf *.o