#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <omp.h>
#include "structDefinitions.h"
#include "inputFunctions.h"
#include "helperFunctions.h"
#include "computeBridgesBetweenBins.h"
#include "computeBridgeYDistribution.h"
#include "computeBridgeCenterDistribution.h"

int main(int argc, char const *argv[])
{
	if (argc != 5)
	{
		printf("\nINCORRECT ARGUMENTS PASSED:\n~~~~~~~~~~~~~~~~~~~~~~~~~\n\n -> Pass the input file as (char *) argv[1]\n -> Enter bin width as (float) argv[2]\n -> Enter the distance cut off for particle-bead as (float) argv[3]\n -> Enter the increment in bin distance as (float) argv[4]\n\nFor example: ./countBridges last.lammpstrj 5.0 3.38 1.0\n~~~~~~~~~~~\n\nTo compile: gcc -o countBridges countBridges.c -lm\n~~~~~~~~~~\n\n");
		exit (1);
	}

	FILE *file_inputTrj;
	file_inputTrj = fopen (argv[1], "r");

	int nAtoms = countNAtoms (file_inputTrj), file_status = 0, nMicelles = countNMicelles (file_inputTrj, nAtoms), nTimeframes = 0;
	BOUNDARY simBoundary;
	simBoundary = getBoundary (file_inputTrj, simBoundary);

	char lineString[2000];
	float binWidth_vertBridges = atof (argv[2]), distanceCutoff_vertBridges = atof (argv[3]), delBinDistance_vertBridges = atof (argv[4]);
	int nBins_vertBridges = ceil ((((simBoundary.yhi - simBoundary.ylo) * 2) - (2 * binWidth_vertBridges)) / delBinDistance_vertBridges);

	BRIDGESBIN *bridgeBetweenBins;
	bridgeBetweenBins = (BRIDGESBIN *) malloc (nBins_vertBridges * sizeof (BRIDGESBIN));

	TRAJECTORY *atoms, *micelles;
	atoms = (TRAJECTORY *) malloc (nAtoms * sizeof (TRAJECTORY));
	micelles = (TRAJECTORY *) malloc (nMicelles * sizeof (TRAJECTORY));

	float maxFeneExtension = 60.0, binWidth_yDist = (maxFeneExtension / 20);
	int nBins_yDist = 20; //Taken arbitrarily for now; 20 bins across the maximum extensible length of 60 sigma
	YDIST *bridgeYDistribution;
	bridgeYDistribution = (YDIST *) malloc (nBins_yDist * sizeof (YDIST));

	int nBonds = (nAtoms - nMicelles) / 2, nBins_bridgeCenter = 20;
	BONDINFO *allBonds;
	allBonds = (BONDINFO *) malloc (nBonds * sizeof (BONDINFO));
	BRIDGESBIN *bridgeCenterDistribution;

	float binWidth_centerDistribution = 3;
	int nBins_centerDistribution = ceil (simBoundary.yLength / binWidth_centerDistribution);
	bridgeCenterDistribution = (BRIDGESBIN *) malloc (nBins_centerDistribution * sizeof (BRIDGESBIN));

	bridgeBetweenBins = assignBinBounds (bridgeBetweenBins, simBoundary, binWidth_vertBridges, delBinDistance_vertBridges, nBins_vertBridges);
	bridgeYDistribution = assignBridgeYDistribution (maxFeneExtension, nBins_yDist, binWidth_yDist, bridgeYDistribution);
	bridgeCenterDistribution = assignBridgeCenterDistribution (bridgeCenterDistribution, nBins_centerDistribution, binWidth_centerDistribution, simBoundary);

	while (file_status != EOF)
	{
		if (nTimeframes%1 == 0) {
			fprintf(stdout, "computing %d timesteps...         \r", nTimeframes);
			fflush (stdout); }

		atoms = getAtoms (atoms, nAtoms, simBoundary, distanceCutoff_vertBridges, file_inputTrj, file_status, &micelles, nMicelles);

		bridgeBetweenBins = countBridgesBetweenBins (&atoms, simBoundary, distanceCutoff_vertBridges, bridgeBetweenBins, nAtoms, micelles, nMicelles, nBins_vertBridges);
		bridgeYDistribution = computeBridgeDistribution (atoms, nAtoms, bridgeYDistribution, nBins_yDist, simBoundary);

		// Finish the bridge center distribution calculations
		allBonds = computeBridgeCenter (atoms, nAtoms, allBonds, simBoundary);
		bridgeCenterDistribution = computeBridgeCenterDistribution (allBonds, nBonds, bridgeCenterDistribution, nBins_centerDistribution);

		file_status = fgetc (file_inputTrj);
		nTimeframes++;
	}

	FILE *file_bridgeBetweenBinsOuptut, *file_bridgeYDistributionOutput, *file_bridgeCenterDistributionOutput;
	file_bridgeBetweenBinsOuptut = fopen ("nBridgesBetweenBins.count", "w");
	file_bridgeYDistributionOutput = fopen ("bridges.distribution", "w");
	file_bridgeCenterDistributionOutput = fopen ("bridgeCenter.distribution", "w");

	fprintf(file_bridgeBetweenBinsOuptut, "# y1lo, y1hi, y2lo, y2hi, avgCounts\n");
	fprintf(file_bridgeYDistributionOutput, "# ylo, yhi, avgCounts\n");
	fprintf(file_bridgeCenterDistributionOutput, "# ylo, yhi, avgCounts\n");

	for (int i = 0; i < nBins_vertBridges; ++i) {
		fprintf(file_bridgeBetweenBinsOuptut, "%f %f %f %f %f\n", 
			bridgeBetweenBins[i].y1lo, 
			bridgeBetweenBins[i].y1hi, 
			bridgeBetweenBins[i].y2lo, 
			bridgeBetweenBins[i].y2hi, 
			((float)bridgeBetweenBins[i].count / (float)nTimeframes)); }

	for (int i = 0; i < nBins_yDist; ++i) {
		fprintf(file_bridgeYDistributionOutput, "%f %f %f\n", 
			bridgeYDistribution[i].ylo, 
			bridgeYDistribution[i].yhi, 
			((float)bridgeYDistribution[i].count / (float)nTimeframes)); }

	for (int i = 0; i < nBins_centerDistribution; ++i) {
		fprintf(file_bridgeCenterDistributionOutput, "%f %f %f\n", 
			bridgeCenterDistribution[i].y1lo, 
			bridgeCenterDistribution[i].y1hi, 
			(float)bridgeCenterDistribution[i].count / (float)nTimeframes); }

	fclose (file_inputTrj);
	return 0;
}