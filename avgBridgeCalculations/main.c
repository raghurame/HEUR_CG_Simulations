#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <omp.h>
#include "headers/structDefinitions.h"
#include "headers/inputFunctions.h"
#include "headers/helperFunctions.h"
#include "headers/computeBridgesBetweenBins.h"
#include "headers/computeBridgeYDistribution.h"
#include "headers/computeBridgeCenterDistribution.h"
#include "headers/inputParameters.h"
#include "headers/computeStates.h"
#include "headers/computeBeadOrientation.h"

/*
LIST OF FUNCTIONS IN HEADER FILES:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. structDefinitions.h - 
		BOUNDARY, 
		TRAJECTORY, 
		BONDINFO, 
		YDIST, 
		BRIDGESBIN

2. inputFunctions.h - 
		BOUNDARY getBoundary (FILE *file_inputTrj, BOUNDARY simBoundary), 
		int countNAtoms (FILE *file_inputTrj), 
		TRAJECTORY *getAtoms (TRAJECTORY *atoms, int nAtoms, BOUNDARY simBoundary, 
		float distanceCutoff, FILE *file_inputTrj, 
		int file_status, TRAJECTORY **micelles, int nMicelles), 
		int countNMicelles (FILE *file_inputTrj, int nAtoms).

3. helperFunctions.h
		float translatePeriodic (float r1, float r2, float simulationBoxLength),
		float computePeriodicDistance (float x1, float y1, float z1, float x2, float y2, float z2, float xLength, float yLength, float zLength).

4. computeBridgesBetweenBins.h
		BRIDGESBIN *assignBinBounds (BRIDGESBIN *bridgeBetweenBins, BOUNDARY simBoundary, float binWidth, float delBinDistance, int nBins),
		BRIDGESBIN *countBridgesBetweenBins (TRAJECTORY **atoms, BOUNDARY simBoundary, float distanceCutoff, BRIDGESBIN *bridgeBetweenBins, int nAtoms, TRAJECTORY *micelles, int nMicelles, int nBins).

5. computeBridgeYDistribution.h
		YDIST *assignBridgeYDistribution (float maxFeneExtension, int nBins, float binWidth, YDIST *bridgeYDistribution),
		YDIST *computeBridgeDistribution (TRAJECTORY *atoms, int nAtoms, YDIST *bridgeYDistribution, int nBins, BOUNDARY simBoundary).

6. computeBridgeCenterDistribution.h
		BRIDGESBIN *assignBridgeCenterDistribution (BRIDGESBIN *bridgeCenterDistribution, int nBins_centerDistribution, float binWidth_centerDistribution, BOUNDARY simBoundary),
		BONDINFO *computeBridgeCenter (TRAJECTORY *atoms, int nAtoms, BONDINFO *allBonds, BOUNDARY simBoundary),
		BRIDGESBIN *computeBridgeCenterDistribution (BONDINFO *allBonds, int nBonds, BRIDGESBIN *bridgeCenterDistribution, int nBins_centerDistribution).
*/

int main(int argc, char const *argv[])
{
	system ("mkdir outputs");

	FILE *file_inputTrj;
	file_inputTrj = fopen (INPUTFILE, "r");

	int nAtoms = countNAtoms (file_inputTrj), file_status = 0, nMicelles = countNMicelles (file_inputTrj, nAtoms), nTimeframes = 0;
	BOUNDARY simBoundary;
	simBoundary = getBoundary (file_inputTrj, simBoundary);

	char lineString[2000];
	float binWidth_vertBridges = BINWIDTH_VERTICALBRIDGES, distanceCutoff_vertBridges = ADSORPTION_DISTANCE_CUTOFF, delBinDistance_vertBridges = DELBINDISTANCE_VERTICALBRIDGES;
	int nBins_vertBridges = ceil ((((simBoundary.yhi - simBoundary.ylo) * 2) - (2 * binWidth_vertBridges)) / delBinDistance_vertBridges);

	BRIDGESBIN *bridgeBetweenBins;
	bridgeBetweenBins = (BRIDGESBIN *) malloc (nBins_vertBridges * sizeof (BRIDGESBIN));

	TRAJECTORY *atoms, *micelles;
	atoms = (TRAJECTORY *) malloc (nAtoms * sizeof (TRAJECTORY));
	micelles = (TRAJECTORY *) malloc (nMicelles * sizeof (TRAJECTORY));

	int nBins_yDist = NBINS_YDISTRIBUTIONBRIDGES;
	float maxFeneExtension = MAX_FENE_EXTENSION, binWidth_yDist = (maxFeneExtension / (float)nBins_yDist);
	YDIST *bridgeYDistribution;
	bridgeYDistribution = (YDIST *) malloc (nBins_yDist * sizeof (YDIST));

	int nBonds = (nAtoms - nMicelles) / 2, nBins_bridgeCenter = NBINS_YDISTRIBUTIONBRIDGES;
	BONDINFO *allBonds;
	allBonds = (BONDINFO *) malloc (nBonds * sizeof (BONDINFO));
	BRIDGESBIN *bridgeCenterDistribution;

	float binWidth_centerDistribution = BINWIDTH_BONDCENTERDISTRIBUTION;
	int nBins_centerDistribution = ceil (simBoundary.yLength / binWidth_centerDistribution);
	bridgeCenterDistribution = (BRIDGESBIN *) malloc (nBins_centerDistribution * sizeof (BRIDGESBIN));

	bridgeBetweenBins = assignBinBounds (bridgeBetweenBins, simBoundary, binWidth_vertBridges, delBinDistance_vertBridges, nBins_vertBridges);
	bridgeYDistribution = assignBridgeYDistribution (maxFeneExtension, nBins_yDist, binWidth_yDist, bridgeYDistribution);
	bridgeCenterDistribution = assignBridgeCenterDistribution (bridgeCenterDistribution, nBins_centerDistribution, binWidth_centerDistribution, simBoundary);

	STATES currentStates, avgStates, stdevStates, *allStates;
	avgStates = initializeStates (avgStates);
	stdevStates = initializeStates (stdevStates);

	FILE *file_printStates;
	file_printStates = fopen (OUTPUT_CURRENT_STATES, "w");
	fprintf(file_printStates, "# free, dangles, loop, bridges\n");

	ANGLE_DISTRIBUTION *bridgeOrientation, *freeOrientation, *dangleOrientation, *loopOrientation;
	int nBins_orientation = ceil (360.0 / BINWIDTH_ANGLEORIENTATION);

	bridgeOrientation = (ANGLE_DISTRIBUTION *) malloc (nBins_orientation * sizeof (ANGLE_DISTRIBUTION));
	freeOrientation = (ANGLE_DISTRIBUTION *) malloc (nBins_orientation * sizeof (ANGLE_DISTRIBUTION));
	dangleOrientation = (ANGLE_DISTRIBUTION *) malloc (nBins_orientation * sizeof (ANGLE_DISTRIBUTION));
	loopOrientation = (ANGLE_DISTRIBUTION *) malloc (nBins_orientation * sizeof (ANGLE_DISTRIBUTION));

	bridgeOrientation = assignBeadOrientation (bridgeOrientation, nBins_orientation);
	freeOrientation = assignBeadOrientation (freeOrientation, nBins_orientation);
	dangleOrientation = assignBeadOrientation (dangleOrientation, nBins_orientation);
	loopOrientation = assignBeadOrientation (loopOrientation, nBins_orientation);

	while (file_status != EOF)
	{
		if (MAXTIMESTEPS != 0 && nTimeframes >= MAXTIMESTEPS) {
			break; }

		if (nTimeframes%NFREQ_PROGRESS_SCREEN == 0) {
			fprintf(stdout, "computing %d timesteps...         \r", nTimeframes);
			fflush (stdout); }

		atoms = getAtoms (atoms, nAtoms, simBoundary, distanceCutoff_vertBridges, file_inputTrj, file_status, &micelles, nMicelles);

		bridgeBetweenBins = countBridgesBetweenBins (&atoms, simBoundary, distanceCutoff_vertBridges, bridgeBetweenBins, nAtoms, micelles, nMicelles, nBins_vertBridges);
		bridgeYDistribution = computeBridgeDistribution (atoms, nAtoms, bridgeYDistribution, nBins_yDist, simBoundary);

		// Finish the bridge center distribution calculations
		allBonds = computeBridgeCenter (atoms, nAtoms, allBonds, simBoundary);
		bridgeCenterDistribution = computeBridgeCenterDistribution (allBonds, nBonds, bridgeCenterDistribution, nBins_centerDistribution);

		// Calculate the overall bridges, loops, dangles, free chains
		currentStates = initializeStates (currentStates);
		currentStates = computeAllStates (currentStates, atoms, nAtoms);
		avgStates = sumAllStates (currentStates, avgStates);
		printCurrentStates (file_printStates, currentStates);

		// Calculate the distribution of centers of bridges, loops, dangles, free chains

		// Calculate the distribution of orientation angles of bridges, loops, dangles, free chains
		computeBeadOrientationDistribution (atoms, nBins_orientation, allBonds, nBonds, &bridgeOrientation, &loopOrientation, &dangleOrientation, &loopOrientation);

		file_status = fgetc (file_inputTrj);
		nTimeframes++;
	}

	fclose (file_printStates);

	avgStates = computeAvgStates (avgStates, nTimeframes);

	allStates = (STATES *) malloc (nTimeframes * sizeof (STATES));
	allStates = readAllStates (allStates, nTimeframes, OUTPUT_CURRENT_STATES);

	stdevStates = initializeStates (stdevStates);
	stdevStates = computeStdevStates (stdevStates, avgStates, allStates, nTimeframes);

	printAverageStates (OUTPUT_AVG_STATES, avgStates, stdevStates);

	FILE *file_bridgeBetweenBinsOuptut, *file_bridgeYDistributionOutput, *file_bridgeCenterDistributionOutput;
	file_bridgeBetweenBinsOuptut = fopen (OUTPUT_BRIDGESBETWEENBINS, "w");
	file_bridgeYDistributionOutput = fopen (OUTPUT_YDISTRIBUTIONBRIDGES, "w");
	file_bridgeCenterDistributionOutput = fopen (OUTPUT_BONDCENTERDISTRIBUTION, "w");

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
	fclose (file_bridgeBetweenBinsOuptut);
	fclose (file_bridgeYDistributionOutput);
	fclose (file_bridgeCenterDistributionOutput);

	return 0;
}