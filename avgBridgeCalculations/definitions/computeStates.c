#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <omp.h>
#include "../headers/structDefinitions.h"
#include "../headers/inputFunctions.h"
#include "../headers/helperFunctions.h"
#include "../headers/computeBridgesBetweenBins.h"
#include "../headers/computeBridgeYDistribution.h"
#include "../headers/computeBridgeCenterDistribution.h"
#include "../headers/inputParameters.h"
#include "../headers/computeStates.h"
#include "../headers/computeBeadOrientation.h"

STATES initializeStates (STATES currentStates)
{
	currentStates.nBridges = 0;
	currentStates.nLoops = 0;
	currentStates.nDangles = 0;
	currentStates.nFreeChains = 0;

	return currentStates;
}

STATES computeAllStates (STATES currentStates, TRAJECTORY *atoms, int nAtoms)
{
	for (int i = 0; i < nAtoms; )
	{
		if (atoms[i].atomType == 1)
		{
			if (atoms[i].adsorbedID == 0 && atoms[i + 1].adsorbedID == 0) {
				currentStates.nFreeChains++; }

			if ((atoms[i].adsorbedID == 0 && atoms[i + 1].adsorbedID > 0) || (atoms[i].adsorbedID > 0 && atoms[i + 1].adsorbedID == 0)) {
				currentStates.nDangles++; }

			if (atoms[i].adsorbedID == atoms[i + 1].adsorbedID && atoms[i].adsorbedID > 0 && atoms[i + 1].adsorbedID > 0) {
				currentStates.nLoops++; }

			if (atoms[i].adsorbedID > 0 && atoms[i + 1].adsorbedID > 0 && (atoms[i].adsorbedID != atoms[i + 1].adsorbedID)) {
				currentStates.nBridges++; }

			i += 2;
		}
		else if (atoms[i].atomType == 2) {
			i += 1; }
	}

/*	currentStates.nFreeChains /= (float)NPOLYMERS;
	currentStates.nDangles /= (float)NPOLYMERS;
	currentStates.nLoops /= (float)NPOLYMERS;
	currentStates.nBridges /= (float)NPOLYMERS;
*/
/*	printf("%d %d %d %d\n", (int)ceil(currentStates.nFreeChains), (int)ceil(currentStates.nDangles), (int)ceil(currentStates.nLoops), (int)ceil(currentStates.nBridges));
*/
	return currentStates;
}

STATES sumAllStates (STATES currentStates, STATES avgStates)
{

	currentStates.nFreeChains /= (float)NPOLYMERS;
	currentStates.nDangles /= (float)NPOLYMERS;
	currentStates.nLoops /= (float)NPOLYMERS;
	currentStates.nBridges /= (float)NPOLYMERS;

	avgStates.nFreeChains += currentStates.nFreeChains;
	avgStates.nDangles += currentStates.nDangles;
	avgStates.nLoops += currentStates.nLoops;
	avgStates.nBridges += currentStates.nBridges;

	return avgStates;
}

void printCurrentStates (FILE *file_printStates, STATES currentStates)
{
	fprintf(file_printStates, "%d %d %d %d\n", (int)currentStates.nFreeChains, (int)currentStates.nDangles, (int)currentStates.nLoops, (int)currentStates.nBridges);
	fflush (file_printStates);
}

STATES computeAvgStates (STATES avgStates, int nTimeframes)
{
	avgStates.nFreeChains /= nTimeframes;
	avgStates.nDangles /= nTimeframes;
	avgStates.nLoops /= nTimeframes;
	avgStates.nBridges /= nTimeframes;

	return avgStates;
}

STATES *readAllStates (STATES *allStates, int nTimeframes, const char *filename_states)
{
	FILE *file_statesReopen;
	file_statesReopen = fopen (OUTPUT_CURRENT_STATES, "r");

	char lineString[2000];
	int index = 0;

	while (fgets (lineString, 2000, file_statesReopen) != NULL)
	{
		if (lineString[0] != '#')
		{
			sscanf (lineString, "%f %f %f %f\n", 
				&allStates[index].nFreeChains, 
				&allStates[index].nDangles, 
				&allStates[index].nLoops, 
				&allStates[index].nBridges);

			index++;
		}
	}

	fclose (file_statesReopen);

	return allStates;
}

STATES computeStdevStates (STATES stdevStates, STATES avgStates, STATES *allStates, int nTimeframes)
{

	stdevStates.nFreeChains = 0;
	stdevStates.nDangles = 0;
	stdevStates.nLoops = 0;
	stdevStates.nBridges = 0;

	for (int i = 0; i < nTimeframes; ++i)
	{
		stdevStates.nFreeChains += ((allStates[i].nFreeChains - avgStates.nFreeChains) * (allStates[i].nFreeChains - avgStates.nFreeChains));
		stdevStates.nDangles += ((allStates[i].nDangles - avgStates.nDangles) * (allStates[i].nDangles - avgStates.nDangles));
		stdevStates.nLoops += ((allStates[i].nLoops - avgStates.nLoops) * (allStates[i].nLoops - avgStates.nLoops));
		stdevStates.nBridges += ((allStates[i].nBridges - avgStates.nBridges) * (allStates[i].nBridges - avgStates.nBridges));
	}

	stdevStates.nFreeChains /= nTimeframes;
	stdevStates.nDangles /= nTimeframes;
	stdevStates.nLoops /= nTimeframes;
	stdevStates.nBridges /= nTimeframes;

	stdevStates.nFreeChains = sqrt (stdevStates.nFreeChains);
	stdevStates.nDangles = sqrt (stdevStates.nDangles);
	stdevStates.nLoops = sqrt (stdevStates.nLoops);
	stdevStates.nBridges = sqrt (stdevStates.nBridges);

	return stdevStates;
}

void printAverageStates (const char *filename_avgStates, STATES avgStates, STATES stdevStates)
{
	FILE *file_avgStates;
	file_avgStates = fopen (filename_avgStates, "w");

	fprintf(file_avgStates, "# State, average, stdev\n");
	fprintf(file_avgStates, "FreeChains, %f, %f\n", avgStates.nFreeChains, stdevStates.nFreeChains);
	fprintf(file_avgStates, "Dangles, %f, %f\n", avgStates.nDangles, stdevStates.nDangles);
	fprintf(file_avgStates, "Loops, %f, %f\n", avgStates.nLoops, stdevStates.nLoops);
	fprintf(file_avgStates, "Bridges, %f, %f\n", avgStates.nBridges, stdevStates.nBridges);

	fclose (file_avgStates);
}