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

YDIST *assignBridgeYDistribution (float maxFeneExtension, int nBins, float binWidth, YDIST *bridgeYDistribution)
{
	omp_set_num_threads (NTHREADS);

	#pragma omp parallel for
	for (int i = 0; i < nBins; ++i)
	{
		if (i == 0) {
			bridgeYDistribution[i].ylo = 0; }
		else {
			bridgeYDistribution[i].ylo = bridgeYDistribution[i - 1].yhi; }

		bridgeYDistribution[i].yhi = bridgeYDistribution[i].ylo + binWidth;
		bridgeYDistribution[i].count = 0;
	}

	return bridgeYDistribution;
}

YDIST *computeBridgeDistribution (TRAJECTORY *atoms, int nAtoms, YDIST *bridgeYDistribution, int nBins, BOUNDARY simBoundary)
{
	omp_set_num_threads (NTHREADS);

	int bridgeCountLocal = 0;
	float yDistance;
	float tempY;

	for (int i = 0; i < nAtoms; )
	{
		if ((atoms[i].adsorbedID != atoms[i + 1].adsorbedID) && atoms[i].adsorbedID != 0 && atoms[i + 1].adsorbedID != 0)
		{
			tempY = translatePeriodic (atoms[i].y, atoms[i + 1].y, simBoundary.yLength);

/*			if (tempY > y2)
			{
				tempY += simBoundary.xy;
			}
			else if (tempY < y2)
			{
				tempY -= simBoundary.xy;
			}
*/
			yDistance = fabs (atoms[i].y - tempY);

			#pragma omp parallel for
			for (int j = 0; j < nBins; ++j)
			{
				if (yDistance > bridgeYDistribution[j].ylo && yDistance <= bridgeYDistribution[j].yhi) {
					bridgeYDistribution[j].count++; }
			}
		}

		if ((atoms[i].adsorbedID != atoms[i + 1].adsorbedID) && atoms[i].atomType == 1 && atoms[i + 1].atomType == 1 && atoms[i].adsorbedID > 0 && atoms[i + 1].adsorbedID > 0) {
			bridgeCountLocal++; }

		if (atoms[i].atomType == 1) {
			i += 2; }
		else if (atoms[i].atomType == 2) {
			i += 1; }
	}

	return bridgeYDistribution;
}