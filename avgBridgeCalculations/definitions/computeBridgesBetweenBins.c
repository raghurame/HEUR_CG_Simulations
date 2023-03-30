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

BRIDGESBIN *assignBinBounds (BRIDGESBIN *bridgeBetweenBins, BOUNDARY simBoundary, float binWidth, float delBinDistance, int nBins)
{
	omp_set_num_threads (NTHREADS);

	#pragma omp parallel for
	for (int i = 0; i < nBins; ++i)
	{
		if (i == 0)
		{
			bridgeBetweenBins[i].y1lo = simBoundary.yhi - binWidth;
			bridgeBetweenBins[i].y1hi = simBoundary.yhi;
			bridgeBetweenBins[i].y2lo = simBoundary.ylo;
			bridgeBetweenBins[i].y2hi = simBoundary.ylo + binWidth;
		}
		else
		{
			bridgeBetweenBins[i].y1lo = bridgeBetweenBins[i - 1].y1lo + delBinDistance;
			bridgeBetweenBins[i].y1hi = bridgeBetweenBins[i - 1].y1hi + delBinDistance;
			bridgeBetweenBins[i].y2lo = bridgeBetweenBins[i - 1].y2lo + delBinDistance;
			bridgeBetweenBins[i].y2hi = bridgeBetweenBins[i - 1].y2hi + delBinDistance;
		}

		bridgeBetweenBins[i].count = 0;
	}

	return bridgeBetweenBins;
}

BRIDGESBIN *countBridgesBetweenBins (TRAJECTORY **atoms, BOUNDARY simBoundary, float distanceCutoff, BRIDGESBIN *bridgeBetweenBins, int nAtoms, TRAJECTORY *micelles, int nMicelles, int nBins)
{
	omp_set_num_threads (NTHREADS);

	float distance1, distance2;
	float tempY1, tempY2;

	for (int i = 0; i < nAtoms; )
	{
		if ((*atoms)[i].atomType == 1)
		{
			#pragma omp parallel for
			for (int j = 0; j < nMicelles; ++j)
			{
				distance1 = computePeriodicDistance ((*atoms)[i].x, (*atoms)[i].y, (*atoms)[i].z, micelles[j].x, micelles[j].y, micelles[j].z, simBoundary.xLength, simBoundary.yLength, simBoundary.zLength);
				distance2 = computePeriodicDistance ((*atoms)[i + 1].x, (*atoms)[i + 1].y, (*atoms)[i + 1].z, micelles[j].x, micelles[j].y, micelles[j].z, simBoundary.xLength, simBoundary.yLength, simBoundary.zLength);

				if (distance1 <= distanceCutoff) {
					(*atoms)[i].adsorbedID = j; 
				}

				if (distance2 <= distanceCutoff) {
					(*atoms)[i + 1].adsorbedID = j; 
				}
			}
		}
		else if ((*atoms)[i].atomType == 2) {
			(*atoms)[i].adsorbedID = -1; }

		if ((*atoms)[i].atomType == 1) {
			i += 2; }
		else if ((*atoms)[i].atomType == 2) {
			i += 1; }
	}

	for (int i = 0; i < nAtoms; )
	{
		if (((*atoms)[i].adsorbedID != (*atoms)[i + 1].adsorbedID) && (*atoms)[i].adsorbedID != 0 && (*atoms)[i + 1].adsorbedID != 0 && (*atoms)[i].atomType == 1 && (*atoms)[i + 1].atomType == 1)
		{
			#pragma omp parallel for
			for (int k = 0; k < nBins; ++k)
			{
				tempY1 = translatePeriodic (bridgeBetweenBins[k].y1hi, (*atoms)[i].y, simBoundary.yLength);
				tempY2 = translatePeriodic (bridgeBetweenBins[k].y1hi, (*atoms)[i + 1].y, simBoundary.yLength);

				if ((tempY1 > bridgeBetweenBins[k].y1lo && 
					tempY1 <= bridgeBetweenBins[k].y1hi && 
					(*atoms)[i + 1].y > bridgeBetweenBins[k].y2lo && 
					(*atoms)[i + 1].y <= bridgeBetweenBins[k].y2hi)
					||
					(tempY2 > bridgeBetweenBins[k].y1lo && 
					tempY2 <= bridgeBetweenBins[k].y1hi && 
					(*atoms)[i].y > bridgeBetweenBins[k].y2lo && 
					(*atoms)[i].y <= bridgeBetweenBins[k].y2hi))
				{
					bridgeBetweenBins[k].count++;
				}
			}
		}

		if ((*atoms)[i].atomType == 1) {
			i += 2; }
		if ((*atoms)[i].atomType == 2) {
			i += 1; }
	}

	return bridgeBetweenBins;
}