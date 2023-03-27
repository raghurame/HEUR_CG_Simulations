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

BRIDGESBIN *assignBinBounds (BRIDGESBIN *bridgeBetweenBins, BOUNDARY simBoundary, float binWidth, float delBinDistance, int nBins)
{
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
	float distance1, distance2;
	float tempX, tempY, tempY1, tempY2, tempZ;

	for (int i = 0; i < nAtoms; )
	{
		if ((*atoms)[i].atomType == 1)
		{
			#pragma omp parallel for
			for (int j = 0; j < (nMicelles); ++j)
			{
				if (fabs ((*atoms)[i].x - micelles[j].x) > (simBoundary.xLength / 2))
				{
					if ((*atoms)[i].x > micelles[j].x) {
						tempX = micelles[j].x + simBoundary.xLength; }
					else if ((*atoms)[i].x < micelles[j].x) {
						tempX = micelles[j].x - simBoundary.xLength; }
				}
				else {
					tempX = (*atoms)[i].x; }

				if (fabs ((*atoms)[i].y - micelles[j].y) > (simBoundary.yLength / 2))
				{
					if ((*atoms)[i].y > micelles[j].y) {
						tempY = micelles[j].y + simBoundary.yLength; }
					else if ((*atoms)[i].y < micelles[j].y) {
						tempY = micelles[j].y - simBoundary.yLength; }
				}
				else {
					tempY = (*atoms)[i].y; }

				if (fabs ((*atoms)[i].z - micelles[j].z) > (simBoundary.zLength / 2))
				{
					if ((*atoms)[i].z > micelles[j].z) {
						tempZ = micelles[j].z + simBoundary.zLength; }
					else if ((*atoms)[i].z < micelles[j].z) {
						tempZ = micelles[j].z - simBoundary.zLength; }
				}
				else {
					tempZ = (*atoms)[i].z; }

				distance1 = sqrt (
					pow (((*atoms)[i].x - tempX), 2) +
					pow (((*atoms)[i].y - tempY), 2) +
					pow (((*atoms)[i].z - tempZ), 2)
					);

				distance2 = sqrt (
					pow (((*atoms)[i + 1].x - tempX), 2) +
					pow (((*atoms)[i + 1].y - tempY), 2) +
					pow (((*atoms)[i + 1].z - tempZ), 2)
					);

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
				if (fabs ((*atoms)[i].y - bridgeBetweenBins[k].y1hi) > (simBoundary.yLength / 2))
				{
					if ((*atoms)[i].y < bridgeBetweenBins[k].y1hi) {
						tempY1 = (*atoms)[i].y + simBoundary.yLength; }
					else if ((*atoms)[i].y > bridgeBetweenBins[k].y1hi) {
						tempY1 = (*atoms)[i].y - simBoundary.yLength; }
				}
				else {
					tempY1 = (*atoms)[i].y; }

				if (fabs ((*atoms)[i + 1].y - bridgeBetweenBins[k].y1hi) > (simBoundary.yLength / 2))
				{
					if ((*atoms)[i + 1].y < bridgeBetweenBins[k].y1hi) {
						tempY2 = (*atoms)[i + 1].y + simBoundary.yLength; }
					else if ((*atoms)[i + 1].y > bridgeBetweenBins[k].y1hi) {
						tempY2 = (*atoms)[i + 1].y - simBoundary.yLength; }
				}
				else {
					tempY2 = (*atoms)[i + 1].y; }

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