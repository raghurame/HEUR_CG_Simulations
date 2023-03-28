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

BRIDGESBIN *assignBridgeCenterDistribution (BRIDGESBIN *bridgeCenterDistribution, int nBins_centerDistribution, float binWidth_centerDistribution, BOUNDARY simBoundary)
{
	omp_set_num_threads (NTHREADS);

	#pragma omp parallel for
	for (int i = 0; i < nBins_centerDistribution; ++i)
	{
		if (i == 0)
		{
			bridgeCenterDistribution[i].y1lo = simBoundary.ylo;
			bridgeCenterDistribution[i].y1hi = simBoundary.ylo + binWidth_centerDistribution;
		}
		else
		{
			bridgeCenterDistribution[i].y1lo = bridgeCenterDistribution[i - 1].y1hi;
			bridgeCenterDistribution[i].y1hi = bridgeCenterDistribution[i].y1lo + binWidth_centerDistribution;
		}

		bridgeCenterDistribution[i].count = 0;
	}

	return bridgeCenterDistribution;
}

BONDINFO *computeBridgeCenter (TRAJECTORY *atoms, int nAtoms, BONDINFO *allBonds, BOUNDARY simBoundary)
{
	int currentBondIndex = 0;
	float tempX, tempY, tempZ;
	float dotProduct, magnitude1, magnitude2, cosTheta;

	for (int i = 0; i < nAtoms; )
	{
		if (atoms[i].atomType == 1 && atoms[i + 1].atomType == 1)
		{
			tempX = translatePeriodic (atoms[i].x, atoms[i + 1].x, simBoundary.xLength);
			tempY = translatePeriodic (atoms[i].y, atoms[i + 1].y, simBoundary.yLength);
			tempZ = translatePeriodic (atoms[i].z, atoms[i + 1].z, simBoundary.zLength);

			allBonds[currentBondIndex].xc = (atoms[i].x + tempX) / 2;
			allBonds[currentBondIndex].yc = (atoms[i].y + tempY) / 2;
			allBonds[currentBondIndex].zc = (atoms[i].z + tempZ) / 2;

			allBonds[currentBondIndex].x1 = atoms[i].x;
			allBonds[currentBondIndex].y1 = atoms[i].y;
			allBonds[currentBondIndex].z1 = atoms[i].z;

			allBonds[currentBondIndex].x2 = tempX;
			allBonds[currentBondIndex].y2 = tempY;
			allBonds[currentBondIndex].z2 = tempZ;

			if (allBonds[currentBondIndex].xc > simBoundary.xhi) {
				allBonds[currentBondIndex].xc -= simBoundary.xLength; }
			else if (allBonds[currentBondIndex].xc < simBoundary.xlo) {
				allBonds[currentBondIndex].xc += simBoundary.xLength; }

			if (allBonds[currentBondIndex].yc > simBoundary.yhi) {
				allBonds[currentBondIndex].yc -= simBoundary.yLength; }
			else if (allBonds[currentBondIndex].yc < simBoundary.ylo) {
				allBonds[currentBondIndex].yc += simBoundary.yLength; }

			if (allBonds[currentBondIndex].zc > simBoundary.zhi) {
				allBonds[currentBondIndex].zc -= simBoundary.zLength; }
			else if (allBonds[currentBondIndex].zc < simBoundary.zlo) {
				allBonds[currentBondIndex].zc += simBoundary.zLength; }

			dotProduct = (allBonds[currentBondIndex].x2 - allBonds[currentBondIndex].x1) * (1 - 0);
			magnitude1 = sqrt (
				pow ((allBonds[currentBondIndex].x2 - allBonds[currentBondIndex].x1), 2) + 
				pow ((allBonds[currentBondIndex].y2 - allBonds[currentBondIndex].y1), 2) + 
				pow ((allBonds[currentBondIndex].z2 - allBonds[currentBondIndex].z1), 2));
			magnitude2 = 1;

			cosTheta = dotProduct / (magnitude1 * magnitude2);
			allBonds[currentBondIndex].xOrientationAngle = (acosf (cosTheta) * 180) / PI;

			i += 2;
			currentBondIndex++;
		}
		else {
			i += 1; }
	}

	return allBonds;
}

BRIDGESBIN *computeBridgeCenterDistribution (BONDINFO *allBonds, int nBonds, BRIDGESBIN *bridgeCenterDistribution, int nBins_centerDistribution)
{
	omp_set_num_threads (NTHREADS);

	#pragma omp parallel for
	for (int i = 0; i < nBonds; ++i)
	{
		for (int j = 0; j < nBins_centerDistribution; ++j)
		{
			if (allBonds[i].yc > bridgeCenterDistribution[j].y1lo && allBonds[i].yc <= bridgeCenterDistribution[j].y1hi) {
				bridgeCenterDistribution[j].count++; }
		}
	}

	return bridgeCenterDistribution;
}
