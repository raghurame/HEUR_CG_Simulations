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

ANGLE_DISTRIBUTION *assignBeadOrientation (ANGLE_DISTRIBUTION *beadOrientation, int nBins_orientation)
{
	for (int i = 0; i < nBins_orientation; ++i)
	{
		if (i == 0) {
			beadOrientation[i].angleLo = 0;
			beadOrientation[i].angleHi = BINWIDTH_ANGLEORIENTATION; }
		else
		{
			beadOrientation[i].angleLo = beadOrientation[i - 1].angleHi;
			beadOrientation[i].angleHi = beadOrientation[i].angleLo + BINWIDTH_ANGLEORIENTATION;
		}

		beadOrientation[i].count = 0;
	}

	return beadOrientation;
}

ANGLE_DISTRIBUTION *computeBridgeOrientationDistribution (ANGLE_DISTRIBUTION *bridgeOrientation, int nBins_orientation, int nBonds, BONDINFO *allBonds)
{
	for (int i = 0; i < nBonds; ++i)
	{
		if (allBonds[i].adsorbedID1 > 0 && allBonds[i].adsorbedID2 > 0 && allBonds[i].adsorbedID1 != allBonds[i].adsorbedID2)
		{
			for (int j = 0; j < nBins_orientation; ++j)
			{
				if (allBonds[i].xOrientationAngle > bridgeOrientation[j].angleLo && allBonds[i].xOrientationAngle <= bridgeOrientation[j].angleHi) {
					bridgeOrientation[j].count++; }
			}
		}
	}

	return bridgeOrientation;
}

ANGLE_DISTRIBUTION *computeLoopOrientationDistribution (ANGLE_DISTRIBUTION *loopOrientation, int nBins_orientation, int nBonds, BONDINFO *allBonds)
{
	for (int i = 0; i < nBonds; ++i)
	{
		if (allBonds[i].adsorbedID1 > 0 && allBonds[i].adsorbedID2 > 0 && allBonds[i].adsorbedID1 == allBonds[i].adsorbedID2)
		{
			for (int j = 0; j < nBins_orientation; ++j)
			{
				if (allBonds[i].xOrientationAngle > loopOrientation[j].angleLo && allBonds[i].xOrientationAngle <= loopOrientation[j].angleHi) {
					loopOrientation[j].count++; }
			}
		}
	}

	return loopOrientation;
}

ANGLE_DISTRIBUTION *computeDangleOrientationDistribution (ANGLE_DISTRIBUTION *dangleOrientation, int nBins_orientation, int nBonds, BONDINFO *allBonds)
{
	for (int i = 0; i < nBonds; ++i)
	{
		if ((allBonds[i].adsorbedID1 > 0 && allBonds[i].adsorbedID2 == 0) || (allBonds[i].adsorbedID1 == 0 && allBonds[i].adsorbedID2 > 0))
		{
			for (int j = 0; j < nBins_orientation; ++j)
			{
				if (allBonds[i].xOrientationAngle > dangleOrientation[j].angleLo && allBonds[i].xOrientationAngle <= dangleOrientation[j].angleHi) {
					dangleOrientation[j].count++; }
			}
		}
	}

	return dangleOrientation;
}

ANGLE_DISTRIBUTION *computeFreeOrientationDistribution (ANGLE_DISTRIBUTION *freeOrientation, int nBins_orientation, int nBonds, BONDINFO *allBonds)
{
	for (int i = 0; i < nBonds; ++i)
	{
		if (allBonds[i].adsorbedID1 == 0 && allBonds[i].adsorbedID2 == 0)
		{
			for (int j = 0; j < nBins_orientation; ++j)
			{
				if (allBonds[i].xOrientationAngle > freeOrientation[j].angleLo && allBonds[i].xOrientationAngle <= freeOrientation[j].angleHi) {
					freeOrientation[j].count++; }
			}
		}
	}

	return freeOrientation;
}

ANGLE_DISTRIBUTION *computeAverageOrientationDistribution (ANGLE_DISTRIBUTION *beadOrientation, int nBins_orientation, int nTimeframes)
{
	for (int i = 0; i < nBins_orientation; ++i)
	{
		beadOrientation[i].count /= nTimeframes;
	}

	return beadOrientation;
}

void printAverageOrientationDistribution (ANGLE_DISTRIBUTION *beadOrientation, int nBins_orientation, const char *filename_orientationDistribution)
{
	FILE *file_orientationDistribution;
	file_orientationDistribution = fopen (filename_orientationDistribution, "w");

	fprintf(file_orientationDistribution, "# angleLo, angleHi, avgCount\n");

	for (int i = 0; i < nBins_orientation; ++i)
	{
		fprintf(file_orientationDistribution, "%f %f %f\n", beadOrientation[i].angleLo, beadOrientation[i].angleHi, beadOrientation[i].count);
	}

	fclose (file_orientationDistribution);
}