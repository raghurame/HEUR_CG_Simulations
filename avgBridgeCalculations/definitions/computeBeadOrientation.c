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

ANGLE_DISTRIBUTION *computeBeadOrientationDistribution (TRAJECTORY *atoms, int nBins_orientation, BONDINFO *allBonds, int nBonds, ANGLE_DISTRIBUTION **bridgeOrientation, ANGLE_DISTRIBUTION **loopOrientation, ANGLE_DISTRIBUTION **dangleOrientation, ANGLE_DISTRIBUTION **loopOrientation)
{
	for (int i = 0; i < nBonds; ++i)
	{
		// allBonds[currentBondIndex].xOrientationAngle
		for (int j = 0; j < nBins_orientation; ++j)
		{
			if (allBonds[i].xOrientationAngle > beadOrientation[j].angleLo && allBonds[i].xOrientationAngle <= beadOrientation[j].angleHi) 
			{
				
			}
		}
	}

	return beadOrientation;
}