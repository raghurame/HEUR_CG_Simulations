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

float translatePeriodic (float r1, float r2, float simulationBoxLength)
{
	if (fabs (r1 - r2) > (simulationBoxLength / 2))
	{
		if (r1 > r2) {
			r2 += simulationBoxLength; }
		else if (r2 > r1) {
			r2 -= simulationBoxLength; }
	}

	return r2;
}

float computePeriodicDistance (float x1, float y1, float z1, float x2, float y2, float z2, float xLength, float yLength, float zLength)
{
	float distance;

	x2 = translatePeriodic (x1, x2, xLength);
	y2 = translatePeriodic (y1, y2, yLength);
	z2 = translatePeriodic (z1, z2, zLength);

	distance = sqrt (pow ((x2 - x1), 2) + pow ((y2 - y1), 2) + pow ((z2 - z1), 2));

	return distance;
}
