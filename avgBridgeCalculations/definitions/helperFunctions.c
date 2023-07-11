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

float computePeriodicDistance (float x1, float y1, float z1, float x2, float y2, float z2, BOUNDARY simBoundary)
{
	float distance;
	float x2_new, y2_new, z2_new;

	x2_new = translatePeriodic (x1, x2, simBoundary.xLength);
	y2_new = translatePeriodic (y1, y2, simBoundary.yLength);
	z2_new = translatePeriodic (z1, z2, simBoundary.zLength);

	// simBoundary.xy, simBoundary.xz, simBoundary.yz
	if (y2_new > y2)
	{
		x2_new += simBoundary.xy;
	}
	else if (y2_new < y2)
	{
		x2_new -= simBoundary.xy;
	}

	distance = sqrt (pow ((x2_new - x1), 2) + pow ((y2_new - y1), 2) + pow ((z2_new - z1), 2));

	return distance;
}
