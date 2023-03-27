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

float computePeriodicDistance (float x1, float y1, float z1, float x2, float y2, float z2, float simulationBoxLength)
{
	float distance;

	x2 = translatePeriodic (x1, x2, simulationBoxLength);
	y2 = translatePeriodic (y1, y2, simulationBoxLength);
	z2 = translatePeriodic (z1, z2, simulationBoxLength);

	distance = sqrt (pow ((x2 - x1), 2) + pow ((y2 - y1), 2) + pow ((z2 - z1), 2));

	return distance;
}
