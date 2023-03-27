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

BOUNDARY getBoundary (FILE *file_inputTrj, BOUNDARY simBoundary)
{
	char lineString[2000];

	for (int i = 0; i < 5; ++i)
	{
		fgets (lineString, 2000, file_inputTrj);
	}

	fgets (lineString, 2000, file_inputTrj);
	sscanf (lineString, "%f %f %f\n", &simBoundary.xlo, &simBoundary.xhi, &simBoundary.xy);

	fgets (lineString, 2000, file_inputTrj);
	sscanf (lineString, "%f %f %f\n", &simBoundary.ylo, &simBoundary.yhi, &simBoundary.xz);

	fgets (lineString, 2000, file_inputTrj);
	sscanf (lineString, "%f %f %f\n", &simBoundary.zlo, &simBoundary.zhi, &simBoundary.yz);

	rewind (file_inputTrj);

	simBoundary.xLength = simBoundary.xhi - simBoundary.xlo;
	simBoundary.yLength = simBoundary.yhi - simBoundary.ylo;
	simBoundary.zLength = simBoundary.zhi - simBoundary.zlo;

	return simBoundary;
}

int countNAtoms (FILE *file_inputTrj)
{
	int nAtoms;
	char lineString[2000];

	for (int i = 0; i < 4; ++i)
	{
		fgets (lineString, 2000, file_inputTrj);
	}

	sscanf (lineString, "%d\n", &nAtoms);
	rewind (file_inputTrj);

	return nAtoms;
}

TRAJECTORY *getAtoms (TRAJECTORY *atoms, int nAtoms, BOUNDARY simBoundary, float distanceCutoff, FILE *file_inputTrj, int file_status, TRAJECTORY **micelles, int nMicelles)
{
	char lineString[2000];
	float xLength, yLength, zLength;
	int currentMicelle = 0;

	for (int i = 0; i < 9; ++i)
	{
		fgets (lineString, 2000, file_inputTrj);

		if (i == 5) {
			sscanf (lineString, "%f %f %f\n", &simBoundary.xlo, &simBoundary.xhi, &simBoundary.xy); }
		else if (i == 6) {
			sscanf (lineString, "%f %f %f\n", &simBoundary.ylo, &simBoundary.yhi, &simBoundary.xz); }
		else if (i == 7) {
			sscanf (lineString, "%f %f %f\n", &simBoundary.zlo, &simBoundary.zhi, &simBoundary.yz); }

		if (i == 8) {
			xLength = (simBoundary.xhi - simBoundary.xlo);
			yLength = (simBoundary.yhi - simBoundary.ylo);
			zLength = (simBoundary.zhi - simBoundary.zlo); }
	}

	for (int i = 0; i < nAtoms; ++i)
	{
		fgets (lineString, 2000, file_inputTrj);
		// sscanf (lineString, "%d %d %f %f %f %d %d %d\n", atoms[i].sino, atoms[i].atomType, atoms[i].x, atoms[i].y, atoms[i].z, atoms[i].ix, atoms[i].iy, atoms[i].iz);
		sscanf (lineString, "%d %d %f %f %f\n", 
			&atoms[i].sino, 
			&atoms[i].atomType, 
			&atoms[i].x, 
			&atoms[i].y, 
			&atoms[i].z);

		atoms[i].adsorbedID = 0;

		if (atoms[i].atomType == 2)
		{
			(*micelles)[currentMicelle].sino = atoms[i].sino;
			(*micelles)[currentMicelle].atomType = atoms[i].atomType;
			(*micelles)[currentMicelle].x = atoms[i].x;
			(*micelles)[currentMicelle].y = atoms[i].y;
			(*micelles)[currentMicelle].z = atoms[i].z;
			(*micelles)[currentMicelle].adsorbedID = -1;

			currentMicelle++;
		}
	}

	return atoms;
}

int countNMicelles (FILE *file_inputTrj, int nAtoms)
{
	rewind (file_inputTrj);
	char lineString[2000];
	int atomType, nMicelles = 0;

	for (int i = 0; i < 9; ++i)
	{
		fgets (lineString, 2000, file_inputTrj);
	}
	for (int i = 0; i < nAtoms; ++i)
	{
		fgets (lineString, 2000, file_inputTrj);
		sscanf (lineString, "%*d %d\n", &atomType);

		if (atomType == 2)
		{
			nMicelles++;
		}
	}

	return nMicelles;
}
