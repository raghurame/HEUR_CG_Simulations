#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <omp.h>

typedef struct boundary
{
	float xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz;
} BOUNDARY;

typedef struct trajectory
{
	int sino, atomType, ix, iy, iz;
	float x, y, z;
	int adsorbedID;
} TRAJECTORY;

typedef struct brigdes
{
	int count;
	float y1lo, y1hi, y2lo, y2hi;
} BRIDGES;

typedef struct yDistribution
{
	float ylo, yhi;
	int count;
} YDIST;

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

	return simBoundary;
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

		atoms[i + nAtoms].sino = atoms[i].sino + nAtoms;
		atoms[i + nAtoms].atomType = atoms[i].atomType;
		atoms[i + nAtoms].x = atoms[i].x;
		atoms[i + nAtoms].y = atoms[i].y - yLength; // translating along negative Y direction
		atoms[i + nAtoms].z = atoms[i].z;
		// atoms[i + nAtoms].ix = atoms[i].ix;
		// atoms[i + nAtoms].iy = atoms[i].iy;
		// atoms[i + nAtoms].iz = atoms[i].iz;
		atoms[i + nAtoms].adsorbedID = 0;

		if (atoms[i].atomType == 2)
		{
			(*micelles)[currentMicelle].sino = atoms[i].sino;
			(*micelles)[currentMicelle].atomType = atoms[i].atomType;
			(*micelles)[currentMicelle].x = atoms[i].x;
			(*micelles)[currentMicelle].y = atoms[i].y;
			(*micelles)[currentMicelle].z = atoms[i].z;
			(*micelles)[currentMicelle].adsorbedID = -1;

			(*micelles)[currentMicelle + nMicelles].sino = atoms[i].sino;
			(*micelles)[currentMicelle + nMicelles].atomType = atoms[i].atomType;
			(*micelles)[currentMicelle + nMicelles].x = atoms[i].x;
			(*micelles)[currentMicelle + nMicelles].y = atoms[i].y - yLength;
			(*micelles)[currentMicelle + nMicelles].z = atoms[i].z;
			(*micelles)[currentMicelle + nMicelles].adsorbedID = -1;

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

BRIDGES *assignBinBounds (BRIDGES *bridgeBetweenBins, BOUNDARY simBoundary, float binWidth, float delBinDistance, int nBins)
{
	for (int i = 0; i < nBins; ++i)
	{
		if (i == 0)
		{
			bridgeBetweenBins[i].y1lo = simBoundary.ylo - (simBoundary.yhi - simBoundary.ylo);
		}
		else
		{
			bridgeBetweenBins[i].y1lo = bridgeBetweenBins[i - 1].y1lo + delBinDistance;
		}

		bridgeBetweenBins[i].y1hi = bridgeBetweenBins[i].y1lo + binWidth;
		bridgeBetweenBins[i].y2lo = bridgeBetweenBins[i].y1hi;
		bridgeBetweenBins[i].y2hi = bridgeBetweenBins[i].y2lo + binWidth;

		bridgeBetweenBins[i].count = 0;
	}

	return bridgeBetweenBins;
}

BRIDGES *countBridgesBetweenBins (TRAJECTORY **atoms, BOUNDARY simBoundary, float distanceCutoff, BRIDGES *bridgeBetweenBins, int nAtoms, TRAJECTORY *micelles, int nMicelles, int nBins)
{
	float distance1, distance2;

	for (int i = 0; i < (nAtoms * 2); )
	{
		if ((*atoms)[i].atomType == 1)
		{
			for (int j = 0; j < (nMicelles * 2); ++j)
			{
				distance1 = sqrt (
					pow (((*atoms)[i].x - micelles[j].x), 2) +
					pow (((*atoms)[i].y - micelles[j].y), 2) +
					pow (((*atoms)[i].z - micelles[j].z), 2)
					);

				distance2 = sqrt (
					pow (((*atoms)[i + 1].x - micelles[j].x), 2) +
					pow (((*atoms)[i + 1].y - micelles[j].y), 2) +
					pow (((*atoms)[i + 1].z - micelles[j].z), 2)
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

	for (int i = 0; i < (nAtoms * 2); )
	{
		if (((*atoms)[i].adsorbedID != (*atoms)[i + 1].adsorbedID) && (*atoms)[i].adsorbedID != 0 && (*atoms)[i + 1].adsorbedID != 0 && (*atoms)[i].atomType == 1 && (*atoms)[i + 1].atomType == 1)
		{
			for (int k = 0; k < nBins; ++k)
			{
				if (((*atoms)[i].y > bridgeBetweenBins[k].y1lo && 
					(*atoms)[i].y <= bridgeBetweenBins[k].y1hi && 
					(*atoms)[i + 1].y > bridgeBetweenBins[k].y2lo && 
					(*atoms)[i + 1].y <= bridgeBetweenBins[k].y2hi)
					||
					((*atoms)[i + 1].y > bridgeBetweenBins[k].y1lo && 
					(*atoms)[i + 1].y <= bridgeBetweenBins[k].y1hi && 
					(*atoms)[i].y > bridgeBetweenBins[k].y2lo && 
					(*atoms)[i].y <= bridgeBetweenBins[k].y2hi))
				{
					bridgeBetweenBins[k].count++;
				}
			}
		}

		i += 2;
	}

	return bridgeBetweenBins;
}

YDIST *assignBridgeDistribution (float maxFeneExtension, int nBins, float binWidth, YDIST *bridgeDistribution)
{
	for (int i = 0; i < nBins; ++i)
	{
		if (i == 0) {
			bridgeDistribution[i].ylo = 0; }
		else {
			bridgeDistribution[i].ylo = bridgeDistribution[i - 1].yhi; }

		bridgeDistribution[i].yhi = bridgeDistribution[i].ylo + binWidth;
		bridgeDistribution[i].count = 0;
	}

	return bridgeDistribution;
}

YDIST *computeBridgeDistribution (TRAJECTORY *atoms, int nAtoms, YDIST *bridgeDistribution, int nBins)
{
	int bridgeCountLocal = 0;

	for (int i = 0; i < nAtoms; )
	{
		if ((atoms[i].adsorbedID != atoms[i + 1].adsorbedID) && atoms[i].adsorbedID != 0 && atoms[i + 1].adsorbedID != 0)
		{
			for (int j = 0; j < nBins; ++j)
			{
				if (abs (atoms[i].y - atoms[i + 1].y) > bridgeDistribution[j].ylo && abs (atoms[i].y - atoms[i + 1].y) <= bridgeDistribution[j].yhi)
				{
					bridgeDistribution[j].count++;
				}
			}
		}

		if ((atoms[i].adsorbedID != atoms[i + 1].adsorbedID) && atoms[i].atomType == 1 && atoms[i + 1].atomType == 1 && atoms[i].adsorbedID > 0 && atoms[i + 1].adsorbedID > 0)
		{
			bridgeCountLocal++;
		}

		if (atoms[i].atomType == 1) {
			i += 2; }
		else if (atoms[i].atomType == 2) {
			i += 1; }
	}

	return bridgeDistribution;
}

int main(int argc, char const *argv[])
{
	if (argc != 5)
	{
		printf("\nINCORRECT ARGUMENTS PASSED:\n~~~~~~~~~~~~~~~~~~~~~~~~~\n\n -> Pass the input file as (char *) argv[1]\n -> Enter bin width as (float) argv[2]\n -> Enter the distance cut off for particle-bead as (float) argv[3]\n -> Enter the increment in bin distance as (float) argv[4]\n\nFor example: ./countBridges last.lammpstrj 5.0 3.38 1.0\n~~~~~~~~~~~\n\nTo compile: gcc -o countBridges countBridges.c -lm\n~~~~~~~~~~\n\n");
		exit (1);
	}

	FILE *file_inputTrj;
	file_inputTrj = fopen (argv[1], "r");

	int nAtoms = countNAtoms (file_inputTrj), file_status = 0, nMicelles = countNMicelles (file_inputTrj, nAtoms), nTimeframes = 0;
	BOUNDARY simBoundary;
	simBoundary = getBoundary (file_inputTrj, simBoundary);

	char lineString[2000];
	float binWidth_vertBridges = atof (argv[2]), distanceCutoff_vertBridges = atof (argv[3]), delBinDistance_vertBridges = atof (argv[4]);
	int nBins_vertBridges = ceil ((((simBoundary.yhi - simBoundary.ylo) * 2) - (2 * binWidth_vertBridges)) / delBinDistance_vertBridges);

	BRIDGES *bridgeBetweenBins;
	bridgeBetweenBins = (BRIDGES *) malloc (nBins_vertBridges * sizeof (BRIDGES));

	TRAJECTORY *atoms, *micelles;
	atoms = (TRAJECTORY *) malloc (nAtoms * 2 * sizeof (TRAJECTORY));
	micelles = (TRAJECTORY *) malloc (nMicelles * 2 * sizeof (TRAJECTORY));

	float maxFeneExtension = 60.0, binWidth_yDist = (maxFeneExtension / 20);
	int nBins_yDist = 20; //Taken arbitrarily for now; 20 bins across the maximum extensible length of 60 sigma
	YDIST *bridgeDistribution;
	bridgeDistribution = (YDIST *) malloc (nBins_yDist * sizeof (YDIST));

	bridgeBetweenBins = assignBinBounds (bridgeBetweenBins, simBoundary, binWidth_vertBridges, delBinDistance_vertBridges, nBins_vertBridges);
	bridgeDistribution = assignBridgeDistribution (maxFeneExtension, nBins_yDist, binWidth_yDist, bridgeDistribution);

	while (file_status != EOF)
	{
		if (nTimeframes%1 == 0) {
			fprintf(stdout, "computing %d timesteps...         \r", nTimeframes);
			fflush (stdout); }

		atoms = getAtoms (atoms, nAtoms, simBoundary, distanceCutoff_vertBridges, file_inputTrj, file_status, &micelles, nMicelles);
		bridgeBetweenBins = countBridgesBetweenBins (&atoms, simBoundary, distanceCutoff_vertBridges, bridgeBetweenBins, nAtoms, micelles, nMicelles, nBins_vertBridges);
		bridgeDistribution = computeBridgeDistribution (atoms, nAtoms, bridgeDistribution, nBins_yDist);

		file_status = fgetc (file_inputTrj);
		nTimeframes++;
	}

	FILE *file_bridgeBetweenBinsOuptut, *file_bridgeDistributionOutput;
	file_bridgeBetweenBinsOuptut = fopen ("nBridgesBetweenBins.count", "w");
	file_bridgeDistributionOutput = fopen ("bridges.distribution", "w");

	fprintf(file_bridgeBetweenBinsOuptut, "# y1lo, y1hi, y2lo, y2hi, avgCounts\n");
	fprintf(file_bridgeDistributionOutput, "# ylo, yhi, avgCounts\n");

	for (int i = 0; i < nBins_vertBridges; ++i) {
		fprintf(file_bridgeBetweenBinsOuptut, "%f %f %f %f %f\n", bridgeBetweenBins[i].y1lo, bridgeBetweenBins[i].y1hi, bridgeBetweenBins[i].y2lo, bridgeBetweenBins[i].y2hi, ((float)bridgeBetweenBins[i].count / (float)nTimeframes)); }

	for (int i = 0; i < nBins_yDist; ++i) {
		fprintf(file_bridgeDistributionOutput, "%f %f %f\n", bridgeDistribution[i].ylo, bridgeDistribution[i].yhi, ((float)bridgeDistribution[i].count / (float)nTimeframes)); }

	fclose (file_inputTrj);
	return 0;
}