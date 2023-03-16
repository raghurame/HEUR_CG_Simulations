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
	float xLength, yLength, zLength;
} BOUNDARY;

typedef struct trajectory
{
	int sino, atomType, ix, iy, iz;
	float x, y, z;
	int adsorbedID;
} TRAJECTORY;

typedef struct bondIndo
{
	float x1, y1, z1, x2, y2, z2;
	float xc, yc, zc;
} BONDINFO;

typedef struct yDistribution
{
	float ylo, yhi;
	int count;
} YDIST;

typedef struct brigdesBin
{
	int count;
	float y1lo, y1hi, y2lo, y2hi;
} BRIDGESBIN;

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

	simBoundary.xLength = simBoundary.xhi - simBoundary.xlo;
	simBoundary.yLength = simBoundary.yhi - simBoundary.ylo;
	simBoundary.zLength = simBoundary.zhi - simBoundary.zlo;

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
			for (int j = 0; j < (nMicelles * 2); ++j)
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

YDIST *assignBridgeYDistribution (float maxFeneExtension, int nBins, float binWidth, YDIST *bridgeYDistribution)
{
	for (int i = 0; i < nBins; ++i)
	{
		if (i == 0) {
			bridgeYDistribution[i].ylo = 0; }
		else {
			bridgeYDistribution[i].ylo = bridgeYDistribution[i - 1].yhi; }

		bridgeYDistribution[i].yhi = bridgeYDistribution[i].ylo + binWidth;
		bridgeYDistribution[i].count = 0;
	}

	return bridgeYDistribution;
}

YDIST *computeBridgeDistribution (TRAJECTORY *atoms, int nAtoms, YDIST *bridgeYDistribution, int nBins, BOUNDARY simBoundary)
{
	int bridgeCountLocal = 0;
	float yDistance;
	float tempX, tempY, tempZ;

	for (int i = 0; i < nAtoms; )
	{
		if ((atoms[i].adsorbedID != atoms[i + 1].adsorbedID) && atoms[i].adsorbedID != 0 && atoms[i + 1].adsorbedID != 0)
		{
			if (fabs (atoms[i].y - atoms[i + 1].y) > (simBoundary.yLength / 2))
			{
				if (atoms[i + 1].y > atoms[i].y) {
					tempY = atoms[i + 1].y - simBoundary.yLength;
					yDistance = fabs (atoms[i].y - tempY); }
				else if (atoms[i + 1].y < atoms[i].y) {
					tempY = atoms[i + 1].y + simBoundary.yLength;
					yDistance = fabs (atoms[i].y - tempY); }
			}
			else {
				yDistance = fabs (atoms[i].y - atoms[i + 1].y); }

			#pragma omp parallel for
			for (int j = 0; j < nBins; ++j)
			{
				if (yDistance > bridgeYDistribution[j].ylo && yDistance <= bridgeYDistribution[j].yhi) {
					bridgeYDistribution[j].count++; }
			}
		}

		if ((atoms[i].adsorbedID != atoms[i + 1].adsorbedID) && atoms[i].atomType == 1 && atoms[i + 1].atomType == 1 && atoms[i].adsorbedID > 0 && atoms[i + 1].adsorbedID > 0) {
			bridgeCountLocal++; }

		if (atoms[i].atomType == 1) {
			i += 2; }
		else if (atoms[i].atomType == 2) {
			i += 1; }
	}

	return bridgeYDistribution;
}

BONDINFO *computeBridgeCenter (TRAJECTORY *atoms, int nAtoms, BONDINFO *allBonds, BOUNDARY simBoundary)
{
	int currentBondIndex = 0;
	float tempX, tempY, tempZ;

	for (int i = 0; i < nAtoms; ++i)
	{
		if (atoms[i].atomType == 1 && atoms[i + 1].atomType == 1)
		{
			if (fabs (atoms[i].x - atoms[i + 1].x) > (simBoundary.xLength / 2))
			{
				if (atoms[i].x > atoms[i + 1].x) {
					tempX = atoms[i + 1].x + simBoundary.xLength; }
				else {
					tempX = atoms[i + 1].x - simBoundary.xLength; }					
			}
			else {
				tempX = atoms[i + 1].x; }

			if (fabs (atoms[i].y - atoms[i + 1].y) > (simBoundary.yLength / 2))
			{
				if (atoms[i].y > atoms[i + 1].y) {
					tempY = atoms[i + 1].y + simBoundary.yLength; }
				else {
					tempY = atoms[i + 1].y - simBoundary.yLength; }
			}
			else {
				tempY = atoms[i + 1].y; }

			if (fabs (atoms[i].z - atoms[i + 1].z) > (simBoundary.zLength / 2))
			{
				if (atoms[i].z > atoms[i + 1].z) {
					tempZ = atoms[i + 1].z + simBoundary.zLength; }
				else {
					tempZ = atoms[i + 1].z - simBoundary.zLength; }
			}
			else {
				tempZ = atoms[i + 1].z; }

			allBonds[currentBondIndex].xc = (atoms[i].x + tempX) / 2;
			allBonds[currentBondIndex].yc = (atoms[i].y + tempY) / 2;
			allBonds[currentBondIndex].zc = (atoms[i].z + tempZ) / 2;

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

			i += 2;
			currentBondIndex++;
		}
		else
		{
			i += 1;
		}
	}

	return allBonds;
}

BRIDGESBIN *assignBridgeCenterDistribution (BRIDGESBIN *bridgeCenterDistribution, int nBins_centerDistribution, float binWidth_centerDistribution, BOUNDARY simBoundary)
{
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

BRIDGESBIN *computeBridgeCenterDistribution (BONDINFO *allBonds, int nBonds, BRIDGESBIN *bridgeCenterDistribution, int nBins_centerDistribution)
{
	#pragma omp parallel for
	for (int i = 0; i < nBonds; ++i)
	{
		for (int j = 0; j < nBins_centerDistribution; ++j)
		{
			if (allBonds[i].yc >= bridgeCenterDistribution[j].y1lo && allBonds[i].yc < bridgeCenterDistribution[j].y1hi) {
				bridgeCenterDistribution[j].count++; }
		}
	}

	return bridgeCenterDistribution;
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

	BRIDGESBIN *bridgeBetweenBins;
	bridgeBetweenBins = (BRIDGESBIN *) malloc (nBins_vertBridges * sizeof (BRIDGESBIN));

	TRAJECTORY *atoms, *micelles;
	atoms = (TRAJECTORY *) malloc (nAtoms * 2 * sizeof (TRAJECTORY));
	micelles = (TRAJECTORY *) malloc (nMicelles * 2 * sizeof (TRAJECTORY));

	float maxFeneExtension = 60.0, binWidth_yDist = (maxFeneExtension / 20);
	int nBins_yDist = 20; //Taken arbitrarily for now; 20 bins across the maximum extensible length of 60 sigma
	YDIST *bridgeYDistribution;
	bridgeYDistribution = (YDIST *) malloc (nBins_yDist * sizeof (YDIST));

	int nBonds = (nAtoms - nMicelles) / 2, nBins_bridgeCenter = 20;
	BONDINFO *allBonds;
	allBonds = (BONDINFO *) malloc (nBonds * sizeof (BONDINFO));
	BRIDGESBIN *bridgeCenterDistribution;

	float binWidth_centerDistribution = 3;
	int nBins_centerDistribution = ceil (simBoundary.yLength / binWidth_centerDistribution);
	bridgeCenterDistribution = (BRIDGESBIN *) malloc (nBins_centerDistribution * sizeof (BRIDGESBIN));

	bridgeBetweenBins = assignBinBounds (bridgeBetweenBins, simBoundary, binWidth_vertBridges, delBinDistance_vertBridges, nBins_vertBridges);
	bridgeYDistribution = assignBridgeYDistribution (maxFeneExtension, nBins_yDist, binWidth_yDist, bridgeYDistribution);
	bridgeCenterDistribution = assignBridgeCenterDistribution (bridgeCenterDistribution, nBins_centerDistribution, binWidth_centerDistribution, simBoundary);

	while (file_status != EOF)
	{
		// if (nTimeframes%1 == 0) {
			// fprintf(stdout, "computing %d timesteps...         \r", nTimeframes);
			// fflush (stdout); }

		atoms = getAtoms (atoms, nAtoms, simBoundary, distanceCutoff_vertBridges, file_inputTrj, file_status, &micelles, nMicelles);

		bridgeBetweenBins = countBridgesBetweenBins (&atoms, simBoundary, distanceCutoff_vertBridges, bridgeBetweenBins, nAtoms, micelles, nMicelles, nBins_vertBridges);
		bridgeYDistribution = computeBridgeDistribution (atoms, nAtoms, bridgeYDistribution, nBins_yDist, simBoundary);

		// Finish the bridge center distribution calculations
		allBonds = computeBridgeCenter (atoms, nAtoms, allBonds, simBoundary);
		bridgeCenterDistribution = computeBridgeCenterDistribution (allBonds, nBonds, bridgeCenterDistribution, nBins_centerDistribution);

		// Plot the distribution of angles between all the bonds and an unit vector along X-axis (along the velocity profile direction.)

		file_status = fgetc (file_inputTrj);
		nTimeframes++;
	}

	FILE *file_bridgeBetweenBinsOuptut, *file_bridgeYDistributionOutput;
	file_bridgeBetweenBinsOuptut = fopen ("nBridgesBetweenBins.count", "w");
	file_bridgeYDistributionOutput = fopen ("bridges.distribution", "w");

	fprintf(file_bridgeBetweenBinsOuptut, "# y1lo, y1hi, y2lo, y2hi, avgCounts\n");
	fprintf(file_bridgeYDistributionOutput, "# ylo, yhi, avgCounts\n");

	for (int i = 0; i < nBins_vertBridges; ++i) {
		fprintf(file_bridgeBetweenBinsOuptut, "%f %f %f %f %f\n", 
			bridgeBetweenBins[i].y1lo, 
			bridgeBetweenBins[i].y1hi, 
			bridgeBetweenBins[i].y2lo, 
			bridgeBetweenBins[i].y2hi, 
			((float)bridgeBetweenBins[i].count / (float)nTimeframes)); }

	for (int i = 0; i < nBins_yDist; ++i) {
		fprintf(file_bridgeYDistributionOutput, "%f %f %f\n", bridgeYDistribution[i].ylo, bridgeYDistribution[i].yhi, ((float)bridgeYDistribution[i].count / (float)nTimeframes)); }

	fclose (file_inputTrj);
	return 0;
}