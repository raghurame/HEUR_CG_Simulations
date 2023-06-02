#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

#define BITS sizeof(float) * 8

typedef struct boundary
{
	float xlo, xhi, ylo, yhi, zlo, zhi;
} BOUNDARY;

typedef struct dumpEntries
{
	int id, type;
	float x, y, z;
} DUMP_ENTRIES;

int countNAtoms (int nAtoms, FILE *input)
{
	char lineString[2000];

	for (int i = 0; i < 4; ++i)
	{
		fgets (lineString, 2000, input);
	}

	sscanf (lineString, "%d\n", &nAtoms);
	rewind (input);

	return nAtoms;
}

BOUNDARY computeSimBoundary (BOUNDARY simBoundary, FILE *input)
{
	char lineString[2000];

	for (int i = 0; i < 5; ++i)
	{
		fgets (lineString, 2000, input);
	}

	fgets (lineString, 2000, input);
	sscanf (lineString, "%f %f\n", &simBoundary.xlo, &simBoundary.xhi);
	fgets (lineString, 2000, input);
	sscanf (lineString, "%f %f\n", &simBoundary.ylo, &simBoundary.yhi);
	fgets (lineString, 2000, input);
	sscanf (lineString, "%f %f\n", &simBoundary.zlo, &simBoundary.zhi);

	rewind (input);

	return simBoundary;
}

DUMP_ENTRIES *readDump (DUMP_ENTRIES *atoms, int nAtoms, FILE *input)
{
	char lineString[2000];

	for (int i = 0; i < 9; ++i)
	{
		fgets (lineString, 2000, input);
	}

	for (int i = 0; i < nAtoms; ++i)
	{
		fgets (lineString, 2000, input);
		sscanf (lineString, "%d %d %f %f %f\n", &atoms[i].id, &atoms[i].type, &atoms[i].x, &atoms[i].y, &atoms[i].z);
	}

	return atoms;
}

float mic (float r1, float r2, float boxLength)
{
	if (fabs (r2 - r1) > (boxLength / 2))
	{
		if (r1 > r2) {
			r1 -= boxLength; }
		else {
			r1 += boxLength; }
	}

	return r1;
}

float computeDistance (DUMP_ENTRIES *atoms, int i, int j, BOUNDARY simBoundary, float distance)
{
	float newX1, newY1, newZ1;

	float xLength = (simBoundary.xhi - simBoundary.xlo), yLength = (simBoundary.yhi - simBoundary.ylo), zLength = (simBoundary.zhi - simBoundary.zlo);

	newX1 = mic (atoms[i].x, atoms[j].x, xLength);
	newY1 = mic (atoms[i].y, atoms[j].y, yLength);
	newZ1 = mic (atoms[i].z, atoms[j].z, zLength);

	distance = sqrt (
		(atoms[j].x - newX1) * (atoms[j].x - newX1) +
		(atoms[j].y - newY1) * (atoms[j].y - newY1) +
		(atoms[j].z - newZ1) * (atoms[j].z - newZ1)
		);

	return distance;
}

int findCoordinationNumber (int coordinationNumber, DUMP_ENTRIES *atoms, int nAtoms, float cutoffDistance, int atomType1, int atomType2, BOUNDARY simBoundary)
{
	float distance;
	coordinationNumber = 0;

	for (int i = 0; i < nAtoms; ++i)
	{
		if (atoms[i].type == atomType1)
		{
			for (int j = (i + 1); j < nAtoms; ++j)
			{
				if (atoms[j].type == atomType2)
				{
					distance = computeDistance (atoms, i, j, simBoundary, distance);

					if (distance < cutoffDistance) {
						coordinationNumber++; }
				}
			}
		}
	}

	return coordinationNumber;
}

int main(int argc, char const *argv[])
{
	if (argc != 5)
	{
		printf("REQUIRED ARGUMENTS:\n~~~~~~~~~~~~~~~~~\n\n{~} argv[0] = program\n{~} argv[1] = input filename\n{~} argv[2] = atom type 1\n{~} argv[3] = atom type 2\n{~} argv[4] = cut-off distance\n\n");
		exit (1);
	}

	FILE *input, *output;
	input = fopen (argv[1], "r");
	char *outputFilename;
	outputFilename = (char *) malloc (200 * sizeof (char));
	snprintf (outputFilename, 200, "%s.coordination", argv[1]);
	output = fopen (outputFilename, "w");

	int fileStatus = fgetc (input), nAtoms = countNAtoms (nAtoms, input);
	int atomType1 = atoi (argv[2]), atomType2 = atoi (argv[3]);
	float cutoffDistance = atof (argv[4]);

	BOUNDARY simBoundary;
	simBoundary = computeSimBoundary (simBoundary, input);

	DUMP_ENTRIES *atoms;
	atoms = (DUMP_ENTRIES *) malloc (nAtoms * sizeof (DUMP_ENTRIES));

	int coordinationNumber, currentTimestep = 0;

	printf("\n");

	while (fileStatus != EOF)
	{
		atoms = readDump (atoms, nAtoms, input);
		coordinationNumber = findCoordinationNumber (coordinationNumber, atoms, nAtoms, cutoffDistance, atomType1, atomType2, simBoundary);

		fprintf(output, "%d %d\n", currentTimestep, coordinationNumber);
		fprintf(stdout, "Curren timestep: %d                         \r", currentTimestep);

		currentTimestep++;
		fileStatus = fgetc (input);
	}

	fclose (input);
	fclose (output);
	return 0;
}