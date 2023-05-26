#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>

typedef struct dumpInformation
{
	int id, type;
	float x, y, z;
} DUMP_INFORMATION;

typedef struct coordinates
{
	float x, y, z;
} COORDINATES;

typedef struct boxDimension
{
	float xlo, xhi, xtilt, ylo, yhi, ytilt, zlo, zhi, ztilt;
} BOX_DIMENSION;

int countNAtoms (int nAtoms, const char filename[])
{
	FILE *pipeNAtoms;
	char *pipeNAtoms_string;
	pipeNAtoms_string = (char *) malloc (500 * sizeof (char));
	snprintf (pipeNAtoms_string, 500, "head -4 %s | tail -1", filename);
	pipeNAtoms = popen (pipeNAtoms_string, "r");
	char lineString[2000];
	fgets (lineString, 2000, pipeNAtoms);
	sscanf (lineString, "%d\n", &nAtoms);
	fclose (pipeNAtoms);
	return nAtoms;
}

int main(int argc, char const *argv[])
{
	FILE *inputDumpfile, *maxMolBoundary;
	inputDumpfile = fopen (argv[1], "r");
	maxMolBoundary = fopen (argv[2], "w");

	DUMP_INFORMATION *inputDumpCoords;
	int nAtoms = countNAtoms (nAtoms, argv[1]), fileStatus = fgetc (inputDumpfile);
	inputDumpCoords = (DUMP_INFORMATION *) malloc (nAtoms * sizeof (DUMP_INFORMATION));
	char lineString[2000];

	COORDINATES maxDimension, currentDimension;
	maxDimension.x = 0;
	maxDimension.y = 0;
	maxDimension.z = 0;
	int currentAtom = 0, currentTimestep = 0;

	BOX_DIMENSION simDimension;

	while (fileStatus != EOF)
	{
		currentTimestep++;
		for (int i = 0; i < 9; ++i) {
			currentAtom = 0;
			fgets (lineString, 2000, inputDumpfile);

			if (i == 6) {
				sscanf (lineString, "%f %f %f\n", &simDimension.xlo, &simDimension.xhi, &simDimension.xtilt); }
			if (i == 7) {
				sscanf (lineString, "%f %f %f\n", &simDimension.ylo, &simDimension.yhi, &simDimension.ytilt); }
			if (i == 8) {
				sscanf (lineString, "%f %f %f\n", &simDimension.zlo, &simDimension.zhi, &simDimension.ztilt); }
		}

		for (int i = 0; i < nAtoms; ++i)
		{
			fgets (lineString, 2000, inputDumpfile);
			sscanf (lineString, "%d %d %f %f %f\n", &inputDumpCoords[i].id, &inputDumpCoords[i].type, &inputDumpCoords[i].x, &inputDumpCoords[i].y, &inputDumpCoords[i].z);
		}

		while ((currentAtom + 1) < nAtoms)
		{
			if (inputDumpCoords[currentAtom].type == inputDumpCoords[currentAtom + 1].type)
			{
				currentDimension.x = inputDumpCoords[currentAtom].x - inputDumpCoords[currentAtom + 1].x;
				currentDimension.y = inputDumpCoords[currentAtom].y - inputDumpCoords[currentAtom + 1].y;
				currentDimension.z = inputDumpCoords[currentAtom].z - inputDumpCoords[currentAtom + 1].z;

				if (currentDimension.x > maxDimension.x) {
					maxDimension.x = currentDimension.x; }

				if (currentDimension.y > maxDimension.y) {
					maxDimension.y = currentDimension.y; }

				if (currentDimension.z > maxDimension.z) {
					maxDimension.z = currentDimension.z; }

				currentAtom += 2;
			}
			else
			{
				currentAtom++;
			}
		}

		fprintf(stdout, "  x: %f; y: %f; z: %f (current timestep: %d)                               \r", maxDimension.x, maxDimension.y, maxDimension.z, currentTimestep);
		fflush (stdout);
	}

	printf("Max dimensions\n\n  x: %f\n  y: %f\n  z: %f\n\n", maxDimension.x, maxDimension.y, maxDimension.z);
	fprintf(maxMolBoundary, "Max dimensions\n\n  x: %f\n  y: %f\n  z: %f", maxDimension.x, maxDimension.y, maxDimension.z);

	fclose (inputDumpfile);
	return 0;
}