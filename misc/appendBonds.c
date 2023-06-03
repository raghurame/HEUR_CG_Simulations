#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <stdbool.h>

typedef struct dataFileInfo
{
	int sino, atomType, molType;
	float x, y, z;
} DATAFILE;

int main(int argc, char const *argv[])
{
	if (argc == 1)
	{
		printf("REQURIED ARGUMENTS:\n~~~~~~~~~~~~~~~~~~~\n\n{~} argv[0] = program\n{~} argv[1] = Input data filename\n\n");
		exit (1);
	}

	FILE *inputDatafile, *writeDatafile;
	inputDatafile = fopen (argv[1], "r");
	char *outputFilename;
	outputFilename = (char *) malloc (200 * sizeof (char));
	snprintf (outputFilename, 200, "%s.withBonds", argv[1]);
	writeDatafile = fopen (outputFilename, "w");

	int fileStatus = fgetc (inputDatafile), nAtoms, atomID;
	char lineString[2000];
	bool atomsSection = false, nAtomsLine = false;
	DATAFILE *inputDataInformation;

	while (fgets (lineString, 2000, inputDatafile) != NULL)
	{
		if (strstr (lineString, "atoms")) {
			nAtomsLine = true; }

		if (strstr (lineString, "Atoms")) {
			atomsSection = true;
			fgets (lineString, 2000, inputDatafile);
			fgets (lineString, 2000, inputDatafile); }

		if (nAtomsLine)
		{
			sscanf (lineString, "%d\n", &nAtoms);
			nAtomsLine = false;
			inputDataInformation = (DATAFILE *) malloc (nAtoms * sizeof (DATAFILE));
		}

		if (atomsSection)
		{
			sscanf (lineString, "%d", &atomID);
			if (atomID <= nAtoms) {
				sscanf (lineString, "%d %d %d %f %f %f\n", &inputDataInformation[atomID - 1].sino, &inputDataInformation[atomID - 1].atomType, &inputDataInformation[atomID - 1].molType, &inputDataInformation[atomID - 1].x, &inputDataInformation[atomID - 1].y, &inputDataInformation[atomID - 1].z); }
		}
	}

	int currentAtom1 = 1, currentAtom2 = 2, bondID = 0;

	fprintf(writeDatafile, "\n\nBonds\n\n");

	while (currentAtom1 <= nAtoms && currentAtom2 <= nAtoms)
	{
		if (inputDataInformation[currentAtom1 - 1].molType == inputDataInformation[currentAtom2 - 1].molType) {
			bondID++;
			fprintf(writeDatafile, "%d %d %d %d\n", bondID, 1, currentAtom1, currentAtom2);
			currentAtom1 += 2;
			currentAtom2 += 2;
		}
		else
		{
			currentAtom1++;
			currentAtom2++;
		}

	}

	fclose (inputDatafile);
	fclose (writeDatafile);
	return 0;
}