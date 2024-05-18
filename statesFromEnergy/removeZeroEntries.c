#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

typedef struct dumpEnergy
{
	int atom1, atom2;
	float distance, energy;
} DUMP_ENERGY;

int countDumpEntries (FILE *inputDump, int nDumpEntries)
{
	char lineString[5000];
	nDumpEntries = 0;

	for (int i = 0; i < 4; ++i)
	{
		fgets (lineString, 3000, inputDump);
	}

	sscanf (lineString, "%d\n", &nDumpEntries);

	for (int i = 0; i < 5; ++i)
	{
		fgets (lineString, 3000, inputDump);
	}

	return nDumpEntries;
}

DUMP_ENERGY *saveAllEnergyEntries (DUMP_ENERGY *energyEntries, FILE*inputDump, int nDumpEntries)
{
	char lineString[5000];

	for (int i = 0; i < nDumpEntries; ++i)
	{
		fgets (lineString, 1000, inputDump);
		sscanf (lineString, "%d %d %f %f\n", &energyEntries[i].atom1, &energyEntries[i].atom2, &energyEntries[i].distance, &energyEntries[i].energy);
	}

	return energyEntries;
}

int countNonZeroEntries (int nNonZeroEntries, DUMP_ENERGY *energyEntries, int nDumpEntries)
{
	nNonZeroEntries = 0;

	for (int i = 0; i < nDumpEntries; ++i)
	{
		if (energyEntries[i].energy < -1)
		{
			nNonZeroEntries++;

/*			if (i > 0)
			{
				if (energyEntries[i].atom1 == energyEntries[i - 1].atom1 && energyEntries[i].atom2 == energyEntries[i - 1].atom2)
				{
					printf("%d %d %f %f\n", energyEntries[i - 1].atom1, energyEntries[i - 1].atom2, energyEntries[i - 1].distance, energyEntries[i - 1].energy);
					printf("%d %d %f %f\n\n", energyEntries[i].atom1, energyEntries[i].atom2, energyEntries[i].distance, energyEntries[i].energy);
					usleep (100000);
				}
			}
*/		}
	}

	return nNonZeroEntries;
}

void printNonZeroEntries (FILE *outputDump, int nNonZeroEntries, DUMP_ENERGY *energyEntries, int nDumpEntries, int currentTimeframe)
{
	fprintf(outputDump, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ENTRIES\n%d\nITEM: BOX BOUNDS xy xz yz pp pp pp\n0.0 0.0 0.0\n0.0 0.0 0.0\n0.0 0.0 0.0\nITEM: ENTRIES c_1[1] c_1[2] c_2[1] c_2[2]\n", currentTimeframe + 1, nNonZeroEntries);

	for (int i = 0; i < nDumpEntries; ++i)
	{
		if (energyEntries[i].energy < -1)
		{
			if (i > 0)
			{
				if (energyEntries[i - 1].atom1 == energyEntries[i].atom1 && energyEntries[i - 1].atom2 == energyEntries[i].atom2)
				{
					// printf("[%d] ignored -> %d %d %f %f\n", currentTimeframe + 1, energyEntries[i].atom1, energyEntries[i].atom2, energyEntries[i].distance, energyEntries[i].energy);
					// usleep (100000);
				}
				else
				{
					fprintf(outputDump, "%d %d %f %f\n", energyEntries[i].atom1, energyEntries[i].atom2, energyEntries[i].distance, energyEntries[i].energy);
					fflush (outputDump);
					// printf("%d %d %f %f\n", energyEntries[i].atom1, energyEntries[i].atom2, energyEntries[i].distance, energyEntries[i].energy);
				}
			}
			else
			{
				fprintf(outputDump, "%d %d %f %f\n", energyEntries[i].atom1, energyEntries[i].atom2, energyEntries[i].distance, energyEntries[i].energy);
				fflush (outputDump);
				// printf("%d %d %f %f\n", energyEntries[i].atom1, energyEntries[i].atom2, energyEntries[i].distance, energyEntries[i].energy);
			}
		}
	}
}

int main(int argc, char const *argv[])
{
	FILE *inputDump, *outputDump;
	char *pipeString, *outputFileName;
	pipeString = (char *) malloc (500 * sizeof (char));
	snprintf (pipeString, 500, "zcat %s", argv[1]);
	outputFileName = (char *) malloc (500 * sizeof (char));
	snprintf (outputFileName, 500, "%s.mini", argv[1]);

	outputDump = fopen (outputFileName, "w");

	if (strstr (argv[1], ".gz"))
	{
		inputDump = popen (pipeString, "r");
	}
	else
	{
		inputDump = fopen (argv[1], "r");
	}

	int file_status = fgetc (inputDump);
	int nDumpEntries = 0, nNonZeroEntries = 0;
	int currentTimeframe = 0;

	while (file_status > 0)
	{
		printf("Scanning timeframe: %d                  \r", currentTimeframe + 1);
		fflush (stdout);

		DUMP_ENERGY *energyEntries;
		nDumpEntries = countDumpEntries (inputDump, nDumpEntries);		
		energyEntries = (DUMP_ENERGY *) malloc (nDumpEntries * sizeof (DUMP_ENERGY));
		energyEntries = saveAllEnergyEntries (energyEntries, inputDump, nDumpEntries);
		nNonZeroEntries = countNonZeroEntries (nNonZeroEntries, energyEntries, nDumpEntries);
		printNonZeroEntries (outputDump, nNonZeroEntries, energyEntries, nDumpEntries, currentTimeframe);

		currentTimeframe++;
		file_status = fgetc (inputDump);

		free (energyEntries);
	}

	char *compressString;
	compressString = (char *) malloc (500 * sizeof (char));
	snprintf (compressString, 500, "gzip %s", outputFileName);
	system (compressString);

	return 0;
}

// mv /scratch/rlarson_root/rlarson0/relanche/newDirMar16/HEUR/sriramRerun/p100/scaling/ghostParticleSize/akshatConcentration2/r0.64conc20_eps10/stats5.0588/wi8e-5/dump.stat.gz .