#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

typedef struct distribution
{
	float minBinValue, maxBinValue;
	int count;
} DIST;

int countNLines (const char *inputFilename, int nLines)
{
	FILE *pipeFile;
	char *pipeString, lineString[1000];
	pipeString = (char *) malloc (500 * sizeof (char));
	snprintf (pipeString, 500, "wc -l %s", inputFilename);

	pipeFile = popen (pipeString, "r");

	fgets (lineString, 1000, pipeFile);
	sscanf (lineString, "%d\n", &nLines);
	printf("Number of lines in the input file: %d\n", nLines);

	pclose (pipeFile);

	return nLines;
}

int findMinValue (FILE *inputFile, int minValue, int nLines)
{
	char lineString[1000];
	rewind (inputFile);

	minValue = 0;
	int tempValue = 0;

	while (fgets (lineString, 1000, inputFile) != NULL)
	{
		sscanf (lineString, "%d\n", &tempValue);

		if (tempValue > 0)
		{
			if (minValue == 0) {
				minValue = tempValue; }

			if (tempValue < minValue) {
				minValue = tempValue; }
		}
	}

	printf("Min value: %d\n", minValue);

	return minValue;
}

int findMaxValue (FILE *input, int maxValue, int nLines)
{
	char lineString[1000];
	rewind (input);

	maxValue = 0;
	int tempValue = 0;

	while (fgets (lineString, 1000, input) != NULL)
	{
		sscanf (lineString, "%d\n", &tempValue);

		if (tempValue > 0)
		{
			if (maxValue == 0) {
				maxValue = tempValue; }

			if (tempValue > maxValue) {
				maxValue = tempValue; }
		}
	}

	printf("Max value: %d\n", maxValue);

	return maxValue;
}

DIST *setBinValues (DIST *bins, int nBins, float binWidth, int minValue)
{
	for (int i = 0; i < nBins; ++i)
	{
		if (i == 0)
		{
			bins[i].minBinValue = minValue;
			bins[i].maxBinValue = minValue + binWidth;
			bins[i].count = 0;
		}
		else
		{
			bins[i].minBinValue = bins[i - 1].maxBinValue;
			bins[i].maxBinValue = bins[i - 1].maxBinValue + binWidth;
			bins[i].count = 0;
		}
	}

	return bins;
}

void checkBinValues (DIST *bins, int nBins)
{
	for (int i = 0; i < nBins; ++i)
	{
		printf("%f %f => %d\n", bins[i].minBinValue, bins[i].maxBinValue, bins[i].count);
	}
}

DIST *computeDistribution (DIST *bins, int nBins, FILE *input)
{
	char lineString[1000];
	rewind (input);
	int tempValue = 0;

	while (fgets (lineString, 1000, input) != NULL)
	{
		sscanf (lineString, "%d\n", &tempValue);

		if (tempValue > 0)
		{
			for (int i = 0; i < nBins; ++i)
			{
				if ((float)tempValue > bins[i].minBinValue && (float)tempValue <= bins[i].maxBinValue)
				{
					bins[i].count++;
					break;
				}
			}
		}
	}

	return bins;
}

void printDistribution (DIST *bins, int nBins, FILE *outputFile)
{
	for (int i = 0; i < nBins; ++i)
	{
		fprintf(outputFile, "%f %f %d\n", bins[i].minBinValue, bins[i].maxBinValue, bins[i].count);
	}
}

int main(int argc, char const *argv[])
{
	if (argc != 4)
	{
		printf("INCORRECT ARGUMENTS PASSED:\n~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n {~} argv[0] = ./program\n {~} argv[1] = Input filename\n {~} argv[2] = Number of bins in the output distribution\n {~} argv[3] = Output filename\n\n");
		exit (1);
	}

	FILE *inputFile, *outputFile;

	inputFile = fopen (argv[1], "r");
	int nBins = atoi (argv[2]);
	outputFile = fopen (argv[3], "w");

	int nLines = 0;
	nLines = countNLines (argv[1], nLines);
	// printf("Number of lines in the input file: %d\n", nLines);

	int maxValue = 0, minValue = 0;  

	maxValue = findMaxValue (inputFile, maxValue, nLines);
	minValue = findMinValue (inputFile, minValue, nLines);

	// printf("Min value: %d\n", minValue);
	// printf("Max value: %d\n", maxValue);

	float rangeOfValue = 0;
	rangeOfValue = (float)(maxValue - minValue);
	float binWidth = 0;
	binWidth = (float)(rangeOfValue / nBins);

	DIST *bins;
	bins = (DIST *) malloc (nBins * sizeof (DIST));

	bins = setBinValues (bins, nBins, binWidth, minValue);
	// checkBinValues (bins, nBins);

	bins = computeDistribution (bins, nBins, inputFile);
	// checkBinValues (bins, nBins);
	printDistribution (bins, nBins, outputFile);

	fclose (inputFile);
	fclose (outputFile);
	return 0;
}