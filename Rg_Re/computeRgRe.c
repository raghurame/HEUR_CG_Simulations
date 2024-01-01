#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

int countNTimeframes (int nTimeframes, const char *filename)
{
	FILE *input;

	char *pipeString;
	pipeString = (char *) malloc (500 * sizeof (char));
	snprintf (pipeString, 500, "zcat %s", filename);

	if (strstr (filename, "gz"))
	{
		input = popen (pipeString, "r");
	}
	else
	{
		input = fopen (filename, "r");
	}

	char lineString[3000];
	nTimeframes = 0;

	while (fgets (lineString, 3000, input) != NULL)
	{
		if (strstr (lineString, "ITEM: TIMESTEP"))
		{
			nTimeframes++;
		}
	}

	printf("Found %d timeframes...\n", nTimeframes);
	fclose (input);
	return nTimeframes;
}

int countNBonds (int nBonds, const char *filename)
{
	FILE *input;

	char *pipeString;
	pipeString = (char *) malloc (500 * sizeof (char));
	snprintf (pipeString, 500, "zcat %s", filename);

	if (strstr (filename, "gz"))
	{
		input = popen (pipeString, "r");
	}
	else
	{
		input = fopen (filename, "r");
	}

	char lineString[3000];

	for (int i = 0; i < 4; ++i)
	{
		fgets (lineString, 3000, input);
	}

	sscanf (lineString, "%d\n", &nBonds);

	printf("Found %d bonds in the input file...\n", nBonds);

	fclose (input);
	return nBonds;
}

typedef struct bondDistHist
{
	float min, max;
	int count;
} BOND_DISTANCE_HISTOGRAM;

typedef struct bins
{
	int nBins;
	float minBondDistance, maxBondDistance, binDistance;
} BINS;

float *storeBondDistances (float *bondDistances, FILE *input, int nBonds)
{
	char lineString[3000];

	for (int i = 0; i < 9; ++i)
	{
		fgets (lineString, 3000, input);
	}

	for (int i = 0; i < nBonds; ++i)
	{
		fgets (lineString, 3000, input);
		sscanf (lineString, "%*d %*d %*f %f\n", &bondDistances[i]);
	}

	return bondDistances;
}

BOND_DISTANCE_HISTOGRAM *initHistogram (BOND_DISTANCE_HISTOGRAM *bondDistancesHistogram, BINS binInfo)
{
	bondDistancesHistogram[0].min = binInfo.minBondDistance;
	bondDistancesHistogram[0].max = binInfo.minBondDistance + binInfo.binDistance;
	bondDistancesHistogram[0].count = 0;

	for (int i = 1; i < binInfo.nBins; ++i)
	{
		bondDistancesHistogram[i].min = bondDistancesHistogram[i - 1].max;
		bondDistancesHistogram[i].max = bondDistancesHistogram[i].min + binInfo.binDistance;
		bondDistancesHistogram[i].count = 0;
	}

	return bondDistancesHistogram;
}

BOND_DISTANCE_HISTOGRAM *computeBondDistanceHistogram (BOND_DISTANCE_HISTOGRAM *bondDistancesHistogram, float *bondDistances, BINS binInfo, int nBonds)
{
	for (int i = 0; i < binInfo.nBins; ++i)
	{
		for (int j = 0; j < nBonds; ++j)
		{
			if ((bondDistances[j] > bondDistancesHistogram[i].min) && (bondDistances[j] <= bondDistancesHistogram[i].max))
			{
				bondDistancesHistogram[i].count++;
			}
		}
	}

	return bondDistancesHistogram;
}

int findMaxHistoCount (int maxHistoCount, BOND_DISTANCE_HISTOGRAM *bondDistancesHistogram, BINS binInfo)
{
	for (int i = 0; i < binInfo.nBins; ++i)
	{
		if (bondDistancesHistogram[i].count > maxHistoCount)
		{
			maxHistoCount = bondDistancesHistogram[i].count;
		}
	}

	return maxHistoCount;
}

void printNormalizedHistogram (int maxHistoCount, BINS binInfo, BOND_DISTANCE_HISTOGRAM *bondDistancesHistogram)
{
	FILE *histoOutput;
	histoOutput = fopen ("bondDistances.histogram", "w");

	for (int i = 0; i < binInfo.nBins; ++i)
	{
		fprintf(histoOutput, "%f\n", (float)bondDistancesHistogram[i].count / (float)maxHistoCount);
	}

	fclose (histoOutput);
}

float computeAverageBondDistance (float averageBondDistance, float *bondDistances, int nBonds)
{
	for (int i = 0; i < nBonds; ++i)
	{
		averageBondDistance += bondDistances[i];
	}

	return averageBondDistance;
}

float computeStandardDeviation (float deviationBondDistance, float averageBondDistance, float *bondDistances, int nBonds)
{
	for (int i = 0; i < nBonds; ++i)
	{
		deviationBondDistance += (bondDistances[i] - averageBondDistance) * (bondDistances[i] - averageBondDistance);
	}

	return deviationBondDistance;
}

int main(int argc, char const *argv[])
{
	FILE *input;

	char *pipeString;
	pipeString = (char *) malloc (500 * sizeof (char));
	snprintf (pipeString, 500, "zcat %s", argv[1]);

	if (strstr (argv[1], "gz"))
	{
		input = popen (pipeString, "r");
	}
	else
	{
		input = fopen (argv[1], "r");
	}

	int nTimeframes = countNTimeframes (nTimeframes, argv[1]), nBonds = countNBonds (nBonds, argv[1]);

	float *bondDistances;
	bondDistances = (float **) malloc (nBonds * sizeof (float *));

	BINS binInfo;
	binInfo.nBins = atoi (argv[2]);
	binInfo.minBondDistance = atof (argv[3]);
	binInfo.maxBondDistance = atof (argv[4]);
	binInfo.binDistance = (binInfo.maxBondDistance - binInfo.minBondDistance) / binInfo.nBins;

	BOND_DISTANCE_HISTOGRAM *bondDistancesHistogram;
	bondDistancesHistogram = (BOND_DISTANCE_HISTOGRAM *) malloc (binInfo.nBins * sizeof (BOND_DISTANCE_HISTOGRAM));

	bondDistancesHistogram = initHistogram (bondDistancesHistogram, binInfo);

	int file_status, currentTimeframe = 0;

	if (strstr (argv[1], "gz"))
	{
		pclose (input);
		input = popen (pipeString, "r");
	}
	else
	{
		rewind (input);
	}

	file_status = fgetc (input);

	float averageBondDistance = 0;

	while (file_status > 0)
	{
		bondDistances = storeBondDistances (bondDistances, input, nBonds);

		bondDistancesHistogram = computeBondDistanceHistogram (bondDistancesHistogram, bondDistances, binInfo, nBonds);

		averageBondDistance = computeAverageBondDistance (averageBondDistance, bondDistances, nBonds);

		currentTimeframe++;
		printf("Reading timeframe %d/%d...                 \r", currentTimeframe, nTimeframes);
		fflush (stdout);

		file_status = fgetc (input);
	}

	averageBondDistance = averageBondDistance / (nBonds * currentTimeframe);

	if (strstr (argv[1], "gz"))
	{
		pclose (input);
		input = popen (pipeString, "r");
	}
	else
	{
		rewind (input);
	}

	file_status = fgetc (input);
	currentTimeframe = 0;

	printf("\n");
	float deviationBondDistance = 0;

	while (file_status > 0)
	{
		bondDistances = storeBondDistances (bondDistances, input, nBonds);

		deviationBondDistance = computeStandardDeviation (deviationBondDistance, averageBondDistance, bondDistances, nBonds);

		currentTimeframe++;
		printf("Calculating standard deviation %d/%d...                 \r", currentTimeframe, nTimeframes);
		fflush (stdout);

		file_status = fgetc (input);
	}

	printf("\n\n");

	deviationBondDistance /= (nBonds * currentTimeframe);
	deviationBondDistance = sqrt (deviationBondDistance);

	printf("Average bond distance: %f\n", averageBondDistance);
	printf("Standard deviation: %f\n", deviationBondDistance);

	int maxHistoCount = 0;
	maxHistoCount = findMaxHistoCount (maxHistoCount, bondDistancesHistogram, binInfo);
	printNormalizedHistogram (maxHistoCount, binInfo, bondDistancesHistogram);

	if (strstr (argv[1], "gz"))
	{
		pclose (input);
	}
	else
	{
		fclose (input);
	}
	return 0;
}

// 9V^fRZGYn4A7Lc4FA4ATfzZBdHtkorbDw4%69oA$3Hhtb