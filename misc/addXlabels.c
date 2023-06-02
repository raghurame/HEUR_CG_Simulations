#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

int countNLines (int nLines, FILE *input)
{
	char lineString[2000];
	nLines = 0;

	while (fgets (lineString, 2000, input) != NULL)
	{
		nLines++;
	}

	rewind (input);
	return nLines;
}

float *readYLabels (float *yLabels, int nLines, FILE *input)
{
	char lineString[2000];
	int currentLine = 0;

	while (fgets (lineString, 2000, input) != NULL)
	{
		if (currentLine < nLines)
		{
			sscanf (lineString, "%*f %f\n", &yLabels[currentLine]);
			currentLine++;
		}
	}

	return yLabels;
}

float *createXLabels (float *xLabels, int nLines, float dt, float printFrequency)
{
	for (int i = 0; i < nLines; ++i)
	{
		xLabels[i] = (float)i * dt * printFrequency;
	}

	return xLabels;
}

void writeNewDatFile (float *xLabels, float *yLabels, int nLines, FILE *output)
{
	for (int i = 0; i < nLines; ++i)
	{
		fprintf(output, "%f %f\n", xLabels[i], yLabels[i]);
	}
}

int main(int argc, char const *argv[])
{
	FILE *input, *output;
	char *outputFilename;
	outputFilename = (char *) malloc (200 * sizeof (char));
	input = fopen (argv[1], "r");
	snprintf (outputFilename, 200, "%s.xLabelsAdded", argv[1]);
	output = fopen (outputFilename, "w");

	float dt = atof (argv[2]), printFrequency = atof (argv[3]);
	int nLines = countNLines (nLines, input);
	float *xLabels, *yLabels;
	yLabels = (float *) malloc (nLines * sizeof (float));
	xLabels = (float *) malloc (nLines * sizeof (float));

	yLabels = readYLabels (yLabels, nLines, input);
	xLabels = createXLabels (xLabels, nLines, dt, printFrequency);

	writeNewDatFile (xLabels, yLabels, nLines, output);

	fclose (input);
	fclose (output);
	return 0;
}