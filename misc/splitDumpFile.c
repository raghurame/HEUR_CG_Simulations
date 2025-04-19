#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

int countNTimeframes (int nTimeframes_inInput, const char *inputFilename)
{
	char *pipeString;
	pipeString = (char *) malloc (500 * sizeof (char));

	if (strstr (inputFilename, ".gz"))
	{
		snprintf (pipeString, 500, "gzcat %s | grep \"ITEM: TIMESTEP\" | wc -l", inputFilename);
	}
	else
	{
		snprintf (pipeString, 500, "cat %s | grep \"ITEM: TIMESTEP\" | wc -l", inputFilename);
	}

	FILE *countLines;
	countLines = popen (pipeString, "r");

	char *lineString;
	lineString = (char *) malloc (500 * sizeof (char));

	fgets (lineString, 500, countLines);
	sscanf (lineString, "%d\n", &nTimeframes_inInput);

	free (pipeString);

	return nTimeframes_inInput;
}

void splitNextSetOfTimeframes (const char *inputFilename, int nTimeframes_inOutput, int i)
{
	char *pipeString1, *outputFilename, *pipeString2, *pipeString3;
	pipeString1 = (char *) malloc (500 * sizeof (char));
	pipeString2 = (char *) malloc (500 * sizeof (char));
	pipeString3 = (char *) malloc (500 * sizeof (char));
	outputFilename = (char *) malloc (500 * sizeof (char));

	int lastTimeframe = (i * nTimeframes_inOutput) + nTimeframes_inOutput, startTimeframe, lastLine, startLine;

	if (i > 0) {
		startTimeframe = ((i - 1) * nTimeframes_inOutput) + nTimeframes_inOutput; }
	else if (i == 0) {
		startTimeframe = 0; }

	FILE *countLineNumberEnd, *countLineNumberStart, *printOutputfile;

	if (strstr (inputFilename, ".gz")) {
		snprintf (pipeString1, 500, "gzcat %s | grep -n \"ITEM: TIMESTEP\" | head -%d | tail -1", inputFilename, lastTimeframe); }
	else {
		snprintf (pipeString1, 500, "cat %s | grep -n \"ITEM: TIMESTEP\" | head -%d | tail -1", inputFilename, lastTimeframe); }

	if (strstr (inputFilename, ".gz")) {
		snprintf (pipeString2, 500, "gzcat %s | grep -n \"ITEM: TIMESTEP\" | head -%d | tail -1", inputFilename, startTimeframe); }
	else {
		snprintf (pipeString2, 500, "cat %s | grep -n \"ITEM: TIMESTEP\" | head -%d | tail -1", inputFilename, startTimeframe); }

	countLineNumberEnd = popen (pipeString1, "r");
	countLineNumberStart = popen (pipeString2, "r");

	char lineString[5000];

	fgets (lineString, 5000, countLineNumberEnd);
	sscanf (lineString, "%d\n", &lastLine);
	lastLine--;

	fgets (lineString, 5000, countLineNumberStart);
	sscanf (lineString, "%d\n", &startLine);

	startLine++;

	if (i == 0) {
		startLine = 0; }

	snprintf (outputFilename, 500, "%s.split.%d", inputFilename, i);

	if (strstr (inputFilename, ".gz")) {
		snprintf (pipeString3, 500, "gzcat %s | head -%d dump.short | tail -%d >> splitDumpFiles/%s", inputFilename, lastLine, (lastLine - startLine + 2), outputFilename); }
	else {
		snprintf (pipeString3, 500, "cat %s | head -%d dump.short | tail -%d >> splitDumpFiles/%s", inputFilename, lastLine, (lastLine - startLine + 2), outputFilename); }

	// printf("Starting line: %d, End line: %d\n", startLine, lastLine);
	// printf("head -%d dump.short | tail -%d | head -10\n", lastLine, (lastLine - startLine + 2));

	printOutputfile = popen (pipeString3, "w");

	free (pipeString1);
	free (pipeString2);
	free (pipeString3);
	free (outputFilename);
}

int main(int argc, char const *argv[])
{
	if (argc != 3) {
		printf("INCORRECT ARGUMENTS PASSED:\n~~~~~~~~~~~~~~~~~~~~~~~~~~\n\nRequired arguments:\n\n  argv[0] = ./program\n  argv[1] = Input dump file\n  argv[2] = Number of timeframes per each output dump files\n\n");
		exit (1); }

	system ("mkdir splitDumpFiles");

	int nTimeframes_inOutput = atoi (argv[2]);
	int nTimeframes_inInput = countNTimeframes (nTimeframes_inInput, argv[1]);

	printf("Number of timeframes detected in the input file: %d\n", nTimeframes_inInput);

	int nOutputFiles = (int) floor ((float) nTimeframes_inInput/(float) nTimeframes_inOutput);
	printf("Number of output files to be generated... %d\n", nOutputFiles);

	for (int i = 0; i < nOutputFiles; ++i)
	{
		splitNextSetOfTimeframes (argv[1], nTimeframes_inOutput, i);
	}

	return 0;
}