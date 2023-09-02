#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

int main(int argc, char const *argv[])
{
	bool **check1;
	check1 = (bool **) malloc (1000 * sizeof (bool *));

	for (int i = 0; i < 1000; ++i)
	{
		check1[i] = (bool *) malloc (5 * sizeof (bool));
	}

	sleep (10000);

	return 0;
}