#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int main(int argc, char const *argv[])
{
	int bondNumber1, bondNumber2, bondNumber3;
	int i = atoi (argv[1]);

	for (int i = 1; i < 1300; ++i)
	{
		bondNumber1 = (i - (i / 80) - 1) / 2;
		bondNumber2 = ceil ((float)(i - (int)floor (i / 80)) / 2) - 1;
		bondNumber3 = floor ((i - (int)floor (i / (80 + 1)) - 1) / 2);

		if (i % 81 == 0)
		{
			printf("==> partID ==> ");
		}

		printf("i: %d ;", i);
		printf("bond number 1: %d; ", bondNumber1);
		printf("bond number 2: %d; ", bondNumber2);
		printf("bond number 3: %d\n", bondNumber3);
		usleep (100000);
	}

	int partID = 0;
	// printing particle IDs
	for (int i = 0; i < 100; ++i)
	{
		partID += (80 + 1);
		printf("%d\n", partID);
		usleep (100000);
	}

	return 0;
}