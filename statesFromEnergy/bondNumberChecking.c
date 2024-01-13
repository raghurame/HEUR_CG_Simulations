#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int main(int argc, char const *argv[])
{
	int bondNumber1, bondNumber2, bondNumber3;

	for (int i = 1; i < 1300; ++i)
	{
		bondNumber1 = (i - (i / (80)) - 1) / 2;
		bondNumber2 = (i - (i / (80 + 1)) - 1) / 2;

		if (i % 81 == 0)
		{
			printf("==> partID ==> ");
		}

		printf("i: %d ;", i);
		printf("bond number 1: %d; ", bondNumber1);
		printf("bond number 2: %d; \n", bondNumber2);
		usleep (10000);
	}

	return 0;
}