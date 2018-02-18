#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>

#define TRYS 5000000

static int throw() {
	double x, y;
	x = (double) rand() / (double) RAND_MAX;
	y = (double) rand() / (double) RAND_MAX;
	if ((x * x + y * y) <= 1.0)
		return 1;

	return 0;
}

int main(int argc, char **argv) {
	int globalCount = 0, globalSamples = TRYS;

	printf("How many threads? : ");

	int nthreads=0;

	// Max length of 99 characters; puts a null terminated string in path, thus 99 chars + null is the max
	scanf( "%d" , &nthreads);

#pragma omp parallel for reduction ( +:globalCount)  num_threads(nthreads)
	for (int i = 0; i < globalSamples; ++i) {
		globalCount += throw();

		// find end of each thread and use it to print last globalCount-value
		if (i
				== ((omp_get_thread_num() + 1) * globalSamples
						/ omp_get_num_threads()) - 1) {
			printf("Thread %d: treffer %d \n", omp_get_thread_num(),
					globalCount);
		}

	}

	double pi = 4.0 * (double) globalCount / (double) (globalSamples);

	printf("pi is %.9lf\n", pi);

	return 0;
}
