#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>

int main(int argc, char **argv) {

#pragma omp parallel
	{
		printf("Hello World from %d thread of %d \n", omp_get_thread_num(),
				omp_get_max_threads());
	}

  return 0;
}
