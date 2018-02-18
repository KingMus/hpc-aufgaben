#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>

int main(int argc, char **argv) {

#pragma omp parallel sections num_threads(4)
{

	#pragma omp section
	{
	printf("Hallo Welt from %d thread of %d \n", omp_get_thread_num(),
			omp_get_max_threads());
	printf("Hej varlden from %d thread of %d \n", omp_get_thread_num(),
			omp_get_max_threads());
	}
	#pragma omp section
	{
	printf("Bonjour tout le monde from %d thread of %d \n",
			omp_get_thread_num(), omp_get_max_threads());
	}
	#pragma omp section
	{
		printf("Hola mundo from %d thread of %d \n", omp_get_thread_num(),	omp_get_max_threads());
	}
	#pragma omp section
	{
	printf("Hello World from %d thread of %d \n", omp_get_thread_num(),
			omp_get_max_threads());
	}
}

  return 0;
}
