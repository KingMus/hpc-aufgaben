#include <endian.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <sys/time.h>

#define calcIndex(width, x,y)  ((y)*(width) + (x))

long TimeSteps = 100;

void writeVTK2(long timestep, double *data, char prefix[1024], int wstart,
		int hstart, long w, long h, char threadnum[1024]) {
	char filename[2048];
	int x, y;

	long offsetX = 0;
	long offsetY = 0;
	float deltax = 1.0;
	float deltay = 1.0;
	long nxy = (w - wstart) * (h - hstart) * sizeof(float);

	snprintf(filename, sizeof(filename), "%s-%s%05ld%s", prefix, threadnum,
			timestep, ".vti");
	FILE* fp = fopen(filename, "w");

	fprintf(fp, "<?xml version=\"1.0\"?>\n");
	fprintf(fp,
			"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(fp,
			"<ImageData WholeExtent=\"%d %d %d %d %d %d\" Origin=\"0 0 0\" Spacing=\"%le %le %le\">\n",
			offsetX, offsetX + (w - wstart), offsetY, offsetY + (h - hstart), 0,
			0, deltax, deltax, 0.0);
	fprintf(fp, "<CellData Scalars=\"%s\">\n", prefix);
	fprintf(fp,
			"<DataArray type=\"Float32\" Name=\"%s\" format=\"appended\" offset=\"0\"/>\n",
			prefix);
	fprintf(fp, "</CellData>\n");
	fprintf(fp, "</ImageData>\n");
	fprintf(fp, "<AppendedData encoding=\"raw\">\n");
	fprintf(fp, "_");
	fwrite((unsigned char*) &nxy, sizeof(long), 1, fp);

	for (y = hstart; y < h; y++) {
		for (x = wstart; x < w; x++) {

			long max = w;

			if (h > max) {
				max = h;
			}

			if(threadnum == "_t1_"){
				max = h*2;
			}

			float value = data[calcIndex(max, x, y)];

			fwrite((unsigned char*) &value, sizeof(float), 1, fp);
		}
	}

	fprintf(fp, "\n</AppendedData>\n");
	fprintf(fp, "</VTKFile>\n");
	fclose(fp);
}

int countLifingsPeriodic(double* currentfield, int x, int y, int w, int h) {
	int neighbours = 0;
	for (int y1 = y - 1; y1 <= y + 1; y1++) {
		for (int x1 = x - 1; x1 <= x + 1; x1++) {
			if (currentfield[calcIndex(w, (x1 + w) % w, (y1 + h) % h)]) {
				neighbours++;
			}
		}
	}
	return neighbours;

}

void show(double* currentfield, int w, int h) {
	printf("\033[H");
	int x, y;
	for (y = 0; y < h; y++) {
		for (x = 0; x < w; x++)
			printf(
					currentfield[calcIndex(w, x, y)] ?
							"\033[07m  \033[m" : "  ");
		//printf("\033[E");
		printf("\n");
	}
	fflush(stdout);
}

void evolve(double* currentfield, double* newfield, int w, int h, long t) {

#pragma omp parallel sections num_threads(4)
	{

		//Thread 1
#pragma omp section
		{

			int x, y;

			int xstart, ystart, xende, yende;

			xstart = 0;
			xende = w / 2;
			ystart = 0;
			yende = h / 2;

			for (y = ystart; y < yende; y++) {
				for (x = xstart; x < xende; x++) {

					int n = countLifingsPeriodic(currentfield, x, y, w, h);
					if (currentfield[calcIndex(w, x, y)])
						n--;

					newfield[calcIndex(w, x, y)] = (n == 3
							|| (n == 2 && currentfield[calcIndex(w, x, y)]));

				}
			}

			writeVTK2(t, currentfield, "gol", xstart, ystart, xende, yende,
					"_t1_");
			printf("Thread %d hat berechnet \n", omp_get_thread_num());
		}

		//Thread 2
#pragma omp section
		{

			int x, y;

			int xstart, ystart, xende, yende;

			xstart = 0;
			xende = w / 2;
			ystart = h / 2;
			yende = h;

			for (y = ystart; y < yende; y++) {
				for (x = xstart; x < xende; x++) {

					int n = countLifingsPeriodic(currentfield, x, y, w, h);
					if (currentfield[calcIndex(w, x, y)])
						n--;

					newfield[calcIndex(w, x, y)] = (n == 3
							|| (n == 2 && currentfield[calcIndex(w, x, y)]));

				}
			}

			writeVTK2(t, currentfield, "gol", xstart, ystart, xende, yende,
					"_t2_");
			printf("Thread %d hat berechnet \n", omp_get_thread_num());

		}

		//Thread 3
#pragma omp section
		{
			int x, y;

			int xstart, ystart, xende, yende;

			xstart = w / 2;
			xende = w;
			ystart = 0;
			yende = h / 2;

			for (y = ystart; y < yende; y++) {
				for (x = xstart; x < xende; x++) {

					int n = countLifingsPeriodic(currentfield, x, y, w, h);
					if (currentfield[calcIndex(w, x, y)])
						n--;

					newfield[calcIndex(w, x, y)] = (n == 3
							|| (n == 2 && currentfield[calcIndex(w, x, y)]));

				}
			}

			writeVTK2(t, currentfield, "gol", xstart, ystart, xende, yende,
					"_t3_");
			printf("Thread %d hat berechnet \n", omp_get_thread_num());
		}

		//Thread 4
#pragma omp section
		{
			int x, y;

			int xstart, ystart, xende, yende;

			xstart = w / 2;
			xende = w;
			ystart = h / 2;
			yende = h;

			for (y = ystart; y < yende; y++) {
				for (x = xstart; x < xende; x++) {

					int n = countLifingsPeriodic(currentfield, x, y, w, h);
					if (currentfield[calcIndex(w, x, y)])
						n--;

					newfield[calcIndex(w, x, y)] = (n == 3
							|| (n == 2 && currentfield[calcIndex(w, x, y)]));

				}
			}

			writeVTK2(t, currentfield, "gol", xstart, ystart, xende, yende,
					"_t4_");
			printf("Thread %d hat berechnet \n", omp_get_thread_num());
		}

	}

}

void filling(double* currentfield, int w, int h) {
	int i;
	for (i = 0; i < h * w; i++) {
		currentfield[i] = (rand() < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
	}
}

void game(int w, int h) {
	double *currentfield = calloc(w * h, sizeof(double));
	double *newfield = calloc(w * h, sizeof(double));

//printf("size unsigned %d, size long %d\n",sizeof(float), sizeof(long));

	filling(currentfield, w, h);
	long t;
	for (t = 0; t < TimeSteps; t++) {
		show(currentfield, w, h);

		evolve(currentfield, newfield, w, h, t);

		printf("%ld timestep\n", t);

		usleep(200000);

		//SWAP
		double *temp = currentfield;
		currentfield = newfield;
		newfield = temp;
	}

	free(currentfield);
	free(newfield);

}

int main(int c, char **v) {
	int w = 0, h = 0;
	if (c > 1)
		w = atoi(v[1]); ///< read width
	if (c > 2)
		h = atoi(v[2]); ///< read height
	if (w <= 0)
		w = 30; ///< default width
	if (h <= 0)
		h = 30; ///< default height

	game(w, h);

}
