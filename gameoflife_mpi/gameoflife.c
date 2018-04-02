#include <endian.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <sys/time.h>

#include "/usr/include/mpi/mpi.h"

#define calcIndex(width, x,y)  ((y)*(width) + (x))

long TimeSteps = 100;

/*
 * Bezueglich des Ghost-Randes:
 * Keine Ausgabe des Ghostrandes, der muss rausgerechnet werden.
 * Sonst sind die Randfelder doppelt vorhanden und jede Datei ist größer als sie sein muss.
 */

void writeVTK2(long timestep, double *data, char prefix[1024], int wstart,
		int hstart, long w, long h, char threadnum[1024]) {
	char filename[2048];
	int x, y;

	long offsetX = 0;
	long offsetY = 0;
	float deltax = 1.0;
	float deltay = 1.0;
	long nxy = (w - wstart) * (h - hstart) * sizeof(float);

	mkdir("vtk_folder", 0777);

	snprintf(filename, sizeof(filename), "%s%s-%s%05ld%s","vtk_folder/", prefix, threadnum,
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

	for (x = wstart; x < w; x++) {
		for (y = hstart; y < h; y++) {

			int correctIndex = getCorrectIndex(w, h, threadnum);

			float value = data[calcIndex(correctIndex, y, x)];

			fwrite((unsigned char*) &value, sizeof(float), 1, fp);
		}
	}

	fprintf(fp, "\n</AppendedData>\n");
	fprintf(fp, "</VTKFile>\n");
	fclose(fp);
}

int getCorrectIndex(long w, long h, char threadnum[1024]) {
	long correctIndex = w;
	if (h > correctIndex) {
		correctIndex = h;
	}
	if (threadnum == "_r0_") {
		correctIndex = h * 2;
	}
	return correctIndex;
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

void evolve(double* currentfield, double* newfield, int w, int h) {
	int x, y;
	for (y = 0; y < h; y++) {
		for (x = 0; x < w; x++) {

			//TODO FIXME impletent rules and assign new value

			newfield[calcIndex(w, x, y)] = !newfield[calcIndex(w, x, y)];
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
	int rank, size;

	MPI_Status statusRight;
	MPI_Status statusLeft;

	MPI_Comm comm_cart;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int numProcessPerDimension[1], periodic[1], reorder;

	periodic[0] = 1;
	numProcessPerDimension[0] = size;
	reorder = 1;

	MPI_Cart_create(MPI_COMM_WORLD, 1, numProcessPerDimension, periodic,
			reorder, &comm_cart);

	if (rank == 0) {

		int xstart, ystart, xende, yende;

		xstart = 0;
		xende = w / 2;
		ystart = 0;
		yende = h / 2;

		printf("my rank: %d\n", rank);

		int *rankleft;
		int *rankright;

		MPI_Cart_shift(comm_cart, 1, 1, &rankleft, &rankright);

		printf("neighbour-ranks: %d (left) - %d (right) \n\n", rankleft,
				rankright);

		 writeVTK2(TimeSteps, currentfield, "gol", xstart, ystart, xende, yende,
		 "_r0_");


	}
	if (rank == 1) {

		int xstart, ystart, xende, yende;

		xstart = 0;
		xende = w / 2;
		ystart = h / 2;
		yende = h;

		printf("my rank: %d\n", rank);

		int * rankleft;
		int *rankright;

		MPI_Cart_shift(comm_cart, 1, 1, &rankleft, &rankright);

		printf("neighbour-ranks: %d (left) - %d (right) \n\n", rankleft,
				rankright);

				writeVTK2(TimeSteps, currentfield, "gol", xstart, ystart, xende, yende,
		 "_r1_");



	}
	if (rank == 2) {

		int xstart, ystart, xende, yende;

		xstart = w / 2;
		xende = w;
		ystart = 0;
		yende = h / 2;

		printf("my rank: %d\n", rank);

		int * rankleft;
		int *rankright;

		MPI_Cart_shift(comm_cart, 1, 1, &rankleft, &rankright);

		printf("neighbour-ranks: %d (left) - %d (right) \n\n", rankleft,
				rankright);

			writeVTK2(TimeSteps, currentfield, "gol", xstart, ystart, xende, yende,
		 "_r2_");



	}
	if (rank == 3) {

		int xstart, ystart, xende, yende;

		xstart = w / 2;
		xende = w;
		ystart = h / 2;
		yende = h;

		printf("my rank: %d\n", rank);

		int * rankleft;
		int *rankright;

		MPI_Cart_shift(comm_cart, 1, 1, &rankleft, &rankright);

		printf("neighbour-ranks: %d (left) - %d (right) \n\n", rankleft,
				rankright);

		writeVTK2(TimeSteps, currentfield, "gol", xstart, ystart, xende, yende,
		 "_r3_");



	}

	/*
	 for (t = 0; t < TimeSteps; t++) {
	 show(currentfield, w, h);
	 evolve(currentfield, newfield, w, h);

	 printf("%ld timestep\n", t);
	 writeVTK2(t, currentfield, "gol", w, h);

	 usleep(200000);

	 //SWAP
	 double *temp = currentfield;
	 currentfield = newfield;
	 newfield = temp;
	 }

	 */

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

	MPI_Init(&c, &v);

	game(w, h);

	MPI_Finalize();
	return 0;

}
