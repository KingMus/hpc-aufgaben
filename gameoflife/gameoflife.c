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

void writeVTK2(long timestep, double *data, char prefix[1024], long w, long h) {
	char filename[2048];
	int x, y;

	long offsetX = 0;
	long offsetY = 0;
	float deltax = 1.0;
	float deltay = 1.0;
	long nxy = w * h * sizeof(float);

	mkdir("vtk_folder", 0777);

	snprintf(filename, sizeof(filename), "%s%s-%05ld%s", "vtk_folder/" , prefix, timestep,
			".vti");

	FILE* fp = fopen(filename, "w");

	fprintf(fp, "<?xml version=\"1.0\"?>\n");
	fprintf(fp,
			"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(fp,
			"<ImageData WholeExtent=\"%d %d %d %d %d %d\" Origin=\"0 0 0\" Spacing=\"%le %le %le\">\n",
			offsetX, offsetX + w, offsetY, offsetY + h, 0, 0, deltax, deltax,
			0.0);
	fprintf(fp, "<CellData Scalars=\"%s\">\n", prefix);
	fprintf(fp,
			"<DataArray type=\"Float32\" Name=\"%s\" format=\"appended\" offset=\"0\"/>\n",
			prefix);
	fprintf(fp, "</CellData>\n");
	fprintf(fp, "</ImageData>\n");
	fprintf(fp, "<AppendedData encoding=\"raw\">\n");
	fprintf(fp, "_");
	fwrite((unsigned char*) &nxy, sizeof(long), 1, fp);

	for (y = 0; y < h; y++) {
		for (x = 0; x < w; x++) {
			float value = data[calcIndex(h, x, y)];
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

void evolve(double* currentfield, double* newfield, int w, int h) {

	int x, y;

#pragma omp parallel for num_threads(12) collapse(2)
	for (y = 0; y < h; y++) {
		for (x = 0; x < w; x++) {

			int n = countLifingsPeriodic(currentfield, x, y, w, h);
			if (currentfield[calcIndex(w, x, y)])
				n--;

			newfield[calcIndex(w, x, y)] = (n == 3
					|| (n == 2 && currentfield[calcIndex(w, x, y)]));

		}
	}

}

void filling(double* currentfield, int w, int h) {
	int i;
	for (i = 0; i < h * w; i++) {
		currentfield[i] = (rand() < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
	}
}

double* readFromFile(char filename[256], int w, int h) {
	FILE* file = fopen(filename, "r"); //read mode

	double* field = calloc((w * h), sizeof(double));

	int size = w * h;
	int symbol;

	size_t len = 0;
	size_t width = 0;
	size_t height = 0;

	while ((symbol = fgetc(file)) != EOF) {
		if (symbol == '\n') {
			if (!width)
				width = len;
			height++;
			continue;
		}
		if (symbol == 'X')
			field[len++] = 1;
		if (symbol == '_')
			field[len++] = 0;

		// resize
		if (len == size) {
			field = realloc(field, sizeof(double) * (size += 10));
		}
	}
	height++;

	field = realloc(field, sizeof(*field) * len);

	fclose(file);
	return field;
}

void game(int w, int h) {
//	double *currentfield = calloc(w * h, sizeof(double));
	double *newfield = calloc(w * h, sizeof(double));

//printf("size unsigned %d, size long %d\n",sizeof(float), sizeof(long));

//	filling(currentfield, w, h);
	double *currentfield = readFromFile("input_gol", w, h);


	long t;
	for (t = 0; t < TimeSteps; t++) {
		show(currentfield, w, h);

		evolve(currentfield, newfield, w, h);

		printf("%ld timestep\n", t);

		writeVTK2(t, currentfield, "gol", w, h);

//		usleep(200000);

		//SWAP
		double *temp = currentfield;
		currentfield = newfield;
		newfield = temp;
	}

	free(currentfield);
	free(newfield);

	printf("DONE");

}

int main(int c, char **v) {
	int w = 0, h = 0;
	if (c > 1)
		w = atoi(v[1]); ///< read width
	if (c > 2)
		h = atoi(v[2]); ///< read height
	if (w <= 0)
		w = 10; ///< default width
	if (h <= 0)
		h = 10; ///< default height
	game(w, h);
}
