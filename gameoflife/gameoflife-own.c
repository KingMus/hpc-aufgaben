/*
 * gameoflife-own.c
 *
 *  Created on: 23.02.2018
 *      Author: marco
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define for_x for (int x = 0; x < w; x++)
#define for_y for (int y = 0; y < h; y++)
#define for_xy for_x for_y

long TimeSteps = 100;

void writeVTK2(long timestep, double *data, char prefix[1024], long w, long h) {
	char filename[2048];
	int x, y;

	long offsetX = 0;
	long offsetY = 0;
	float deltax = 1.0;
	float deltay = 1.0;
	long nxy = w * h * sizeof(float);

	snprintf(filename, sizeof(filename), "%s-%05ld%s", prefix, timestep,
			".vti");
	FILE* fp = fopen(filename, "w");

	fprintf(fp, "<?xml version=\"1.0\"?>\n");
	fprintf(fp,
			"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(fp,
			"<ImageData WholeExtent=\"%d %d %d %d %d %d\" Origin=\"0 0 0\" Spacing=\"%le %le %le\">\n",
			offsetX, offsetX + w - 1, offsetY, offsetY + h - 1, 0, 0, deltax,
			deltax, 0.0);
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


void show(void *u, int w, int h) {
	int (*univ)[w] = u;
	printf("\033[H");
	for_y
	{
		for_x
			printf(univ[y][x] ? "\033[07m  \033[m" : "  ");
		printf("\033[E");
	}
	fflush(stdout);
}

void evolve(void *u, int w, int h) {
	unsigned (*univ)[w] = u;
	unsigned new[h][w];

	for_y
		for_x
		{
			int n = 0;
			for (int y1 = y - 1; y1 <= y + 1; y1++)
				for (int x1 = x - 1; x1 <= x + 1; x1++)
					if (univ[(y1 + h) % h][(x1 + w) % w])
						n++;

			if (univ[y][x])
				n--;
			new[y][x] = (n == 3 || (n == 2 && univ[y][x]));
		}
	for_y
		for_x
			univ[y][x] = new[y][x];
}

void game(int w, int h) {
	unsigned univ[h][w];
	for_xy
		univ[y][x] = rand() < RAND_MAX / 10 ? 1 : 0;
	long t;
	for (t = 0; t < TimeSteps; t++) {
		show(univ, w, h);
		evolve(univ, w, h);
		writeVTK2(t, univ, "gol", w, h);
		usleep(200000);
	}
}

int main(int c, char **v) {
	int w = 0, h = 0;
	if (c > 1)
		w = atoi(v[1]);
	if (c > 2)
		h = atoi(v[2]);
	if (w <= 0)
		w = 30;
	if (h <= 0)
		h = 30;
	game(w, h);
}

