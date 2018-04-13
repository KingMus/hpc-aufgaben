#include <endian.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <sys/time.h>

#include "/usr/include/mpi/mpi.h"

#define calcIndex(width, x,y)  ((x)*(width) + (y))

int finish[3];

void writeVTK2(long timestep, int *data, char prefix[1024], long w, long h,
		char threadnum[1024], long indexForCalc, int x1, int y1) {
	char filename[2048];
	int x, y;

	long offsetX = 0;
	long offsetY = 0;
	float deltax = 1.0;
	float deltay = 1.0;
	long nxy = w * h * sizeof(float);

	mkdir("vtk_folder", 0777);

	snprintf(filename, sizeof(filename), "%s%s-%s%05ld%s", "vtk_folder/",
			prefix, threadnum, timestep, ".vti");

	FILE* fp = fopen(filename, "w");

	fprintf(fp, "<?xml version=\"1.0\"?>\n");
	fprintf(fp,
			"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");

	fprintf(fp,
			"<ImageData WholeExtent=\"%d %d %d %d %d %d\" Origin=\"%d %d 0\" Spacing=\"%le %le %le\">\n",
			offsetX, offsetX + w, offsetY, offsetY + h, 0, 0, x1, y1, deltax,
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

	for (x = 0; x < w; x++) {
		for (y = 0; y < h; y++) {

			float value = data[calcIndex(indexForCalc, x, y)];

			fwrite((unsigned char*) &value, sizeof(float), 1, fp);
		}
	}

	fprintf(fp, "\n</AppendedData>\n");
	fprintf(fp, "</VTKFile>\n");
	fclose(fp);
}

int countLifings(int* field, int x, int y, int w, int h) {
	int neighbours = 0;

	for (int x1 = x - 1; x1 <= x + 1; x1++) {
		for (int y1 = y - 1; y1 <= y + 1; y1++) {
			if (field[calcIndex((w - 2) * 2, x1, y1)]) {
				neighbours++;
			}
		}
	}
	return neighbours;

}

void evolve(int* oldfield, int* newfield, int w, int h) {

	int x, y;
	for (x = 1; x < w - 1; x++) {
		for (y = 1; y < h - 1; y++) {

			int n = countLifings(oldfield, x, y, w, h);
			if (oldfield[calcIndex((w - 2) * 2, x, y)])
				n--;

			newfield[calcIndex((w - 2) * 2, x, y)] = (n == 3
					|| (n == 2 && oldfield[calcIndex((w - 2) * 2, x, y)]));
		}
	}
}

int* readFromFile(char filename[256], int w, int h) {
	FILE* file = fopen(filename, "r"); //read mode

	int* field = calloc((w * h), sizeof(double));

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

/*
 * --------------------------------------------------------- MPI STUFF ---------------------------------------------------------
 */

void game(int w, int h) {

	int *currentfield = readFromFile("input_gol", w, h);

	int rank, size;

	int t = 0;

	int half_w = w / 2;
	int half_h = h / 2;

	MPI_Status mpi_status_up;
	MPI_Status mpi_status_down;
	MPI_Status mpi_status_right;
	MPI_Status mpi_status_left;

	MPI_Status mpi_status_nbone;
	MPI_Status mpi_status_nbtwo;
	MPI_Status mpi_status_nbthree;

	MPI_Comm comm_cart;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int numProcessPerDimension[1], periodic[1], reorder;

	periodic[0] = 1;
	numProcessPerDimension[0] = size;
	reorder = 1;

	MPI_Cart_create(MPI_COMM_WORLD, 1, numProcessPerDimension, periodic,
			reorder, &comm_cart);

	int *part_field_with_ghost = calloc((half_w + 2) * (half_h + 2),
			sizeof(double));
	int *part_field = calloc(half_w * half_h, sizeof(double));

	/*Setze Start- und End-Values
	 * ------------------------------------------------------------------------------------------------------------------------------------
	 */

	int xStart, yStart, xEnd, yEnd;

	if (rank == 0) {
		xStart = 0;
		xEnd = half_w;
		yStart = 0;
		yEnd = half_h;
	} else if (rank == 1) {

		xStart = 0;
		xEnd = half_w;
		yStart = half_h;
		yEnd = h;
	} else if (rank == 2) {

		xStart = half_w;
		xEnd = w;
		yStart = 0;
		yEnd = half_h;
	} else if (rank == 3) {
		xStart = half_w;
		xEnd = w;
		yStart = half_h;
		yEnd = h;
	}

	/*Für jeden Prozess wird das passende Teilfeld aus dem großen currentfield geschrieben (ohne Ghostrand-Befüllung)
	 * ------------------------------------------------------------------------------------------------------------------------------------
	 */

	for (int x = 0; x < half_w; x++) {
		for (int y = 0; y < half_h; y++) {
			part_field_with_ghost[calcIndex(w, x + 1, y + 1)] =
					currentfield[calcIndex(w, x + xStart, y + yStart)];
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	//GoL begins
	//for (t = 0; t < TimeSteps; t++) {

	while (finish[0] == 0 || finish[1] == 0 || finish[2] == 0 || finish[3] == 0) {

		t++;

		if (rank == 0) {
			printf("aktueller Timestep: %d - finish: %d%d%d%d\n", t, finish[0],
					finish[1], finish[2], finish[3]);
		}

		/*Das aktuelle part_field_with_ghost wird in ein Feld ohne ghost geschrieben und dann in die VTK geschrieben
		 * ------------------------------------------------------------------------------------------------------------------------------------
		 */

		for (int x = 0; x < (w / 2); x++) {
			for (int y = 0; y < (h / 2); y++) {
				part_field[calcIndex(w, x, y)] =
						part_field_with_ghost[calcIndex(w, x + 1, y + 1)];
			}
		}

		char specific_filename[1024];
		sprintf(specific_filename, "r%d-", rank);

		if (rank == 0) {
			writeVTK2(t, part_field, "gol", half_w, half_h, "r0-", w, 0, 0);
		}
		if (rank == 1) {
			writeVTK2(t, part_field, "gol", half_w, half_h, "r1-", w, 3, 0);
		}
		if (rank == 2) {
			writeVTK2(t, part_field, "gol", half_w, half_h, "r2-", w, 0, 3);
		}
		if (rank == 3) {
			writeVTK2(t, part_field, "gol", half_w, half_h, "r3-", w, 3, 3);
		}

		/*Für den Randaustausch benötigte Variablen.
		 * ------------------------------------------------------------------------------------------------------------------------------------
		 */

		int right_ghost_to_recieve[half_h];
		int left_ghost_to_recieve[half_h];
		int up_ghost_to_recieve[half_w];
		int down_ghost_to_recieve[half_w];

		int down_left_ghost_to_recieve;
		int down_right_ghost_to_recieve;
		int up_left_ghost_to_recieve;
		int up_right_ghost_to_recieve;

		int right_ghost_to_send[half_h];
		int left_ghost_to_send[half_h];
		int up_ghost_to_send[half_w];
		int down_ghost_to_send[half_w];

		int down_left_ghost_to_send;
		int down_right_ghost_to_send;
		int up_left_ghost_to_send;
		int up_right_ghost_to_send;

		MPI_Barrier(MPI_COMM_WORLD);

		/*Der von diesem Prozess zu sendende Ghostrand wird aus dem prozesseigenen part_field geschrieben. Das gilt auch für die Ghost-Ecken.
		 * ------------------------------------------------------------------------------------------------------------------------------------
		 */

		for (int y = 0; y < half_h; y++) {
			up_ghost_to_send[y] = part_field_with_ghost[calcIndex(w, 1, y + 1)];

			down_ghost_to_send[y] = part_field_with_ghost[calcIndex(w, half_h,
					y + 1)];
		}

		for (int x = 0; x < half_w; x++) {
			left_ghost_to_send[x] =
					part_field_with_ghost[calcIndex(w, x + 1, 1)];
			right_ghost_to_send[x] = part_field_with_ghost[calcIndex(w, x + 1,
					half_h)];
		}

		down_left_ghost_to_send =
				part_field_with_ghost[calcIndex(w, half_h, 1)];
		down_right_ghost_to_send = part_field_with_ghost[calcIndex(w, half_w,
				half_h)];
		up_left_ghost_to_send = part_field_with_ghost[calcIndex(w, 1, 1)];
		up_right_ghost_to_send = part_field_with_ghost[calcIndex(w, 1, half_w)];

		MPI_Barrier(MPI_COMM_WORLD);

		/* Jeder Prozess sendet jetzt seinen gerade erstellten Ghostrand seines part_fieldes.
		 * Außerdem recieved jeder Prozess den zu ihm passenden Ghostrand und speichert ihn in ghost_to_recieve
		 * ------------------------------------------------------------------------------------------------------------------------------------
		 */

		if (rank == 0) {

			//Recieve Ghost-Rand
			MPI_Recv(&up_ghost_to_recieve, half_w, MPI_INT, 2, 10,
			MPI_COMM_WORLD, &mpi_status_up);
			MPI_Recv(&left_ghost_to_recieve, half_h, MPI_INT, 1, 11,
			MPI_COMM_WORLD, &mpi_status_left);
			MPI_Recv(&down_ghost_to_recieve, half_w, MPI_INT, 2, 12,
			MPI_COMM_WORLD, &mpi_status_down);
			MPI_Recv(&right_ghost_to_recieve, half_h, MPI_INT, 1, 13,
			MPI_COMM_WORLD, &mpi_status_right);

			//Send Ghost-Rand
			MPI_Send(&down_ghost_to_send, half_w, MPI_INT, 2, 10,
			MPI_COMM_WORLD);
			MPI_Send(&right_ghost_to_send, half_h, MPI_INT, 1, 11,
			MPI_COMM_WORLD);
			MPI_Send(&up_ghost_to_send, half_w, MPI_INT, 2, 12, MPI_COMM_WORLD);
			MPI_Send(&left_ghost_to_send, half_h, MPI_INT, 1, 13,
			MPI_COMM_WORLD);

			//Recieve Ghost-Ecken
			MPI_Recv(&up_left_ghost_to_recieve, 1, MPI_INT, 3, 10,
			MPI_COMM_WORLD, &mpi_status_up);
			MPI_Recv(&up_right_ghost_to_recieve, 1, MPI_INT, 3, 11,
			MPI_COMM_WORLD, &mpi_status_left);
			MPI_Recv(&down_left_ghost_to_recieve, 1, MPI_INT, 3, 12,
			MPI_COMM_WORLD, &mpi_status_down);
			MPI_Recv(&down_right_ghost_to_recieve, 1, MPI_INT, 3, 13,
			MPI_COMM_WORLD, &mpi_status_right);

			//Send Ghost-Ecken
			MPI_Send(&down_right_ghost_to_send, 1, MPI_INT, 3, 10,
			MPI_COMM_WORLD);
			MPI_Send(&down_left_ghost_to_send, 1, MPI_INT, 3, 11,
			MPI_COMM_WORLD);
			MPI_Send(&up_right_ghost_to_send, 1, MPI_INT, 3, 12,
			MPI_COMM_WORLD);
			MPI_Send(&up_left_ghost_to_send, 1, MPI_INT, 3, 13,
			MPI_COMM_WORLD);

		}
		if (rank == 1) {

			//Send Ghost-Rand
			MPI_Send(&down_ghost_to_send, half_w, MPI_INT, 3, 10,
			MPI_COMM_WORLD);
			MPI_Send(&right_ghost_to_send, half_h, MPI_INT, 0, 11,
			MPI_COMM_WORLD);
			MPI_Send(&up_ghost_to_send, half_w, MPI_INT, 3, 12, MPI_COMM_WORLD);
			MPI_Send(&left_ghost_to_send, half_h, MPI_INT, 0, 13,
			MPI_COMM_WORLD);

			//Recieve Ghost-Rand
			MPI_Recv(&up_ghost_to_recieve, half_w, MPI_INT, 3, 10,
			MPI_COMM_WORLD, &mpi_status_up);
			MPI_Recv(&left_ghost_to_recieve, half_h, MPI_INT, 0, 11,
			MPI_COMM_WORLD, &mpi_status_left);
			MPI_Recv(&down_ghost_to_recieve, half_w, MPI_INT, 3, 12,
			MPI_COMM_WORLD, &mpi_status_down);
			MPI_Recv(&right_ghost_to_recieve, half_h, MPI_INT, 0, 13,
			MPI_COMM_WORLD, &mpi_status_right);

			//Recieve Ghost-Ecken
			MPI_Recv(&up_left_ghost_to_recieve, 1, MPI_INT, 2, 10,
			MPI_COMM_WORLD, &mpi_status_up);
			MPI_Recv(&up_right_ghost_to_recieve, 1, MPI_INT, 2, 11,
			MPI_COMM_WORLD, &mpi_status_left);
			MPI_Recv(&down_left_ghost_to_recieve, 1, MPI_INT, 2, 12,
			MPI_COMM_WORLD, &mpi_status_down);
			MPI_Recv(&down_right_ghost_to_recieve, 1, MPI_INT, 2, 13,
			MPI_COMM_WORLD, &mpi_status_right);

			//Send Ghost-Ecken
			MPI_Send(&down_right_ghost_to_send, 1, MPI_INT, 2, 10,
			MPI_COMM_WORLD);
			MPI_Send(&down_left_ghost_to_send, 1, MPI_INT, 2, 11,
			MPI_COMM_WORLD);
			MPI_Send(&up_right_ghost_to_send, 1, MPI_INT, 2, 12,
			MPI_COMM_WORLD);
			MPI_Send(&up_left_ghost_to_send, 1, MPI_INT, 2, 13,
			MPI_COMM_WORLD);

		}
		if (rank == 2) {

			//Send Ghost-Rand
			MPI_Send(&down_ghost_to_send, half_w, MPI_INT, 0, 10,
			MPI_COMM_WORLD);
			MPI_Send(&right_ghost_to_send, half_h, MPI_INT, 3, 11,
			MPI_COMM_WORLD);
			MPI_Send(&up_ghost_to_send, half_w, MPI_INT, 0, 12, MPI_COMM_WORLD);
			MPI_Send(&left_ghost_to_send, half_h, MPI_INT, 3, 13,
			MPI_COMM_WORLD);

			//Recieve Ghost-Rand
			MPI_Recv(&up_ghost_to_recieve, half_w, MPI_INT, 0, 10,
			MPI_COMM_WORLD, &mpi_status_up);
			MPI_Recv(&left_ghost_to_recieve, half_h, MPI_INT, 3, 11,
			MPI_COMM_WORLD, &mpi_status_left);
			MPI_Recv(&down_ghost_to_recieve, half_w, MPI_INT, 0, 12,
			MPI_COMM_WORLD, &mpi_status_down);
			MPI_Recv(&right_ghost_to_recieve, half_h, MPI_INT, 3, 13,
			MPI_COMM_WORLD, &mpi_status_right);

			//Send Ghost-Ecken
			MPI_Send(&down_right_ghost_to_send, 1, MPI_INT, 1, 10,
			MPI_COMM_WORLD);
			MPI_Send(&down_left_ghost_to_send, 1, MPI_INT, 1, 11,
			MPI_COMM_WORLD);
			MPI_Send(&up_right_ghost_to_send, 1, MPI_INT, 1, 12,
			MPI_COMM_WORLD);
			MPI_Send(&up_left_ghost_to_send, 1, MPI_INT, 1, 13,
			MPI_COMM_WORLD);

			//Recieve Ghost-Ecken
			MPI_Recv(&up_left_ghost_to_recieve, 1, MPI_INT, 1, 10,
			MPI_COMM_WORLD, &mpi_status_up);
			MPI_Recv(&up_right_ghost_to_recieve, 1, MPI_INT, 1, 11,
			MPI_COMM_WORLD, &mpi_status_left);
			MPI_Recv(&down_left_ghost_to_recieve, 1, MPI_INT, 1, 12,
			MPI_COMM_WORLD, &mpi_status_down);
			MPI_Recv(&down_right_ghost_to_recieve, 1, MPI_INT, 1, 13,
			MPI_COMM_WORLD, &mpi_status_right);

		}
		if (rank == 3) {

			//Recieve Ghost-Rand
			MPI_Recv(&up_ghost_to_recieve, half_w, MPI_INT, 1, 10,
			MPI_COMM_WORLD, &mpi_status_up);
			MPI_Recv(&left_ghost_to_recieve, half_h, MPI_INT, 2, 11,
			MPI_COMM_WORLD, &mpi_status_left);
			MPI_Recv(&down_ghost_to_recieve, half_w, MPI_INT, 1, 12,
			MPI_COMM_WORLD, &mpi_status_down);
			MPI_Recv(&right_ghost_to_recieve, half_h, MPI_INT, 2, 13,
			MPI_COMM_WORLD, &mpi_status_right);

			//Send Ghost-Rand
			MPI_Send(&down_ghost_to_send, half_w, MPI_INT, 1, 10,
			MPI_COMM_WORLD);
			MPI_Send(&right_ghost_to_send, half_h, MPI_INT, 2, 11,
			MPI_COMM_WORLD);
			MPI_Send(&up_ghost_to_send, half_w, MPI_INT, 1, 12, MPI_COMM_WORLD);
			MPI_Send(&left_ghost_to_send, half_h, MPI_INT, 2, 13,
			MPI_COMM_WORLD);

			//Send Ghost-Ecken
			MPI_Send(&down_right_ghost_to_send, 1, MPI_INT, 0, 10,
			MPI_COMM_WORLD);
			MPI_Send(&down_left_ghost_to_send, 1, MPI_INT, 0, 11,
			MPI_COMM_WORLD);
			MPI_Send(&up_right_ghost_to_send, 1, MPI_INT, 0, 12,
			MPI_COMM_WORLD);
			MPI_Send(&up_left_ghost_to_send, 1, MPI_INT, 0, 13,
			MPI_COMM_WORLD);

			//Recieve Ghost-Ecken
			MPI_Recv(&up_left_ghost_to_recieve, 1, MPI_INT, 0, 10,
			MPI_COMM_WORLD, &mpi_status_up);
			MPI_Recv(&up_right_ghost_to_recieve, 1, MPI_INT, 0, 11,
			MPI_COMM_WORLD, &mpi_status_left);
			MPI_Recv(&down_left_ghost_to_recieve, 1, MPI_INT, 0, 12,
			MPI_COMM_WORLD, &mpi_status_down);
			MPI_Recv(&down_right_ghost_to_recieve, 1, MPI_INT, 0, 13,
			MPI_COMM_WORLD, &mpi_status_right);

		}

		MPI_Barrier(MPI_COMM_WORLD);

		/*Der erhaltene Ghostrand wird dem prozesseigenen part_field_with_ghost hinzugefügt.
		 * ------------------------------------------------------------------------------------------------------------------------------------
		 */

		for (int i = 0; i < half_w; i++) {

			part_field_with_ghost[calcIndex(w, i + 1, 0)] =
					left_ghost_to_recieve[i];
			part_field_with_ghost[calcIndex(w, i + 1, (half_w + 1))] =
					right_ghost_to_recieve[i];
			part_field_with_ghost[calcIndex(w, 0, i + 1)] =
					up_ghost_to_recieve[i];
			part_field_with_ghost[calcIndex(w, (half_w + 1), i + 1)] =
					down_ghost_to_recieve[i];
		}

		part_field_with_ghost[calcIndex(w, 0, 0)] = up_left_ghost_to_recieve;
		part_field_with_ghost[calcIndex(w, 0, half_w + 1)] =
				up_right_ghost_to_recieve;
		part_field_with_ghost[calcIndex(w, half_h + 1, 0)] =
				down_left_ghost_to_recieve;
		part_field_with_ghost[calcIndex(w, half_w + 1, half_h + 1)] =
				down_right_ghost_to_recieve;

		MPI_Barrier(MPI_COMM_WORLD);

		/*Das aktuelle part_field_with_ghost wird eine Generation weiterentwickelt.
		 * ------------------------------------------------------------------------------------------------------------------------------------
		 */

		int *part_field_next_gen_with_ghost = calloc(
				(half_w + 2) * (half_h + 2), sizeof(double));

		evolve(part_field_with_ghost, part_field_next_gen_with_ghost,
				(half_w + 2), (half_h + 2));

		MPI_Barrier(MPI_COMM_WORLD);

		/*Prüfen, ob Veränderung im Feld vorliegt und mit anderen austauschen
		 * ------------------------------------------------------------------------------------------------------------------------------------
		 */

		finish[rank] = 1;

		for (int x = 0; x < (w / 2); x++) {
			for (int y = 0; y < (h / 2); y++) {

				if (part_field_with_ghost[calcIndex(w, x + 1, y + 1)]
						!= part_field_next_gen_with_ghost[calcIndex(w, x + 1,
								y + 1)]) {
					finish[rank] = 0;
				}
			}
		}

		if (rank == 0) {

			MPI_Send(&finish[0], 1, MPI_INT, 1, 10, MPI_COMM_WORLD);
			MPI_Recv(&finish[1], 1, MPI_INT, 1, 10, MPI_COMM_WORLD,
					&mpi_status_nbone);

			MPI_Send(&finish[0], 1, MPI_INT, 2, 11, MPI_COMM_WORLD);
			MPI_Recv(&finish[2], 1, MPI_INT, 2, 11, MPI_COMM_WORLD,
					&mpi_status_nbtwo);

			MPI_Send(&finish[0], 1, MPI_INT, 3, 12, MPI_COMM_WORLD);
			MPI_Recv(&finish[3], 1, MPI_INT, 3, 12, MPI_COMM_WORLD,
					&mpi_status_nbthree);

		}
		if (rank == 1) {

			MPI_Recv(&finish[0], 1, MPI_INT, 0, 10, MPI_COMM_WORLD,
					&mpi_status_nbone);
			MPI_Send(&finish[1], 1, MPI_INT, 0, 10, MPI_COMM_WORLD);

			MPI_Recv(&finish[3], 1, MPI_INT, 3, 11, MPI_COMM_WORLD,
					&mpi_status_nbthree);
			MPI_Send(&finish[1], 1, MPI_INT, 3, 11, MPI_COMM_WORLD);

			MPI_Recv(&finish[2], 1, MPI_INT, 2, 12, MPI_COMM_WORLD,
					&mpi_status_nbtwo);
			MPI_Send(&finish[1], 1, MPI_INT, 2, 12, MPI_COMM_WORLD);

		}
		if (rank == 2) {

			MPI_Send(&finish[2], 1, MPI_INT, 3, 10, MPI_COMM_WORLD);
			MPI_Recv(&finish[3], 1, MPI_INT, 3, 10, MPI_COMM_WORLD,
					&mpi_status_nbthree);

			MPI_Recv(&finish[0], 1, MPI_INT, 0, 11, MPI_COMM_WORLD,
					&mpi_status_nbone);
			MPI_Send(&finish[2], 1, MPI_INT, 0, 11, MPI_COMM_WORLD);

			MPI_Send(&finish[2], 1, MPI_INT, 1, 12, MPI_COMM_WORLD);
			MPI_Recv(&finish[1], 1, MPI_INT, 1, 12, MPI_COMM_WORLD,
					&mpi_status_nbtwo);

		}
		if (rank == 3) {

			MPI_Recv(&finish[2], 1, MPI_INT, 2, 10, MPI_COMM_WORLD,
					&mpi_status_nbthree);
			MPI_Send(&finish[3], 1, MPI_INT, 2, 10, MPI_COMM_WORLD);

			MPI_Send(&finish[3], 1, MPI_INT, 1, 11, MPI_COMM_WORLD);
			MPI_Recv(&finish[1], 1, MPI_INT, 1, 11, MPI_COMM_WORLD,
					&mpi_status_nbtwo);

			MPI_Recv(&finish[0], 1, MPI_INT, 0, 12, MPI_COMM_WORLD,
					&mpi_status_nbone);
			MPI_Send(&finish[3], 1, MPI_INT, 0, 12, MPI_COMM_WORLD);

		}

		MPI_Barrier(MPI_COMM_WORLD);

		//SWAP
		int *temp = part_field_with_ghost;
		part_field_with_ghost = part_field_next_gen_with_ghost;
		part_field_next_gen_with_ghost = temp;

		usleep(200000);

	}

	printf("Process %d: DONE\n", rank);
	free(currentfield);

}

int main(int c, char **v) {
	int w = 0, h = 0;
	if (c > 1)
		w = atoi(v[1]); ///< read width
	if (c > 2)
		h = atoi(v[2]); ///< read height
	if (w <= 0)
		w = 8; ///< default width
	if (h <= 0)
		h = 8; ///< default height

	MPI_Init(&c, &v);

	game(w, h);

	MPI_Finalize();
	return 0;

}
