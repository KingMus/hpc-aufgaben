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

long TimeSteps = 100;
long outputRank = 1;

void evolve(double* currentfield, double* newfield, int w, int h) {
	int x, y;
	for (y = 0; y < h; y++) {
		for (x = 0; x < w; x++) {

			//TODO FIXME impletent rules and assign new value

			newfield[calcIndex(w, x, y)] = !newfield[calcIndex(w, x, y)];
		}
	}
}

/*
 * --------------------------------------------------------- MAP GENERATION ---------------------------------------------------------
 */

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

/*
 * --------------------------------------------------------- MPI STUFF ---------------------------------------------------------
 */

void game(int w, int h) {
//	double *currentfield = calloc(w * h, sizeof(double));
	double *newfield = calloc(w * h, sizeof(double));

//	printf("size unsigned %d, size long %d\n",sizeof(float), sizeof(long));

//	filling(currentfield, w, h);

	double *currentfield = readFromFile("input_gol", w, h);

	int rank, size;

	int half_w = w / 2;
	int half_h = h / 2;

	MPI_Status mpi_status_up;
	MPI_Status mpi_status_down;
	MPI_Status mpi_status_right;
	MPI_Status mpi_status_left;

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

	//Output to check field, only for one process
	if (rank == outputRank) {
		print_Partfield(rank, w, h, part_field_with_ghost);
	}

	/*Für den Randaustausch benötigte Variablen.
	 * ------------------------------------------------------------------------------------------------------------------------------------
	 */

	int right_ghost_to_recieve[half_h];
	int left_ghost_to_recieve[half_h];
	int up_ghost_to_recieve[half_w];
	int down_ghost_to_recieve[half_w];

	int right_ghost_to_send[half_h];
	int left_ghost_to_send[half_h];
	int up_ghost_to_send[half_w];
	int down_ghost_to_send[half_w];

	MPI_Barrier(MPI_COMM_WORLD);

	/*Der von diesem Prozess zu sendende Ghostrand wird aus dem prozesseigenen part_field geschrieben.
	 * ------------------------------------------------------------------------------------------------------------------------------------
	 */

	for (int y = 0; y < half_h; y++) {
		up_ghost_to_send[y] = part_field_with_ghost[calcIndex(w, 1, y + 1)];

		down_ghost_to_send[y] = part_field_with_ghost[calcIndex(w, half_h,
				y + 1)];
	}

	for (int x = 0; x < half_w; x++) {
		left_ghost_to_send[x] = part_field_with_ghost[calcIndex(w, x + 1, 1)];
		right_ghost_to_send[x] = part_field_with_ghost[calcIndex(w, x + 1,
				half_h)];
	}

	//Output to check ghost layer, only for one process
	if (rank == outputRank) {
		printf("Die vier zu sendenen GhostLayer von %d rank\n", rank);
		for (int y = 0; y < half_h; y++) {
			printf("l%d ", left_ghost_to_send[y]);
			printf("r%d ", right_ghost_to_send[y]);
			printf("o%d ", up_ghost_to_send[y]);
			printf("u%d\n", down_ghost_to_send[y]);
		}

	}

	MPI_Barrier(MPI_COMM_WORLD);

	/* Jeder Prozess sendet jetzt seinen gerade erstellen Ghostrand seines part_fieldes.
	 * Außerdem recieved jeder Prozess den zu ihm passenden Ghostrand und speichert ihn in ...
	 * ------------------------------------------------------------------------------------------------------------------------------------
	 */

	if (rank == 0) {

//		print_ProcessInformation(rank, xStart, yStart, xEnd, yEnd, comm_cart);

		MPI_Recv(&up_ghost_to_recieve, half_w, MPI_INT, 2, 96, MPI_COMM_WORLD,
				&mpi_status_up);
		MPI_Recv(&left_ghost_to_recieve, half_h, MPI_INT, 1, 97,
		MPI_COMM_WORLD, &mpi_status_left);
		MPI_Recv(&down_ghost_to_recieve, half_w, MPI_INT, 2, 98,
		MPI_COMM_WORLD, &mpi_status_down);
		MPI_Recv(&right_ghost_to_recieve, half_h, MPI_INT, 1, 99,
		MPI_COMM_WORLD, &mpi_status_right);

		MPI_Send(&down_ghost_to_send, half_w, MPI_INT, 2, 96, MPI_COMM_WORLD);
		MPI_Send(&right_ghost_to_send, half_h, MPI_INT, 1, 97, MPI_COMM_WORLD);
		MPI_Send(&up_ghost_to_send, half_w, MPI_INT, 2, 98, MPI_COMM_WORLD);
		MPI_Send(&left_ghost_to_send, half_h, MPI_INT, 1, 99, MPI_COMM_WORLD);

	}
	if (rank == 1) {

//		print_ProcessInformation(rank, xStart, yStart, xEnd, yEnd, comm_cart);

		MPI_Send(&down_ghost_to_send, half_w, MPI_INT, 3, 96, MPI_COMM_WORLD);
		MPI_Send(&right_ghost_to_send, half_h, MPI_INT, 0, 97, MPI_COMM_WORLD);
		MPI_Send(&up_ghost_to_send, half_w, MPI_INT, 3, 98, MPI_COMM_WORLD);
		MPI_Send(&left_ghost_to_send, half_h, MPI_INT, 0, 99, MPI_COMM_WORLD);

		MPI_Recv(&up_ghost_to_recieve, half_w, MPI_INT, 3, 96, MPI_COMM_WORLD,
				&mpi_status_up);
		MPI_Recv(&left_ghost_to_recieve, half_h, MPI_INT, 0, 97,
		MPI_COMM_WORLD, &mpi_status_left);
		MPI_Recv(&down_ghost_to_recieve, half_w, MPI_INT, 3, 98,
		MPI_COMM_WORLD, &mpi_status_down);
		MPI_Recv(&right_ghost_to_recieve, half_h, MPI_INT, 0, 99,
		MPI_COMM_WORLD, &mpi_status_right);

	}
	if (rank == 2) {

//		print_ProcessInformation(rank, xStart, yStart, xEnd, yEnd, comm_cart);

		MPI_Send(&down_ghost_to_send, half_w, MPI_INT, 0, 96, MPI_COMM_WORLD);
		MPI_Send(&right_ghost_to_send, half_h, MPI_INT, 3, 97, MPI_COMM_WORLD);
		MPI_Send(&up_ghost_to_send, half_w, MPI_INT, 0, 98, MPI_COMM_WORLD);
		MPI_Send(&left_ghost_to_send, half_h, MPI_INT, 3, 99, MPI_COMM_WORLD);

		MPI_Recv(&up_ghost_to_recieve, half_w, MPI_INT, 0, 96, MPI_COMM_WORLD,
				&mpi_status_up);
		MPI_Recv(&left_ghost_to_recieve, half_h, MPI_INT, 3, 97,
		MPI_COMM_WORLD, &mpi_status_left);
		MPI_Recv(&down_ghost_to_recieve, half_w, MPI_INT, 0, 98,
		MPI_COMM_WORLD, &mpi_status_down);
		MPI_Recv(&right_ghost_to_recieve, half_h, MPI_INT, 3, 99,
		MPI_COMM_WORLD, &mpi_status_right);

	}
	if (rank == 3) {

//		print_ProcessInformation(rank, xStart, yStart, xEnd, yEnd, comm_cart);

		MPI_Recv(&up_ghost_to_recieve, half_w, MPI_INT, 1, 96, MPI_COMM_WORLD,
				&mpi_status_up);
		MPI_Recv(&left_ghost_to_recieve, half_h, MPI_INT, 2, 97,
		MPI_COMM_WORLD, &mpi_status_left);
		MPI_Recv(&down_ghost_to_recieve, half_w, MPI_INT, 1, 98,
		MPI_COMM_WORLD, &mpi_status_down);
		MPI_Recv(&right_ghost_to_recieve, half_h, MPI_INT, 2, 99,
		MPI_COMM_WORLD, &mpi_status_right);

		MPI_Send(&down_ghost_to_send, half_w, MPI_INT, 1, 96, MPI_COMM_WORLD);
		MPI_Send(&right_ghost_to_send, half_h, MPI_INT, 2, 97, MPI_COMM_WORLD);
		MPI_Send(&up_ghost_to_send, half_w, MPI_INT, 1, 98, MPI_COMM_WORLD);
		MPI_Send(&left_ghost_to_send, half_h, MPI_INT, 2, 99, MPI_COMM_WORLD);

	}

	MPI_Barrier(MPI_COMM_WORLD);

	//Output to check recieved ghost, only for one process
	if (rank == outputRank) {
		printf("Die vier erhaltenen GhostLayer von %d rank\n", rank);
		for (int y = 0; y < half_h; y++) {
			printf("l%d ", left_ghost_to_recieve[y]);
			printf("r%d ", right_ghost_to_recieve[y]);
			printf("o%d ", up_ghost_to_recieve[y]);
			printf("u%d\n", down_ghost_to_recieve[y]);
		}

	}

	/*Der erhaltene Ghostrand wird dem prozesseigenen part_field_with_ghost hinzugefügt.
	 * ------------------------------------------------------------------------------------------------------------------------------------
	 */

	for (int i = 0; i < half_w; i++) {

		part_field_with_ghost[calcIndex(w, i + 1, 0)] =
				left_ghost_to_recieve[i];
		part_field_with_ghost[calcIndex(w, i + 1, (half_w + 1))] =
				right_ghost_to_recieve[i];
		part_field_with_ghost[calcIndex(w, 0, i + 1)] = up_ghost_to_recieve[i];
		part_field_with_ghost[calcIndex(w, (half_w + 1), i + 1)] =
				down_ghost_to_recieve[i];
	}

	//Output to check field, only for one process
	if (rank == outputRank) {
		printf("\nDANACH\n");
		print_Partfield(rank, w, h, part_field_with_ghost);
	}

	/*

	 if (rank == 0) {
	 writeVTK2(TimeSteps, part_field_with_ghost, "gol", xEnd - xStart,
	 yEnd - yStart, "_r0_");
	 } else if (rank == 1) {
	 writeVTK2(TimeSteps, part_field_with_ghost, "gol", xEnd - xStart,
	 yEnd - yStart, "_r1_");
	 } else if (rank == 2) {
	 writeVTK2(TimeSteps, part_field_with_ghost, "gol", xEnd - xStart,
	 yEnd - yStart, "_r2_");
	 } else if (rank == 3) {
	 writeVTK2(TimeSteps, part_field_with_ghost, "gol", xEnd - xStart,
	 yEnd - yStart, "_r3_");
	 }

	 printf("\n Current Part field with Exchange \n");

	 //         evolve(current_part_field, new_part_field, w, h, 1, 15, 1, 15);

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
		w = 8; ///< default width
	if (h <= 0)
		h = 8; ///< default height

	MPI_Init(&c, &v);

	game(w, h);

	MPI_Finalize();
	return 0;

}

/*
 * --------------------------------------------------------- OUTPUT ---------------------------------------------------------
 */

void print_Partfield(int rank, int w, int h, int* part_field_with_ghost) {
	printf("Partfield-Ghost for %d rank\n", rank);
	for (int x = 0; x < ((w / 2) + 2); x++) {
		for (int y = 0; y < ((h / 2) + 2); y++) {
			printf("%d ", part_field_with_ghost[calcIndex(w, x, y)]);
		}
		printf("\n");
	}
	printf("Partfield ohne Ghost for %d rank\n", rank);
	for (int x = 0; x < (w / 2); x++) {
		for (int y = 0; y < (h / 2); y++) {
			printf("%d ", part_field_with_ghost[calcIndex(w, x + 1, y + 1)]);
		}
		printf("\n");
	}
}

void print_ProcessInformation(int rank, int xstart, int ystart, int xende,
		int yende, MPI_Comm comm_cart) {

	printf("--------------------------------------------\n");
	printf("My rank: %d\n", rank);
	printf("Section: von P(%d|%d) bis P(%d|%d) \n\n", xstart, ystart, xende,
			yende);

	printf("Position:\n");
	if (rank == 0) {
		printf("X0\n");
		printf("00\n");
	} else if (rank == 1) {
		printf("00\n");
		printf("X0\n");
	} else if (rank == 2) {
		printf("0X\n");
		printf("00\n");
	} else if (rank == 3) {
		printf("00\n");
		printf("0X\n");
	}

	int *rankleft;
	int *rankright;

	MPI_Cart_shift(comm_cart, 1, 1, &rankleft, &rankright);

	printf("Neighbour-ranks: %d (left) - %d (right) \n\n", rankleft, rankright);

}

/*
 * Bezueglich des Ghost-Randes:
 * Keine Ausgabe des Ghostrandes, der muss rausgerechnet werden.
 * Sonst sind die Randfelder doppelt vorhanden und jede Datei ist größer als sie sein muss.
 */

void writeVTK2(long timestep, double *data, char prefix[1024], long w, long h,
		char threadnum[1024]) {
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

	for (x = 0; x < w; x++) {
		for (y = 0; y < h; y++) {

			float value = data[calcIndex(h, y, x)];

			fwrite((unsigned char*) &value, sizeof(float), 1, fp);
		}
	}

	fprintf(fp, "\n</AppendedData>\n");
	fprintf(fp, "</VTKFile>\n");
	fclose(fp);
}

