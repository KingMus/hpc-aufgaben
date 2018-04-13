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

/* 0-3 für den Prozess für den der Output geschrieben werden soll.
 * dabei ist die Anordnung der Prozesse wie folgt:
 *
 * 0|1
 * 2|3
 */

long outputRank = 3;
// 1 für Prozessinformationen
int printProcessInformation = 0;
//1 für Output in der Konsole, der aufzeigt wie genau die einzelnen Schritte auf die Felder wirken
int printDebug = 1;

//Zählt alle Nachbarn, nichtperiodisch (es wird also nicht über den Rand geschaut)
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

//Oldfield wird nach den Regeln von Game of Life in Newfield entwickelt.
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

/*
 * --------------------------------------------------------- MAP GENERATION ---------------------------------------------------------
 */

//füllt das initiale currentfield
void filling(int* currentfield, int w, int h) {
	int i;
	for (i = 0; i < h * w; i++) {
		currentfield[i] = (rand() < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
	}
}

//liest eine Datei ein. w und h müssen mit der Dateigröße übereinstimmen
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

/*hier passiert die Action, zur Information: alles in game() wird von vier Prozessen gleichzeitig durchgeführt.
 * Die Anzahl 4 kommt aus dem Makefile bzw. dem Ausführen-Befehl. Mittels "rank" kann zwischen den Prozessen
 * unterschieden werden und es können individuelle Code-Segmente ausgeführt werden.
 *
 */
void game(int w, int h) {

//	int *currentfield = calloc(w * h, sizeof(double));
//	filling(currentfield, w, h); // aktivieren für Befüllung von Programm und nicht von Datei

	int *currentfield = readFromFile("input_gol", w, h);

	// Prozessrank und Size, MPI-Variablen
	int rank, size;

	int t = 0;

	//jeder Prozess bekommt ein Viertel des Feldes, zur Vereinfachung wurden deshalb diese Variablen eingeführt
	int half_w = w / 2;
	int half_h = h / 2;

	//Status, werden später im Austausch des Ghostlayers benötigt
	MPI_Status mpi_status_up;
	MPI_Status mpi_status_down;
	MPI_Status mpi_status_right;
	MPI_Status mpi_status_left;

	MPI_Status mpi_status_nbone;
	MPI_Status mpi_status_nbtwo;
	MPI_Status mpi_status_nbthree;

	//der im Code verwendete MPI-Communicator wird hier deklariert
	MPI_Comm comm_cart;

	//muss von MPI aus gemacht werden
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//Variablendeklarationen, die für Erzeugung unseres eigenen Communicators benötigt werden, siehe MPI_Cart_create
	int numProcessPerDimension[1], periodic[1], reorder;

	periodic[0] = 1;
	numProcessPerDimension[0] = size;
	reorder = 1;

	MPI_Cart_create(MPI_COMM_WORLD, 1, numProcessPerDimension, periodic,
			reorder, &comm_cart);

	//verschiedene, im Programmfluss benötigte Felder. Siehe im Code wofür genau
	int *part_field_with_ghost = calloc((half_w + 2) * (half_h + 2),
			sizeof(double));
	int *part_field = calloc(half_w * half_h, sizeof(double));
	int *part_field_next_gen_with_ghost = calloc((half_w + 2) * (half_h + 2),
			sizeof(double));

	/*Setze Start- und End-Values
	 * ------------------------------------------------------------------------------------------------------------------------------------
	 */

	int xStart, yStart, xEnd, yEnd;

	//Prozessspezifische Startvariablen, damit nicht jeder Prozess im selben Viertel der GoL-Map arbeitet.
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

	//Prozessinformationen für AB4-Teilaufgabe b)
	if (printProcessInformation) {
		print_ProcessInformation(rank, xStart, yStart, xEnd, yEnd, comm_cart);
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

		printf("Complete field in %ld timestep\n", t);
		for (int x = 0; x < w; x++) {
			for (int y = 0; y < h; y++) {
				printf("%d ", currentfield[calcIndex(w, x, y)]);
			}
			printf("\n");
		}

	}

	/*mit diesem Befehl warten alle Prozesse aufeinander, damit keiner weitermacht.
	 * Das macht Sinn, wenn erst eine Aufgabe von allen erledigt sein soll, damit keine Probleme bei der nächsten entstehen können.
	 * In diesem Programm nicht perfekt eingesetzt, geht bestimmt geschickter.
	 */
	MPI_Barrier(MPI_COMM_WORLD);

	//GoL begins und endet wenn alle finish auf 1 gesetzt sind.
	while (finish[0] == 0 || finish[1] == 0 || finish[2] == 0 || finish[3] == 0) {

		t++;

		//Output von finish
		if (rank == outputRank) {
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

		//Output to check field, only for one process
		if (rank == outputRank && printDebug) {
			printf("Partfield für VTK for %d rank\n", rank);
			for (int x = 0; x < (w / 2); x++) {
				for (int y = 0; y < (h / 2); y++) {
					printf("%d ", part_field[calcIndex(w, x, y)]);
				}
				printf("\n");
			}

		}

		//schreibe VTK, jeder Prozess schreibt seine eigene Datei, unterschieden wird durch den String
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

		//Output to check field, only for one process
		if (rank == outputRank && printDebug) {
			print_Partfield(rank, w, h, part_field_with_ghost);
		}

		/*Für den Randaustausch benötigte Variablen.
		 * ------------------------------------------------------------------------------------------------------------------------------------
		 */

		//Die Ghost-Ränder die erhalten werden (von anderen Prozessen)
		int right_ghost_to_recieve[half_h];
		int left_ghost_to_recieve[half_h];
		int up_ghost_to_recieve[half_w];
		int down_ghost_to_recieve[half_w];

		//Die Ghost-Ecken, die erhalten werden sollen (von anderen Prozessen)
		int down_left_ghost_to_recieve;
		int down_right_ghost_to_recieve;
		int up_left_ghost_to_recieve;
		int up_right_ghost_to_recieve;

		//Die Ghost-Ränder die gesendet werden sollen (von diesem Prozess)
		int right_ghost_to_send[half_h];
		int left_ghost_to_send[half_h];
		int up_ghost_to_send[half_w];
		int down_ghost_to_send[half_w];

		//Die Ghost-Ecken die gesendet werden sollen (von diesem Prozess)
		int down_left_ghost_to_send;
		int down_right_ghost_to_send;
		int up_left_ghost_to_send;
		int up_right_ghost_to_send;

		//warum das Ganze? Siehe explanation.md im Repository

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

		//Output to check ghost layer, only for one process
		if (rank == outputRank && printDebug) {
			printf("Die vier zu sendenen GhostLayer von %d rank\n", rank);
			for (int y = 0; y < half_h; y++) {
				printf("l%d ", left_ghost_to_send[y]);
				printf("r%d ", right_ghost_to_send[y]);
				printf("o%d ", up_ghost_to_send[y]);
				printf("u%d\n", down_ghost_to_send[y]);
			}

			printf("Die vier zu sendenen Ghost-Ecken von %d rank\n", rank);
			printf("lo%d ro%d\n", up_left_ghost_to_send,
					up_right_ghost_to_send);
			printf("lu%d ru%d\n", down_left_ghost_to_send,
					down_right_ghost_to_send);
		}

		MPI_Barrier(MPI_COMM_WORLD);

		/* Jeder Prozess sendet jetzt seinen gerade erstellten Ghostrand seines part_fieldes.
		 * Außerdem recieved jeder Prozess den zu ihm passenden Ghostrand und speichert ihn in ghost_to_recieve
		 * ------------------------------------------------------------------------------------------------------------------------------------
		 */

		if (rank == 0) {

			//Recieve Ghost-Rand
			MPI_Recv(&up_ghost_to_recieve, half_w, MPI_INT, 2, 96,
			MPI_COMM_WORLD, &mpi_status_up);
			MPI_Recv(&left_ghost_to_recieve, half_h, MPI_INT, 1, 97,
			MPI_COMM_WORLD, &mpi_status_left);
			MPI_Recv(&down_ghost_to_recieve, half_w, MPI_INT, 2, 98,
			MPI_COMM_WORLD, &mpi_status_down);
			MPI_Recv(&right_ghost_to_recieve, half_h, MPI_INT, 1, 99,
			MPI_COMM_WORLD, &mpi_status_right);

			//Send Ghost-Rand
			MPI_Send(&down_ghost_to_send, half_w, MPI_INT, 2, 96,
			MPI_COMM_WORLD);
			MPI_Send(&right_ghost_to_send, half_h, MPI_INT, 1, 97,
			MPI_COMM_WORLD);
			MPI_Send(&up_ghost_to_send, half_w, MPI_INT, 2, 98, MPI_COMM_WORLD);
			MPI_Send(&left_ghost_to_send, half_h, MPI_INT, 1, 99,
			MPI_COMM_WORLD);

			//Recieve Ghost-Ecken
			MPI_Recv(&up_left_ghost_to_recieve, 1, MPI_INT, 3, 96,
			MPI_COMM_WORLD, &mpi_status_up);
			MPI_Recv(&up_right_ghost_to_recieve, 1, MPI_INT, 3, 97,
			MPI_COMM_WORLD, &mpi_status_left);
			MPI_Recv(&down_left_ghost_to_recieve, 1, MPI_INT, 3, 98,
			MPI_COMM_WORLD, &mpi_status_down);
			MPI_Recv(&down_right_ghost_to_recieve, 1, MPI_INT, 3, 99,
			MPI_COMM_WORLD, &mpi_status_right);

			//Send Ghost-Ecken
			MPI_Send(&down_right_ghost_to_send, 1, MPI_INT, 3, 96,
			MPI_COMM_WORLD);
			MPI_Send(&down_left_ghost_to_send, 1, MPI_INT, 3, 97,
			MPI_COMM_WORLD);
			MPI_Send(&up_right_ghost_to_send, 1, MPI_INT, 3, 98,
			MPI_COMM_WORLD);
			MPI_Send(&up_left_ghost_to_send, 1, MPI_INT, 3, 99,
			MPI_COMM_WORLD);

		}
		if (rank == 1) {

			//Send Ghost-Rand
			MPI_Send(&down_ghost_to_send, half_w, MPI_INT, 3, 96,
			MPI_COMM_WORLD);
			MPI_Send(&right_ghost_to_send, half_h, MPI_INT, 0, 97,
			MPI_COMM_WORLD);
			MPI_Send(&up_ghost_to_send, half_w, MPI_INT, 3, 98, MPI_COMM_WORLD);
			MPI_Send(&left_ghost_to_send, half_h, MPI_INT, 0, 99,
			MPI_COMM_WORLD);

			//Recieve Ghost-Rand
			MPI_Recv(&up_ghost_to_recieve, half_w, MPI_INT, 3, 96,
			MPI_COMM_WORLD, &mpi_status_up);
			MPI_Recv(&left_ghost_to_recieve, half_h, MPI_INT, 0, 97,
			MPI_COMM_WORLD, &mpi_status_left);
			MPI_Recv(&down_ghost_to_recieve, half_w, MPI_INT, 3, 98,
			MPI_COMM_WORLD, &mpi_status_down);
			MPI_Recv(&right_ghost_to_recieve, half_h, MPI_INT, 0, 99,
			MPI_COMM_WORLD, &mpi_status_right);

			//Recieve Ghost-Ecken
			MPI_Recv(&up_left_ghost_to_recieve, 1, MPI_INT, 2, 96,
			MPI_COMM_WORLD, &mpi_status_up);
			MPI_Recv(&up_right_ghost_to_recieve, 1, MPI_INT, 2, 97,
			MPI_COMM_WORLD, &mpi_status_left);
			MPI_Recv(&down_left_ghost_to_recieve, 1, MPI_INT, 2, 98,
			MPI_COMM_WORLD, &mpi_status_down);
			MPI_Recv(&down_right_ghost_to_recieve, 1, MPI_INT, 2, 99,
			MPI_COMM_WORLD, &mpi_status_right);

			//Send Ghost-Ecken
			MPI_Send(&down_right_ghost_to_send, 1, MPI_INT, 2, 96,
			MPI_COMM_WORLD);
			MPI_Send(&down_left_ghost_to_send, 1, MPI_INT, 2, 97,
			MPI_COMM_WORLD);
			MPI_Send(&up_right_ghost_to_send, 1, MPI_INT, 2, 98,
			MPI_COMM_WORLD);
			MPI_Send(&up_left_ghost_to_send, 1, MPI_INT, 2, 99,
			MPI_COMM_WORLD);

		}
		if (rank == 2) {

			//Send Ghost-Rand
			MPI_Send(&down_ghost_to_send, half_w, MPI_INT, 0, 96,
			MPI_COMM_WORLD);
			MPI_Send(&right_ghost_to_send, half_h, MPI_INT, 3, 97,
			MPI_COMM_WORLD);
			MPI_Send(&up_ghost_to_send, half_w, MPI_INT, 0, 98, MPI_COMM_WORLD);
			MPI_Send(&left_ghost_to_send, half_h, MPI_INT, 3, 99,
			MPI_COMM_WORLD);

			//Recieve Ghost-Rand
			MPI_Recv(&up_ghost_to_recieve, half_w, MPI_INT, 0, 96,
			MPI_COMM_WORLD, &mpi_status_up);
			MPI_Recv(&left_ghost_to_recieve, half_h, MPI_INT, 3, 97,
			MPI_COMM_WORLD, &mpi_status_left);
			MPI_Recv(&down_ghost_to_recieve, half_w, MPI_INT, 0, 98,
			MPI_COMM_WORLD, &mpi_status_down);
			MPI_Recv(&right_ghost_to_recieve, half_h, MPI_INT, 3, 99,
			MPI_COMM_WORLD, &mpi_status_right);

			//Send Ghost-Ecken
			MPI_Send(&down_right_ghost_to_send, 1, MPI_INT, 1, 96,
			MPI_COMM_WORLD);
			MPI_Send(&down_left_ghost_to_send, 1, MPI_INT, 1, 97,
			MPI_COMM_WORLD);
			MPI_Send(&up_right_ghost_to_send, 1, MPI_INT, 1, 98,
			MPI_COMM_WORLD);
			MPI_Send(&up_left_ghost_to_send, 1, MPI_INT, 1, 99,
			MPI_COMM_WORLD);

			//Recieve Ghost-Ecken
			MPI_Recv(&up_left_ghost_to_recieve, 1, MPI_INT, 1, 96,
			MPI_COMM_WORLD, &mpi_status_up);
			MPI_Recv(&up_right_ghost_to_recieve, 1, MPI_INT, 1, 97,
			MPI_COMM_WORLD, &mpi_status_left);
			MPI_Recv(&down_left_ghost_to_recieve, 1, MPI_INT, 1, 98,
			MPI_COMM_WORLD, &mpi_status_down);
			MPI_Recv(&down_right_ghost_to_recieve, 1, MPI_INT, 1, 99,
			MPI_COMM_WORLD, &mpi_status_right);

		}
		if (rank == 3) {

			//Recieve Ghost-Rand
			MPI_Recv(&up_ghost_to_recieve, half_w, MPI_INT, 1, 96,
			MPI_COMM_WORLD, &mpi_status_up);
			MPI_Recv(&left_ghost_to_recieve, half_h, MPI_INT, 2, 97,
			MPI_COMM_WORLD, &mpi_status_left);
			MPI_Recv(&down_ghost_to_recieve, half_w, MPI_INT, 1, 98,
			MPI_COMM_WORLD, &mpi_status_down);
			MPI_Recv(&right_ghost_to_recieve, half_h, MPI_INT, 2, 99,
			MPI_COMM_WORLD, &mpi_status_right);

			//Send Ghost-Rand
			MPI_Send(&down_ghost_to_send, half_w, MPI_INT, 1, 96,
			MPI_COMM_WORLD);
			MPI_Send(&right_ghost_to_send, half_h, MPI_INT, 2, 97,
			MPI_COMM_WORLD);
			MPI_Send(&up_ghost_to_send, half_w, MPI_INT, 1, 98, MPI_COMM_WORLD);
			MPI_Send(&left_ghost_to_send, half_h, MPI_INT, 2, 99,
			MPI_COMM_WORLD);

			//Send Ghost-Ecken
			MPI_Send(&down_right_ghost_to_send, 1, MPI_INT, 0, 96,
			MPI_COMM_WORLD);
			MPI_Send(&down_left_ghost_to_send, 1, MPI_INT, 0, 97,
			MPI_COMM_WORLD);
			MPI_Send(&up_right_ghost_to_send, 1, MPI_INT, 0, 98,
			MPI_COMM_WORLD);
			MPI_Send(&up_left_ghost_to_send, 1, MPI_INT, 0, 99,
			MPI_COMM_WORLD);

			//Recieve Ghost-Ecken
			MPI_Recv(&up_left_ghost_to_recieve, 1, MPI_INT, 0, 96,
			MPI_COMM_WORLD, &mpi_status_up);
			MPI_Recv(&up_right_ghost_to_recieve, 1, MPI_INT, 0, 97,
			MPI_COMM_WORLD, &mpi_status_left);
			MPI_Recv(&down_left_ghost_to_recieve, 1, MPI_INT, 0, 98,
			MPI_COMM_WORLD, &mpi_status_down);
			MPI_Recv(&down_right_ghost_to_recieve, 1, MPI_INT, 0, 99,
			MPI_COMM_WORLD, &mpi_status_right);

		}

		MPI_Barrier(MPI_COMM_WORLD);

		//Output to check recieved ghost, only for one process
		if (rank == outputRank && printDebug) {
			printf("Die vier erhaltenen GhostLayer von %d rank\n", rank);
			for (int y = 0; y < half_h; y++) {
				printf("l%d ", left_ghost_to_recieve[y]);
				printf("r%d ", right_ghost_to_recieve[y]);
				printf("o%d ", up_ghost_to_recieve[y]);
				printf("u%d\n", down_ghost_to_recieve[y]);
			}

			printf("Die vier erhaltenen Ghost-Ecken von %d rank\n", rank);
			printf("lo%d ro%d\n", up_left_ghost_to_recieve,
					up_right_ghost_to_recieve);
			printf("lu%d ru%d\n", down_left_ghost_to_recieve,
					down_right_ghost_to_recieve);

		}

		/*Der erhaltene Ghostrand (und die Ecken) werden dem prozesseigenen part_field_with_ghost hinzugefügt.
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

		//Output to check field, only for one process
		if (rank == outputRank && printDebug) {
			printf("\nAFTER EXCHANGE\n\n");
			print_Partfield(rank, w, h, part_field_with_ghost);
		}

		MPI_Barrier(MPI_COMM_WORLD);

		/*Das aktuelle part_field_with_ghost wird eine Generation weiterentwickelt.
		 * ------------------------------------------------------------------------------------------------------------------------------------
		 */

		evolve(part_field_with_ghost, part_field_next_gen_with_ghost,
				(half_w + 2), (half_h + 2));

		//Output to check field, only for one process
		if (rank == outputRank && printDebug) {
			printf("\nAFTER EVOLVING\n\n");
			print_Partfield(rank, w, h, part_field_next_gen_with_ghost);

			printf(
					"\n%ld timestep -------------------------------------------------------------\n\n",
					t);
		}

		/*Prüfen, ob Veränderung im Feld vorliegt und mit anderen austauschen
		 * ------------------------------------------------------------------------------------------------------------------------------------
		 */

		finish[rank] = 1;

		//prüfe auf Veränderung
		for (int x = 0; x < (w / 2); x++) {
			for (int y = 0; y < (h / 2); y++) {

				if (part_field_with_ghost[calcIndex(w, x + 1, y + 1)]
						!= part_field_next_gen_with_ghost[calcIndex(w, x + 1,
								y + 1)]) {
					finish[rank] = 0;
				}
			}
		}

		//jeder Prozess schreibt seine eigene finish-Variable, wenn auch an der richtigen Stelle.
		//die anderen erhält er erneut durch send und recieve
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

		MPI_Barrier(MPI_COMM_WORLD);

		//SWAP
		int *temp = part_field_with_ghost;
		part_field_with_ghost = part_field_next_gen_with_ghost;
		part_field_next_gen_with_ghost = temp;

//		usleep(200000);

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

/*
 * --------------------------------------------------------- OUTPUT METHODS ---------------------------------------------------------
 */

//damit der Code nicht andauernd doppelt vorkommt wurde er ausgelagert
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

//für Teilaufgabe b) auf AB4
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
 *
 * Komischerweise sind die erzeugten VTK-Dateien horizontal gespiegelt. Ich hab aber keine Lust mehr das zu fixen.
 */

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

