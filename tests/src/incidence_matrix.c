#include "../../src/graphs.h"
#include "../../src/eig.h"
#include <stdio.h>
#include <stdlib.h>

static void print_matrix(const double* B, size_t n, size_t m) {
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < m; j++) {
			printf("%4.0f ", B[IDX(i, j, m)]);
		}
		putchar('\n');
	}
}

void simulate() {
	double A[] = {
		0, 1, 0, 0, 0, 0,
		1, 0, 0, 0, 0, 0,
		0, 0, 0, 1, 0, 0,
		0, 0, 1, 0, 0, 0,
		0, 0, 0, 0, 0, 1,
		0, 0, 0, 0, 1, 0
	};

	Graph g = {6, true, A};

	double* B = graph_incidence_matrix(&g);
	print_matrix(B, g.n, graph_num_edges(&g));

	unsigned int rank = matrix_rank(B, g.n, graph_num_edges(&g), 10e-9);
	printf("rank de B: %d\n", rank);
	printf("componentes conexas: %d\n", (int) g.n - rank);

	free(B);
}

int main(void) {
	simulate();
	return 0;
}