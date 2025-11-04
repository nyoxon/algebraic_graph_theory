#include "../../src/graphs.h"
#include "../../src/eig.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

void simulate(size_t n1, size_t n2, double p, size_t maxit) {
	srand(time(NULL));

	for (size_t i = 0; i < maxit; i++) {
		Graph* g = graph_random_bipartite(n1, n2, p);

		double* w = malloc(g->n * sizeof(double));

		graph_spec_adj(g, w);

		printf("\n==== n1 = %ld n2 = %ld ====\n", n1, n2);
		graph_print(g, "Matriz de adjacência");
		printf("\n");
		eigenvalues_print(w, g->n, "Autovalores");

		free(w);
		free(g);
	}
}

int main() {
	simulate(5, 3, 0.5, 10);

	return 0;
}
