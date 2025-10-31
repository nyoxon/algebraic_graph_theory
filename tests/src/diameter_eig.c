#include "../../src/graphs.h"
#include "../../src/eig.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

void simulate(size_t n, size_t n_times, double p) {
	srand(time(NULL));

	for (size_t i = 0; i < n_times; i++) {
		Graph* g = graph_random_connected(n, p);

		assert(graph_is_connected(g)); // importante

		double* w = malloc(n * sizeof(double));
		graph_spec_adj(g, w);

		int diameter = graph_diameter(g);

		printf("\n==== X ====\n");
		printf("diâmetro = %d", diameter);
		printf("\n");
		eigenvalues_print(w, n, NULL);
		printf("\n");

		free(w);
		free(g);
	}
}

int main() {
	simulate(10, 1000, 0.7);

	return 0;
}