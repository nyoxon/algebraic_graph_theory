#include "../../src/graphs.h"
#include "../../src/eig.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void simulate(size_t n) {
	srand(time(NULL));

	for (size_t k = 1; k <= n; k++) {
		Graph* g = graph_random_regular(n, k);

		if (!g) {
			continue;
		}

		double* w = malloc(n * sizeof(double));

		graph_spec_adj(g, w);

		printf("\n==== REGULAR %ld VÉRTICES, grau %ld ====\n", n, k);
		printf("%s\n", graph_is_connected(g)? "CONEXO" : "DISCONEXO");
		eigenvalues_print(w, n, NULL);
		free(w);
		free(g);
	}
}

int main() {
	simulate(10);

	return 0;
}