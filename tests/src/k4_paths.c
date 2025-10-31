#include "../../src/graphs.h"
#include <stdio.h>
#include <stdlib.h>

int main() {
	Graph* g = graph_kn(4);

	for (size_t i = 1; i < g->n; i++) {
		int paths = (int) graph_paths_length(g, 0, 2, i);
		printf("número de passos de tamanho %ld = %d \n", i, paths);
	}

	printf("\n");

	free(g);
	return 0;
}