#include "graphs.h"
#include "eig.h"
#include <stdio.h>
#include <stdbool.h>

int main(void) {
	Graph* g = graph_kn(10);

	double w[g->n];
	graph_spec_adj(g, w);

	eigenvalues_print(w, g->n, NULL);

	int diameter = graph_diameter(g);

	printf("diameter: %d\n", diameter);

	return 0;
}