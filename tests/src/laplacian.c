#include "../../src/graphs.h"
#include "../../src/eig.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <math.h>

static inline int feq(double a, double b, double rtol, double atol) {
    double diff = fabs(a - b);
    double scale = fmax(fabs(a), fabs(b));
    return diff <= (atol + rtol * scale);
}

static unsigned int count_mult
(
        double* eigs, 
        double eig, 
        size_t n,
        double rtol,
        double atol
) {
        unsigned int count = 0;

        for (size_t i = 0; i < n; i++) {
                if (feq(eigs[i], eig, rtol, atol)) {
                        count++;
                }
        }

        return count;
}


void simulate(size_t n, size_t maxit, double p) {
	srand(time(NULL));

	for (size_t i = 0; i < maxit; i++) {
		Graph* g = graph_random(n, p);
		double* x = malloc(g->n * sizeof(double));

		printf("\n%s\n", (graph_is_connected(g))? "CONEXO" : "DISCONEXO");
		graph_spec_lap(g, x);
		eigenvalues_print(x, g->n, "Autovalores da laplaciana");
		double rtol = 1e-9;
		double atol = 1e-12;
		unsigned int count_zero = count_mult(x, 0.0, g->n, rtol, atol);
		printf("número de componentes conexas: %d\n", count_zero);

		free(g);
		free(x);
	}
}

int main() {
	simulate(10, 10, 0.5);

	return 0;
}