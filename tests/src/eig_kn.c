#include "../../src/graphs.h"
#include "../../src/eig.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
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

void simulate(size_t maxit) {
	assert(maxit >= 3);
	double rtol = 1e-9;
	double atol = 1e-12;

	for (size_t i = 3; i <= maxit; i++) {
		Graph* g = graph_kn(i);
		double* w = malloc(i * sizeof(double));

		graph_spec_adj(g, w);

		printf("\n==== n = %ld ====\n", i);
		printf("λ = %f, m(λ) = %d\n", w[i-1], count_mult(w,
			(double) i - 1, i, rtol, atol));
		printf("λ = %f, m(λ) = %d", w[0], count_mult(w, 
			(double) -1, i, rtol, atol));
		printf("\n");

		free(w);
		free(g);
	}
}

int main() {
	simulate(20);

	return 0;
}