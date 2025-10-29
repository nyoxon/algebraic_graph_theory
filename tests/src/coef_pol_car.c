#include "../../src/graphs.h"
#include "../../src/eig.h"
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <assert.h>

unsigned long long binomial(unsigned int n, unsigned int k) {
    if (k > n) return 0;
    if (k == 0 || k == n) return 1;

    if (k > n - k)
        k = n - k;

    unsigned long long res = 1;
    for (unsigned int i = 1; i <= k; i++) {
        res = res * (n - i + 1) / i;
    }

    return res;
}

long long graph_count_triangles_naive(const Graph *g) {
    size_t n = g->n;
    long long cnt = 0;

    if (g->directed) {
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j) if (i!=j && g->A[i*n+j])
                for (size_t k = 0; k < n; ++k) if (k!=i && k!=j && g->A[j*n+k] && g->A[k*n+i])
                    cnt++;
        return cnt;
    } else {
        for (size_t i = 0; i+2 < n; ++i)
            for (size_t j = i+1; j+1 < n; ++j) if (g->A[i*n+j])
                for (size_t k = j+1; k < n; ++k)
                    if (g->A[i*n+k] && g->A[j*n+k]) cnt++;
        return cnt;
    }
}

void simulate(size_t n, size_t n_times, double p) {
	srand(time(NULL));

	for (size_t i = 0; i < n_times; i++) {
		Graph* g = graph_random(n, 0, p);
		double* w = malloc(g->n * sizeof(double));
		matrix_char_coeffs(g->A, g->n, w);

		assert(w[1] == 0.0);
		assert(-w[2] == (double)graph_num_edges(g));
		assert(-w[3] / 2.0 == (double)graph_count_triangles_naive(g));

		free(w);
	}
}

int main() {
	simulate(20, 1000, 0.5);

	return 0;
}