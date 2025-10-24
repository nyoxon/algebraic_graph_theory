#include "graphs.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <limits.h>

static inline void bounds_check(size_t n, size_t u, size_t v) {
	if (u >= n || v >= n) {
		die("vértice fora do intervalo");
	}
}

void die(const char* msg) {
	fprintf(stderr, "erro: %s\n", msg);
	exit(EXIT_FAILURE);
}

Graph* graph_new(size_t	n, bool directed) {
	Graph* g = (Graph*) calloc(1, sizeof(Graph));

	if (!g) {
		die("malloc error (Graph)");
	}

	g->n = n;
	g->directed = directed;
	g->A = (double*) calloc(n*n, sizeof(double));

	if (!g->A) {
		free(g);
		die("malloc error (A)");
	}

	return g;
}

void graph_free(Graph* g) {
	if (!g) {
		return;
	}

	free(g->A);
	free(g);
}

void graph_clear(Graph* g) {
	if (!g) {
		return;
	}

	memset(g->A, 0, g->n * g->n * sizeof(double));
}

void graph_add_edge(Graph* g, size_t u, size_t v, double w) {
	bounds_check(g->n, u, v);
	g->A[IDX(u, v, g->n)] = w;

	/*Se g não for direcionado, g->A será simétrica*/
	if (!g->directed && u != v) {
		g->A[IDX(v, u, g->n)] = w;
	}
}

void graph_remove_edge(Graph* g, size_t u, size_t v) {
	bounds_check(g->n, u, v);
	g->A[IDX(u, v, g->n)] = 0.0;

	if (!g->directed && u != v) {
		g->A[IDX(v, u, g->n)] = 0.0;
	}
}

double graph_get(const Graph* g, size_t u, size_t v) {
	bounds_check(g->n, u, v);

	return g->A[IDX(u, v, g->n)];
}

size_t graph_num_vertices(const Graph* g) {
	return g->n;
}

bool graph_is_directed(const Graph* g) {
	return g->directed;
}

void graph_print(const Graph* g, const char* title) {
	if (title) {
		printf("%s\n", title);
	}

	size_t n = g->n;

	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			printf("%g%c", g->A[IDX(i, j, n)], (j+1<n?' ':'\n'));
		}
	}
}

void matrix_print(const double* A, size_t n, const char* name) {
	if (name) {
		printf("%s:\n", name);
	}

	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			printf("%6.2f ", A[i * n + j]);
		}

		printf("\n");
	}
}

/*
Se g é não direcionado, então deg_out = deg_in = soma da linha.
Caso contrário, deg_out = soma da linha e deg_in = soma da coluna.*/
void graph_degree(const Graph* g, double* deg_out, double* deg_in) {
	size_t n =  g->n;

	if (deg_out) {
		memset(deg_out, 0, n * sizeof(double));
	}

	if (deg_in) {
		memset(deg_in, 0, n * sizeof(double));
	}

	if (!g->directed) {
		for (size_t i = 0; i < n; i++) {
			double s = 0.0;
			const double* row = &g->A[IDX(i, 0, n)];

			for (size_t j = 0; j < n; j++) {
				s += row[j];
			}

			if (deg_out) {
				deg_out[i] = s;
			}

			if (deg_in) {
				deg_in[i] = s;
			}
		}
	} else {
		if (deg_out) {
			for (size_t i = 0; i < n; i++) {
				double s = 0.0;
				const double* row = &g->A[IDX(i, 0, n)];

				for (size_t j = 0; j < n; j++) {
					s += row[j];
				}

				deg_out[i] = s;
			}
		}

		if (deg_in) {
			for (size_t j = 0; j < n; j++) {
				double s = 0.0;

				for (size_t i = 0; i < n; i++) {
					s += g->A[IDX(i, j, n)];
				}

				deg_in[j] = s;
			}
		}
	}
}

void graph_degree_matrix(const Graph* g, double* D) {
	size_t n = g->n;

	memset(D, 0, n * n * sizeof(double));
	double* deg = (double*) calloc(n, sizeof(double));

	if (!deg) {
		die("malloc error (deg)");
	}

	graph_degree(g, deg, NULL);

	for (size_t i = 0; i < n; i++) {
		D[IDX(i, i, n)] = deg[i];
	}

	free(deg);
}

void graph_laplacian(const Graph* g, double* L) {
	size_t n = g->n;

	for (size_t i = 0; i < n * n; i++) {
		L[i] = -g->A[i];
	}

	double* deg = (double*) calloc(n, sizeof(double));

	if (!deg) {
		die("malloc error (deg)");
	}

	graph_degree(g, deg, NULL);

	for (size_t i = 0; i < n; i++) {
		L[IDX(i, i, n)] += deg[i];
	}

	free(deg);
}

void graph_normalized_laplacian(const Graph* g, double* Ln) {
	size_t n = g->n;
	double* deg = (double*) calloc(n, sizeof(double));

	if (!deg) {
		die("malloc error (deg)");
	}

	graph_degree(g, deg, NULL);
	memset(Ln, 0, n * n * sizeof(double));

	for (size_t i = 0; i < n; i++) {
		Ln[IDX(i, i, n)] = 1.0;
	}

	for (size_t i = 0; i < n; i++) {
		double di = deg[i];
		double dinv_i = (di > 0.0) ? 1.0/sqrt(di) : 0.0;

		for (size_t j = 0; j < n; j++) {
			double dj = deg[j];
			double dinv_j = (dj > 0.0) ? 1.0/sqrt(dj) : 0.0;
			double a = g->A[IDX(i, j, n)];

			Ln[IDX(i, j, n)] -= dinv_i * a * dinv_j;
		}
	}

	free(deg);
}

void matrix_vecmult(const double* A, size_t n, const double* x, double* y) {
	for (size_t i = 0; i < n; i++) {
		const double* row = &A[IDX(i, 0, n)];
		double s = 0.0;

		for (size_t j = 0; j < n; j++) {
			s += row[j] * x[j];
		}

		y[i] = s;
	}
}

void matrix_mult(const double* A, const double* B, double* C, size_t n) {
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			double s = 0.0;

			for (size_t k = 0; k < n; k++) {
				s += A[IDX(i, k, n)] * B[IDX(k, j, n)];
			}

			C[IDX(i, j, n)] = s;
		}
	}
}

void graph_ax(const Graph* g, const double* x, double* y) {
	matrix_vecmult(g->A, g->n, x, y);
}

Graph* graph_read_from_file(const char* path) {
	FILE* f = fopen(path, "r");

	if (!f) {
		die("não consegui abrir o arquivo");
	}

	size_t n, m;
	int dir;

	if (fscanf(f, "%zu %d", &n, &dir) != 2) {
		die("cabeçalho inválido");
	}

	if (fscanf(f, "%zu", &m) != 1) {
		die("m inválido");
	}

	Graph* g = graph_new(n, dir != 0);

	for (size_t k = 0; k < m; k++) {
		size_t u, v;
		double w;

		if (fscanf(f, "%zu %zu %lf", &u, &v, &w) != 3) {
			die("linha de aresta inválida");
		}

		graph_add_edge(g, u, v, w);
	}

	fclose(f);
	return g;
}

double graph_paths_length
(
	const Graph* g, 
	size_t v, 
	size_t u, 
	unsigned int k
)
{
	size_t n = g->n;

	if (k == 0) {
		return (v == u) ? 1.0 : 0.0; /*A⁰ = I*/;
	}

	if (k == 1) {
		return g->A[IDX(v, u, n)];
	}

	double sum = 0.0;

	for (size_t w = 0; w < n; w++) {
		if (g->A[IDX(v, w, n)] != 0.0) {
			sum += g->A[IDX(v, w, n)] * graph_paths_length(g,
				w, u, k - 1);
		}
	}

	return sum;
}

Graph* graph_kn(size_t n) {
	Graph* Kn = graph_new(n, false);

	/*Trabalha apenas na parte triangular superior*/
	for (size_t i = 0; i < n; i++) {
		for (size_t j = i + 1; j < n; j++) {
			graph_add_edge(Kn, i, j, 1.0);
		}
	}

	return Kn;
}

/*BFS simples que retorna a maior distância a partir um vértice src
de um grafo CONEXO. Não verifico se o grafo é conexo ou não,
então é melhor passar um grafo conexo para a função :D*/
static int graph_bfs_longest_from(const Graph* g, size_t src) {
	size_t n = g->n;

	int* dist = malloc(n * sizeof(int));
	bool* visited = calloc(n, sizeof(bool));

	if (!dist || !visited) {
		die("malloc error (dist || visited");
	}

	for (size_t i = 0; i < n; i++) {
		dist[i] = INT_MAX;
	}

	dist[src] = 0;
	visited[src] = true;

	size_t* queue = malloc(n * sizeof(size_t));
	size_t front = 0, back = 0;
	queue[back++] = src;

	while (front < back) {
		size_t u = queue[front++];

		for (size_t v = 0; v < n; v++) {
			if (g->A[IDX(u, v, n)] != 0.0 && !visited[v]) {
				visited[v] = true;
				dist[v] = dist[u] + 1;
				queue[back++] = v;
			}
		}
	}

	int max_d = 0;

	for (size_t i = 0; i < n; i++) {
		if (dist[i] != INT_MAX && dist[i] > max_d) {
			max_d = dist[i];
		}
	}

	free(queue);
	free(visited);
	free(dist);

	return max_d;
}

int graph_diameter(const Graph *g) {
	int diameter = 0;

	for (size_t i = 0; i < g->n; i++) {
		int d = graph_bfs_longest_from(g, i);

		if (d > diameter) {
			diameter = d;
		}
	}

	return diameter;
}