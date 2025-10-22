#include "graphs.h"
#include <stdlib.h>
#include <stdio.h>

static int** create_matrix(int n) {
  int** mat = malloc(n * sizeof(int*));

  for (int i = 0; i < n; i++) {
    mat[i] = calloc(n, sizeof(int));
  }

  return mat;
}

double** int_to_double_matrix(int** M, int n) {
  double** R = malloc(n * sizeof(double*));

  if (!R) {
    return NULL;
  }

  for (int i = 0; i < n; i++) {
    R[i] = malloc(n * sizeof(double));

    if (!R[i]) {
      for (int k = 0; k < i; k++) {
        free(R[k]);
      }

      free(R);
      return NULL;
    }

    for (int j = 0; j < n; j++) {
      R[i][j] = (double)M[i][j];
    }
  }

  return R;
}

void free_dmatrix(double **M, int n) {
    for (int i = 0; i < n; i++)
        free(M[i]);
    free(M);
}

void free_matrix(int** mat, int n) {
  for (int i = 0; i < n; i++) {
    free(mat[i]);
  }

  free(mat);
}


Graph* graph_new(int n) {
  Graph* g = malloc(sizeof(Graph));
  g->n = n;
  g->adj = create_matrix(n);

  return g;
}

void graph_free(Graph* g) {
  free_matrix(g->adj, g->n);
  free(g);
}


void graph_add_edge(Graph* G, int i, int j, int directed) {
  if (i >= G->n || j >= G->n) {
    return;
  }

  G->adj[i][j] = 1;

  if (!directed) {
    G->adj[j][i] = 1;
  }
}

Graph* read_graph(int directed) {
  int n, m;
  printf("'número de vértices' 'número de arestas'\n");
  scanf("%d %d", &n, &m);

  Graph* g = graph_new(n);

  for (int k = 0; k < m; k++) {
    int u, v;
    scanf("%d %d", &u, &v);
    graph_add_edge(g, u, v, directed);
  }

  return g;
}

void graph_print(int** adj, int n) {
  printf("%d vértices\n", n);
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      printf("%d ", adj[i][j]);
    }

    printf("\n");
  }
}

int** graph_laplacian(Graph* g) {
  int n = g->n;
  int** L = create_matrix(n);

  for (int i = 0; i < n; i++) {
    int deg = 0;

    for (int j = 0; j < n; j++) {
      if (g->adj[i][j]) {
        deg++;
        L[i][j] = -1;
      }
    }

    L[i][i] = deg;
  }

  return L;
}
