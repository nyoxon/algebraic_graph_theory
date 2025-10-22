#ifndef GRAPHS_H
#define GRAPHS_H

typedef struct {
  int n; // |V|, número de vértices
  int** adj; // A, matriz de adjacência
} Graph;


/* Retorna um ponteiro para um grafo com
   matriz preenchida por zeros. */
Graph* graph_new(int n);

double** int_to_double_matrix(int **M, int n);
void free_dmatrix(double **M, int n);
void free_matrix(int** mat, int n);
void graph_free(Graph* g);

/* Adiciona uma aresta com extremos vi e vj no grafo.
   Se o grafo não for direcionado, então a matriz de adjacência
   será simétrica. */
void graph_add_edge(Graph* g, int i, int j, int directed);

/* Cria um grafo lendo arestas da entrada padrão. */
Graph* read_graph(int directed);


/* Printa o conteúdo de um grafo */
void graph_print(int** adj, int n);


/* Constrói a matriz laplaciana do grafo dado. */
int** graph_laplacian(Graph *g);

#endif
