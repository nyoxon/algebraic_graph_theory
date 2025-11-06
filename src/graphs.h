#ifndef GRAPHS_H
#define GRAPHS_H

/* --- Declarações de funções "não-espectrais" usadas neste projeto. ---*/
/* Uma função espectral é uma função que se relaciona com o espectro
de uma matriz de alguma forma. As declarações desse tipo estão
em eig.h*/

/*
### NOTA IMPORTANTE ###

Para usar funções que criam objetos mas que são do tipo void,
como graph_degree, graph_degree_matrix, graph_laplacian etc,
você deve primeiro manualmente alocar espaço suficiente para o
objeto a ser criado e passar um ponteiro VÁLIDO para esse espaço
para a função.

Exemplo com graph_degree:
			  malloc(g->n * sizeof(double));
double* deg = calloc(g->n, sizeof(double)); // #include <stdlib.h> para calloc

if (!deg) { // garante que o ponteiro é válido. opcional, mas recomendado
	perror("calloc"); // #include <stdio.h> para perror
	exit(1);
}

graph_degree(g, deg, NULL);

--- faz trabalho aqui ---

free(deg); // caso o programa não saia imediatamente
*/

#include <stdbool.h>
#include <stdlib.h>

/*
Toda matriz considerada em qualquer código deste projeto
é quadrada, então as libs foram implementadas pressupondo
que qualquer matriz passada como argumento será quadrada.
A matriz de adjacência é implementada como uma
row-major order matrix:
https://en.wikipedia.org/wiki/Row-_and_column-major_order

O elemento (i, j) é acessado usando a fórmula i * n + j onde 
n é o tamanho da matriz.*/
typedef struct {
	size_t n;  			/* N° de vértices*/
	bool directed;		/* Grafo direcionado?*/
	double* A;			/* Matriz de adjacência*/
} Graph;



/*Macro para acessar o elemento (i, j)*/
#define IDX(i, j, n) ((i) * (n) + j)


/*Função auxiliar para erros*/
void die(const char* msg);


/*Cria um grafo em memória dinâmica (free-after-use) e 
retorna um ponteiro para ele.
A matriz g->A é preenchida com zeros, então o acesso à memória
é seguro mesmo se você não adicionar nenhuma aresta logo
após ter chamado a função.*/
Graph* graph_new(size_t n, bool directed);


/*Cria um grafo aleatório e retorna um ponteiro para ele*/
Graph* graph_random(size_t n, double p) ;


/*Cria um grafo conexo aleatório*/
Graph* graph_random_connected(int n, double p);


/* Cria um grafo regular aleatório */
Graph* graph_random_regular(size_t n, size_t k);


/*Cria um grafo bipartido (de partições de tamanho n1 e n2) aleatório*/
Graph* graph_random_bipartite(size_t n1, size_t n2, double p);


/*Verifica se g é conexo*/
bool graph_is_connected(const Graph* g);


/*Retorna o número de arestas de um grafo qualquer*/
size_t graph_num_edges(const Graph* g);


/*Libera o conteúdo de um grafo*/
void graph_free(Graph* g);


/*Preenche g->A com zeros (ao usar graph_new(), a matriz já
é preenchida com zeros)*/
void graph_clear(Graph* g);


/*Adiciona a aresta no grafo colocando o valor w nas posições
correspondentes de g->A.
Se o grafo for não direcionado, é de se esperar que w = 1.0;
*/
void graph_add_edge(Graph *g, size_t u, size_t v, double w);


/*Remove a aresta uv do grafo*/
void graph_remove_edge(Graph* g, size_t u, size_t v);


/*Retorna g->A[u, v]. Equivalente a g->A[IDX(u, v, g->n)], mas
com bound check para u e v*/
double graph_get(const Graph* g, size_t u, size_t v);


/*Retorna o número de vértices do grafo*/
size_t graph_num_vertice(const Graph* g);


/*Verifica se o grafo é direcionado*/
bool graph_is_directed(const Graph* g);


/*Printa g->A*/
void graph_print(const Graph* g, const char* title);


/*Printa uma matriz qualquer*/
void matrix_print(const double* A, size_t n, const char* name);


/*Coloca nos vetores deg_out e deg_in o grau (out/in) de cada vértice
se directed = 0, então deg_out = deg_in*/
void graph_degree(const Graph* g, double* deg_out, double* deg_in);


/*Calcula D = matriz de graus*/
void graph_degree_matrix(const Graph* g, double* d);


/*Calcula L = D - A laplaciana*/
void graph_laplacian(const Graph* g, double *l);


/*Calcula a matriz de incidência de um grafo
Diferentemente das duas últimas funções, esta retorna
diretamente a matriz. Isso porque precisamos calcular
o número de arestas de um grafo, o que achei mais
conveniente a função fazer.*/
double* graph_incidence_matrix(const Graph* g);


/*Calcula Ln = (D)^1/2 * L * (D)^1/2 laplaciana normalizada*/
void graph_normalized_laplacian(const Graph* g, double *ln);


/*y = Ax*/
void matrix_vecmult(const double* A, size_t n, const double* x, double* y);


/*C = AB*/
void matrix_mult(const double* A, const double* B, double* C, size_t n);


/*y = g->A * x*/
void graph_ax(const Graph* g, const double* x, double* y);


/*Cria um grafo lendo informações de uma arquivo.

O formato do arquivo deve ser da seguinte forma:

n dir
m
u1 v1 w1
u2 v2 w1
.
.
.
um vm wm

n = numero de vértices
dir = direcionado
ui, vi, wi = aresta com extremidadades ui e vi e valor wi

Veja o exemplo (example_file_graph_format) no diretório para 
entender melhor.
*/
Graph* graph_read_from_file(const char* path);


/*Acha o número de caminhos de tamanho k entre u e v*/
double graph_paths_length(const Graph* g, size_t v, size_t u, unsigned int k);


/*Cria o Kn*/
Graph* graph_kn(size_t n);

/*Acha o diâmetro do grafo*/
int graph_diameter(const Graph *g);


#endif