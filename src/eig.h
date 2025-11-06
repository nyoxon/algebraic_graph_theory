#ifndef EIG_H
#define EIG_H

/* --- Declarações de funções "espectrais" usadas neste projeto. --- */

/*
### NOTA IMPORTANTE ###

Leia a nota importante em graphs.h :D
*/


#include "graphs.h"


/*Acha o espectro (conjunto de autovalores) de A e
retorna um int indicando erro se ele (o int) for > 0*/
int matrix_spec(const double* A, size_t n, double* x);


/*Acha o espectro da matriz de adjacência de g e
retorna um int indicando erro se ele (o int) for > 0*/
int graph_spec_adj(const Graph* g, double* x);


/*Acha o espectro da matriz laplaciana de g e
retorna um int indicando erro se ele (o int) for > 0*/
int graph_spec_lap(const Graph* g, double* x);


/*Printa bonitinho o conteúdo de um vetor (double* x passado
à matrix_spec/graph_spec*) contendo o espectro de uma matriz :D*/
void eigenvalues_print(const double* x, size_t n, const char* title);


/*Acha os coeficientes do polinômio característico de A.
O algoritmo implementado tem complexidade O(n⁴), então use
esta função com moderação.
Eu poderia ter implementado de forma a achar apenas os k primeiros
coeficientes para o algoritmo ter complexidade até O(kn³), mas preferi
fazer assim mesmo*/
void matrix_char_coeffs(const double* A, size_t n, double* coeffs);

/*Acha o rank de uma matriz*/
unsigned int matrix_rank(const double* A, size_t n, size_t m, double tol);

#endif