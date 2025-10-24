#include <stdio.h>
#include <lapacke.h>
#include "eig.h"

static double* matrix_cpy(const double* A, size_t n) {
	double* A_cpy = malloc(n * n * sizeof(double));

	if (!A_cpy) {
		die("malloc error (A_cpy)");
	}

	for (size_t i = 0; i < n * n; i++) {
		A_cpy[i] = A[i];
	}

	return A_cpy;
}

static double trace(const double* A, size_t n) {
	double t = 0.0;

	for (size_t i = 0; i < n; i++) {
		t += A[IDX(i, i, n)];
	}

	return t;
}


/*Usa lapacke e cblas para achar o espectro de uma matriz.
'N' indica que não queremos achar os autovetores de A.
Caso seja de seu desejo (e caso você esteja com paciência)
achar os autovetores, altere a função (isso mesmo, você altera)
trocando 'N' para 'V' e retornando a matriz A_cpy, pois ela
conterá os autovetores (ou troque A_cpy por A e tire o const
do argumento double* A, mas aí você vai perder o conteúdo da
matriz passada)*/
int matrix_spec(const double* A, size_t n, double* x) {
	double* A_cpy = matrix_cpy(A, n);

	int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'N', 'U', (int) n, 
		A_cpy, (int) n, x);

	free(A_cpy);
	return info;
}


int graph_spec_adj(const Graph* g, double* x) {
	double* A_cpy = matrix_cpy(g->A, g->n);

	int info = matrix_spec(A_cpy, g->n, x);

	free(A_cpy);
	return info;
}


int graph_spec_lap(const Graph* g, double* x) {
	double* l = (double*) calloc(g->n * g->n, sizeof(double));

	int info = matrix_spec(l, g->n, x);

	free(l);
	return info;
}

void eigenvalues_print(const double* x, size_t n, const char* title) {
	if (title) {
		printf("%s\n", title);
	}

	for (size_t i = 0; i < n; i++) {
		printf("λ[%zu] = %.6f\n", i, x[i]);
	}	
}


/*Usa as fórmulas/algoritmo de Faddeev–LeVerrier para achar 
os coeficientes:
https://en.wikipedia.org/wiki/Faddeev%E2%80%93LeVerrier_algorithm
*/
void matrix_char_coeffs(const double* A, size_t n, double* coeffs) {
	double* Ak = matrix_cpy(A, n);
	double* Tmp = calloc(n * n, sizeof(double));

	if (!Ak || !Tmp) {
		die("malloc error (Ak || Tmp)");
	}

	coeffs[0] = 1.0;

	double* S = calloc(n + 1, sizeof(double));

	for (size_t k = 1; k <= n; k++) {
		S[k] = trace(Ak, n);

		if (k < n) {
			matrix_mult(A, Ak, Tmp, n);
			double* swap = Ak;
			Ak = Tmp;
			Tmp = swap;
		}
	}

	for (size_t k = 1; k <= n; k++) {
		double sum = 0.0;

		for (size_t i = 1; i <= k; i++) {
			sum += coeffs[k - i] * S[i];
		}

		coeffs[k] = -sum / (double) k;
	}

	free(Ak);
	free(Tmp);
	free(S);
}