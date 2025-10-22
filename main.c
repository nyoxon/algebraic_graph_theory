#include <lapacke.h>
#include <stdio.h>

int main(void) {
    int n = 4;
    double A[16] = {
        1,-1, 0, 0,
       -1, 2,-1, 0,
        0,-1, 2,-1,
        0, 0,-1, 1
    }; // Laplaciana do P4 (armazenada por linha)
    double w[4]; // autovalores

    // dsyev: eigenvalues of real symmetric matrix
    // 'V' = compute eigenvectors, 'U' = upper triangle
    int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n, A, n, w);

    if (info > 0) {
        fprintf(stderr, "Decomposição falhou.\n");
        return 1;
    }

    printf("Autovalores:\n");
    for (int i = 0; i < n; i++)
        printf("%.10f\n", w[i]);

    return 0;
}
