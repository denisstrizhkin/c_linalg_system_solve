#define _POSIX_C_SOURCE 200809L

#include "mat.h"

#include <float.h>
#include <stdio.h>

int
main(void)
{
	unsigned n = 3;

	mat *T1, *T2, *T3, *T4, *T5, *A;
	
	T1 = mat_rnd(n, n, 100);
	mat_scale(T1, 5);
	T2 = mat_T(T1);

	A = mat_dot(T1, T2);
	mat_free(T1);
	mat_free(T2);


	double eps = DBL_MAX;
	double eps_adm = 0.001;

	double norm;

	double eval;

	mat *x_prev = mat_rnd(n, 1, 100);
	mat *x_cur  = mat_rnd(n, 1, 100);

	unsigned i = 0;
	while (eps > eps_adm) {
		mat_free(x_prev);
		x_prev = mat_copy(x_cur);
		mat_free(x_cur);

		norm = mat_norm1(x_prev);
		x_cur = mat_dot(A, x_prev);
		mat_scale(x_cur, 1.0/norm);

		T1 = mat_scale_o(x_cur, norm);
		T2 = mat_T(T1);
		mat_free(T1);
		T1 = mat_T(x_prev);
		T3 = mat_dot(T2, x_prev);
		T4 = mat_dot(T1, x_prev);

		eval = mat_get(T3, 0, 0) / mat_get(T4, 0, 0);

		mat_free(T1);
		mat_free(T2);
		mat_free(T3);
		mat_free(T4);

		T1 = mat_dot(A, x_prev);
		T2 = mat_scale_o(x_prev, eval);
		T3 = mat_sub_o(T1, T2);

		mat_free(T1);
		mat_free(T2);
		mat_free(T3);

		eps = mat_norm1(T3);
		i++;
	}

	mat_print(A);
	mat_print(x_prev);
	printf("iteration: %u, eval: %lf\n\n", i, eval);

	mat_free(A);

	return 0;
}
