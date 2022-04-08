#define _POSIX_C_SOURCE 200809L

#include "mat.h"

#include <float.h>
#include <stdio.h>


mat *
get_A0(unsigned n)
{
	mat *T1 = mat_ltrig_rnd(n, 100);
	mat_scale(T1, 100);

	mat *T2 = mat_I(n);
	mat_scale(T2, 100 * n);

	mat *T3 = mat_rnd(n, n, 15);
	mat_scale(T3, 100);

	mat *T4 = mat_add_o(T3, T2);
	mat *A0 = mat_add_o(T4, T1);

	mat_free(T1);
	mat_free(T2);
	mat_free(T3);
	mat_free(T4);
	/*
	mat *A0 = mat_new(n,n);
	mat_set(A0, 0, 0, 13.0);
	mat_set(A0, 0, 1, 0.0);
	mat_set(A0, 0, 2, 3.0);
	mat_set(A0, 1, 0, 7.0);
	mat_set(A0, 1, 1, 15.0);
	mat_set(A0, 1, 2, 0.0);
	mat_set(A0, 2, 0, 4.0);
	mat_set(A0, 2, 1, 9.0);
	mat_set(A0, 2, 2, 11.0);
*/
	return A0;
}

mat *
get_b0(unsigned n)
{
	mat *b0 = mat_rnd(n, 1, 100);
	mat_scale(b0, 100);
/*
	mat *b0 = mat_new(n,1);
	mat_set(b0, 0, 0, 5.0);
	mat_set(b0, 1, 0, 2.0);
	mat_set(b0, 2, 0, 6.0);*/
	return b0;
}

mat *
get_S(unsigned n, const mat *A)
{
	mat *T = mat_I(n);
	mat *S = mat_sub_o(T, A);

	mat_free(T);

	return S;
}

int
main(void)
{
	unsigned n = 10;

	mat *A0 = get_A0(n);
	mat_print(A0);
	mat *b0 = get_b0(n);
	mat_print(b0);

	double t = mat_norm1(A0);
	printf("t: %lf\n", t);
	mat *A = mat_scale_o(A0, 1/t);
	mat *b = mat_scale_o(b0, 1/t);
	mat_print(A);
	mat_print(b);

	mat *S = get_S(n, A);
	mat_print(S);
	double normS = mat_norm1(S);
	printf("normS: %lf\n", normS);
	if (normS >= 1)
		exit(1);

	mat *x_prev, *x_cur;
	double eps = DBL_MAX;
	double eps_adm = 0.003;

	x_cur = mat_zeros(n, 1);
	x_prev = mat_zeros(n, 1);

	unsigned j = 0;
	while (eps > eps_adm) {
		mat_free(x_prev);
		x_prev = mat_copy(x_cur);
		
		mat *T1 = mat_dot(S, x_prev);
		x_cur = mat_add_o(T1, b);

		mat *T2 = mat_sub_o(x_cur, x_prev);
		eps = normS / (1 - normS ) * mat_norm1(T2);

		mat_free(T1);
		mat_free(T2);

		j++;
	}
	printf("iteration: %u\n", j);
	mat_print(x_cur);

	mat_free(x_prev);
	mat_free(x_cur);
	mat_free(S);
	mat_free(A0);
	mat_free(b0);	
	mat_free(A);
	mat_free(b);

	return 0;
}
