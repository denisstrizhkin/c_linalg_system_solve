#define _POSIX_C_SOURCE 200809L

#include "mat.h"

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

	return A0;
}

int
main(void)
{
	unsigned n = 3;

	mat *A0 = get_A0(n);

	mat *b0 = mat_rnd(n, 1, 100);
	mat_scale(b0, 100);

	double t = mat_norm1(A0);

	mat_scale(A0, 1/t);
	mat_scale(b0, 1/t);

	mat_print(A0);
	mat_print(b0);

	free(A0);
	free(b0);
	return 0;
}
