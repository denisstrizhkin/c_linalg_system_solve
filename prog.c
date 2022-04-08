#define _POSIX_C_SOURCE 200809L

#include "mat.h"

#include <float.h>
#include <stdio.h>


mat *
get_A0_part1(unsigned n)
{
	mat *T1 = mat_rnd(n, n, 100);
	mat_scale(T1, 10);

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
get_A0_part2(unsigned n)
{
	mat *T1 = mat_D_rnd(n, 100);
	mat_scale(T1, 100 * n);

	mat *T2 = mat_rnd(n, n, 10);
	mat_scale(T2, 10);

	mat *A0 = mat_add_o(T1, T2);

	mat_free(T1);
	mat_free(T2);

	return A0;
}

mat *
get_A0_part3(unsigned n)
{
	mat *T1 = mat_ltrig_rnd(n, 100);
	mat_scale(T1, 10);

	mat *T2 = mat_I(n);
	mat_scale(T2, 100 * n);

	mat *T3 = mat_rnd(n, n, 20);
	mat_scale(T3, 10);

	mat *T4 = mat_add_o(T3, T2);
	mat *A0 = mat_add_o(T4, T1);

	mat_free(T1);
	mat_free(T2);
	mat_free(T3);
	mat_free(T4);

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
pi_get_S(unsigned n, const mat *A)
{
	mat *T = mat_I(n);
	mat *S = mat_sub_o(T, A);

	mat_free(T);

	return S;
}

void
plain_iteration(const mat *A0, const mat *b0, unsigned n, double eps_adm)
{
	puts("\n### Plain Iteration method ###");

	double t = mat_norm1(A0);
	mat *A = mat_scale_o(A0, 1/t);
	mat *b = mat_scale_o(b0, 1/t);

	mat *S = pi_get_S(n, A);

	double normS = mat_norm1(S);
	printf("\n === normS: %lf\n", normS);
	if (normS >= 1)
	{
		printf(" divergent\n");
		return;
	}

	mat *x_prev, *x_cur;
	double eps = DBL_MAX;

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
	printf(" === iteration: %u\n", j);
	printf("\n === x: ");
	mat_print(x_cur);

	mat_free(x_prev);
	mat_free(x_cur);
	mat_free(S);
	mat_free(A);
	mat_free(b);	
}

mat *
j_get_S(const mat *L, const mat *U, const mat *Dr)
{
	mat *T1 = mat_add_o(L, U);
	mat *T2 = mat_scale_o(Dr, -1.0);
	mat *S = mat_dot(T2, T1);

	mat_free(T1);
	mat_free(T2);

	return S;
}

void
jakobi(const mat *A0, const mat *b0, unsigned n, double eps_adm)
{
	puts("\n### Jakobi method ###");

	double t = mat_norm1(A0);
	mat *A = mat_scale_o(A0, 1/t);
	mat *b = mat_scale_o(b0, 1/t);

	mat *L = mat_get_L(A);
	mat *U = mat_get_U(A);
	mat *Dr = mat_get_Dr(A);

	mat *S = j_get_S(L, U, Dr);
	mat *c = mat_dot(Dr, b);

	double normS = mat_norm1(S);
	printf("\n === normS: %lf\n", normS);
	if (normS >= 1)
	{
		printf(" divergent\n");
		return;
	}

	mat *x_prev, *x_cur;
	double eps = DBL_MAX;

	x_cur = mat_zeros(n, 1);
	x_prev = mat_zeros(n, 1);

	unsigned j = 0;
	while (eps > eps_adm) {
		mat_free(x_prev);
		x_prev = mat_copy(x_cur);
		
		mat *T1 = mat_dot(S, x_prev);
		x_cur = mat_add_o(T1, c);

		mat *T2 = mat_sub_o(x_cur, x_prev);
		eps = normS / (1 - normS ) * mat_norm1(T2);

		mat_free(T1);
		mat_free(T2);

		j++;
	}
	printf(" === iteration: %u\n", j);
	printf("\n === x: ");
	mat_print(x_cur);

	mat_free(x_prev);
	mat_free(x_cur);
	mat_free(S);
	mat_free(A);
	mat_free(b);	
	mat_free(c);
	mat_free(L);
	mat_free(U);
	mat_free(Dr);
}

mat *
z_get_S(const mat *L, const mat *U, const mat *D)
{
	mat *T1 = mat_add_o(L, D);
	mat_ltrig_rev(T1);
	mat_scale(T1, -1.0);

	mat *S = mat_dot(T1, U);

	mat_free(T1);

	return S;
}

mat *
z_get_c(const mat *L, const mat *D, const mat *b)
{
	mat *T1 = mat_add_o(L, D);
	mat_ltrig_rev(T1);

	mat *c = mat_dot(T1, b);

	mat_free(T1);

	return c;
}

void
zeidel(const mat *A0, const mat *b0, unsigned n, double eps_adm)
{
	puts("\n### Zeidel' method ###");

	double t = mat_norm1(A0);
	mat *A = mat_scale_o(A0, 1/t);
	mat *b = mat_scale_o(b0, 1/t);

	mat *L = mat_get_L(A);
	mat *U = mat_get_U(A);
	mat *D = mat_get_D(A);

	mat *S = z_get_S(L, U, D);
	mat *c = z_get_c(L, D, b);

	double normS = mat_norm1(S);
	printf("\n === normS: %lf\n", normS);
	if (normS >= 1)
	{
		printf(" divergent\n");
		return;
	}

	mat *x_prev, *x_cur;
	double eps = DBL_MAX;

	x_cur = mat_zeros(n, 1);
	x_prev = mat_zeros(n, 1);

	unsigned j = 0;
	while (eps > eps_adm) {
		mat_free(x_prev);
		x_prev = mat_copy(x_cur);
		
		mat *T1 = mat_dot(S, x_prev);
		x_cur = mat_add_o(T1, c);

		mat *T2 = mat_sub_o(x_cur, x_prev);
		eps = normS / (1 - normS ) * mat_norm1(T2);

		mat_free(T1);
		mat_free(T2);

		j++;
	}
	printf(" === iteration: %u\n", j);
	printf("\n === x: ");
	mat_print(x_cur);

	mat_free(x_prev);
	mat_free(x_cur);
	mat_free(S);
	mat_free(A);
	mat_free(b);	
	mat_free(c);
	mat_free(L);
	mat_free(U);
	mat_free(D);
}

int
main(void)
{
	/* Задание №1 */
	unsigned n = 9;
	double eps = 0.00003;
	
	mat *A0 = get_A0_part1(n);
	mat *b0 = get_b0(n);

	printf("\n### PART 1 | n: %u, eps_adm: %lf ###", n, eps);
	printf("\n === A: ");
	mat_print(A0);
	printf("\n === b: ");
	mat_print(b0);

	plain_iteration(A0, b0, n, eps);
	jakobi(A0, b0, n, eps);
	zeidel(A0, b0, n, eps);

	mat_free(A0);
	mat_free(b0);
	/* ====== */

	/* Задание №2 */
	n = 8;
	eps = 0.00003;
	
	A0 = get_A0_part2(n);
	b0 = get_b0(n);

	printf("\n### PART 2 | n: %u, eps_adm: %lf ###", n, eps);
	printf("\n === A: ");
	mat_print(A0);
	printf("\n === b: ");
	mat_print(b0);

	plain_iteration(A0, b0, n, eps);
	jakobi(A0, b0, n, eps);
	zeidel(A0, b0, n, eps);

	mat_free(A0);
	mat_free(b0);
	/* ====== */

	/* Задание №3 */
	n = 8;
	eps = 0.00005;
	
	A0 = get_A0_part3(n);
	b0 = get_b0(n);

	printf("\n### PART 3 | n: %u, eps_adm: %lf ###", n, eps);
	printf("\n === A: ");
	mat_print(A0);
	printf("\n === b: ");
	mat_print(b0);

	plain_iteration(A0, b0, n, eps);
	jakobi(A0, b0, n, eps);
	zeidel(A0, b0, n, eps);

	mat_free(A0);
	mat_free(b0);
	/* ====== */

	return 0;
}
