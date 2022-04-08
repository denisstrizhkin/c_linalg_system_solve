#define _POSIX_C_SOURCE 200809L

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

typedef struct mat {
	unsigned m;
	unsigned n;
	double *data;
} mat;

void
die(char *fmt, ...)
{
	va_list ap;

	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
	va_end(ap);
	fputc('\n', stderr);

	exit(0);
}

void
die_dim_match(char *func_name, const mat *A, const mat *B)
{
	if (A->m != B->m || A->n != B->n)
		die("%s: dims do not match (%u x %u) and (%u x %u)",
		    func_name, A->m, A->n, B->m, B->n);
}

mat *
mat_new(unsigned m, unsigned n)
{
	mat *M = malloc(sizeof *M);

	M->m = m;
	M->n = n;
	M->data = malloc((sizeof *M->data) * m * n);

	return M;
}

void
mat_free(mat *M)
{
	free(M->data);
	free(M);
}

mat *
mat_copy(const mat *M)
{
	mat *O = mat_new(M->m, M->n);

	memcpy(O->data, M->data, (sizeof *M->data) * M->m * M->n);

	return O;
}

double
mat_get(const mat *M, unsigned m, unsigned n)
{
	unsigned pos = m * M->n + n;
	if (pos >= M->m * M->n)
		die("mat_get: wrong position m: %u, n: %u", m, n);

	return M->data[pos];
}

void
mat_set(mat *M, unsigned m, unsigned n, double val)
{
	unsigned pos = m * M->n + n;
	if (pos >= M->m * M->n)
		die("mat_set: wrong position m: %u, n: %u", m, n);
	
	M->data[pos] = val;
}

mat *
mat_zeros(unsigned m, unsigned n)
{
	mat *M = mat_new(m, n);
	unsigned i = 0;
	for (; i < m * n; i++)
		M->data[i] = 0.0;
	return M;
}

mat *
mat_I(unsigned n)
{
	mat *M = mat_zeros(n, n);
	unsigned i, j;
	for(i = 0, j = 0; i < n; i++, j++)
		mat_set(M, i, j, 1.0);
	return M;
}

void
mat_print(const mat *M)
{
	if (M->m > 10 || M->n > 10)
		printf("mat_print: matrix is to big to be printed\n");
	printf("mat (%u x %u) %u %u\n", M->m, M->n, M, M->data);
	unsigned i, j;
	for (i = 0; i < M->m; i++) {
		for (j = 0; j < M->n; j++) {
			if (j == 0)
				printf("%lf", mat_get(M, i, j));
			else
				printf("\t%lf", mat_get(M, i, j));
		}
	printf("\n");
	}
}

mat *
mat_ltrig_rnd(unsigned n, unsigned dens)
{
	if (dens > 100)
		die("mat_ltrig_rnd: density is number between 0 and 100");
	mat *M = mat_zeros(n, n);
	unsigned i, j;
	unsigned tmp;

	srand(time(0));
	for (i = 0; i < M->m; i++) {
		for (j = 0; j < i + 1; j++) {
			tmp = rand() % (101);
			if (tmp <= dens)
				mat_set(M, i, j, (double) rand() / RAND_MAX);
		}
	}

	return M;
}

mat *
mat_utrig_rnd(unsigned n, unsigned dens)
{
	if (dens > 100)
		die("mat_ltrig_rnd: density is number between 0 and 100");
	mat *M = mat_zeros(n, n);
	unsigned i, j;
	unsigned tmp;

	srand(time(0));
	for (i = 0; i < M->m; i++) {
		for (j = i; j < n; j++) {
			tmp = rand() % (101);
			if (tmp <= dens)
				mat_set(M, i, j, (double) rand() / RAND_MAX);
		}
	}

	return M;
}

mat *
mat_rnd(unsigned m, unsigned n, unsigned dens)
{
	if (dens > 100)
		die("mat_ltrig_rnd: density is number between 0 and 100");
	mat *M = mat_zeros(m, n);
	unsigned i, j;
	unsigned tmp;

	srand(time(0));
	for (i = 0; i < M->m; i++) {
		for (j = 0; j < M->n; j++) {
			tmp = rand() % (101);
			if (tmp <= dens)
				mat_set(M, i, j, (double) rand() / RAND_MAX);
		}
	}

	return M;
}

void
mat_scale(mat *M, double a)
{
	unsigned pos;

	for (pos = 0; pos < M->m * M->n; pos++)
		M->data[pos] *= a;
}

mat *
mat_scale_o(const mat *M, double a)
{
	mat *O = mat_new(M->m, M->n);
	unsigned pos;

	for (pos = 0; pos < M->m * M->n; pos++)
		O->data[pos] = M->data[pos] * a;

	return O;
}

void
mat_add(mat *A, const mat *B)
{
	die_dim_match("mat_add", A, B);

	double *pa = A->data;
	double *pb = B->data;
	double *pe = &A->data[A->m * A->n];

	for (; pa != pe; pa++, pb++)
		*pa += *pb;
}

mat *
mat_add_o(const mat *A, const mat *B)
{
	die_dim_match("mat_add_o", A, B);

	mat *O = mat_copy(A);
	mat_add(O, B);

	return O;
}

double
mat_norm1(const mat *M)
{
	unsigned i, j;
	double sum, norm = 0.0;

	for (j = 0; j < M->n; j++) {
		sum = 0.0;
		for (i = 0; i < M->n; i++) {
			sum += fabs(mat_get(M, i, j));
		}
		norm = fmax(norm, sum);
	}

	return norm;
}

void
mat_sub(mat *A, const mat *B)
{
	die_dim_match("mat_sub", A, B);

	double *pa = A->data;
	double *pb = B->data;
	double *pe = &A->data[A->m * A->n];

	for (; pa != pe; pa++, pb++)
		*pa -= *pb;
}

mat *
mat_sub_o(const mat *A, const mat *B)
{
	die_dim_match("mat_sub_o", A, B);

	mat *O = mat_copy(A);
	mat_sub(O, B);

	return O;
}

mat *
mat_dot(const mat *A, const mat *B)
{
	if(A->n != B->m)
		die("mat_dot: (%lu x %lu) & (%lu x %lu)", A->m, A->n, B->m, B->n);

	mat *O = mat_new(A->m, B->n);
	unsigned i, j, k;
	for (i = 0; i < O->m; i++) {
		for (j = 0; j < O->n; j++) {
			O->data[i * O->n + j] = 0.0;
			for (k = 0; k < A->n; k++) {
				O->data[i * O->n + j] += A->data[i * A->n + k] *
				                         B->data[k * B->n + j];
			}
		}
	}
	
	return O;
}
