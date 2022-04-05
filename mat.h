#define _POSIX_C_SOURCE 200809L

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
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

	exit(1);
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
	printf("mat (%u x %u)\n", M->m, M->n);
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
	if (A->m != B->m || A->n != B->n)
		die("mat_add: dims do not match (%u x %u) and (%u x %u)",
		    A->m, A->n, B->m, B->n);

	unsigned pos;

	for (pos = 0; pos < A->m * A->n; pos++)
		A->data[pos] += B->data[pos];
}

mat *
mat_add_o(const mat *A, const mat *B)
{
	if (A->m != B->m || A->n != B->n)
		die("mat_add: dims do not match (%u x %u) and (%u x %u)",
		    A->m, A->n, B->m, B->n);

	mat *O = mat_new(A->m, A->n);
	unsigned pos;

	for (pos = 0; pos < A->m * A->n; pos++)
		O->data[pos] = A->data[pos] + B->data[pos];

	return O;
}

double
mat_norm1(const mat *M)
{
	if (M->n != M->m)
		die("mat_norm1: (%u x %u) matrix not square", M->m, M->n);

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
