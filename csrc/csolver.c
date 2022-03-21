#include <stdio.h>
#include <stdlib.h>
#include "csolver.h"


void shm(double t, double* y, double* dy) {
	dy[0] = y[1];
	dy[1] = -1*y[0];
}

 /*void eulerstep(func f, double t, double *y, int d, double h, double *yn) {
	(*f)(t,y,dy);
	return a;
}*/

struct Solution eulersolve(func f, double h, double t0, double tf, double* y0, int d) {
	struct Solution solution;
	double t=t0;
	int solsize = (int)(tf-t0) / h;
	solution.ts = calloc(solsize, sizeof(double));
	solution.yss = calloc(solsize*d, sizeof(double));
	double *ys  = calloc(d, sizeof(double));
	double *dy  = calloc(d, sizeof(double));

	int i,j;

	for (i=0; i<d; i++) {
		ys[i] = y0[i];
	}

	for (i=0; i<solsize; i++) {
		solution.ts[i] = t;

		f(t,ys,dy);
		for (j=0; j<d; j++) {
			*(solution.yss + i*d + j) = ys[j];
			ys[j] += h*dy[j];
		}

		t += h;
	}

	free(ys);
	free(dy);

	return solution;
}


struct Solution rk4solve(func f, double h, double t0, double tf, double* y0, int d) {
	struct Solution solution;
	double t=t0;
	int solsize = (int)(tf-t0) / h;
	solution.ts = calloc(solsize, sizeof(double));
	solution.yss = calloc(solsize*d, sizeof(double));
	double *ys  = calloc(d, sizeof(double));
	double *dy  = calloc(d, sizeof(double));
	double *tmp = calloc(d, sizeof(double));

	double *k1, *k2, *k3, *k4;
	k1 = calloc(d, sizeof(double));
	k2 = calloc(d, sizeof(double));
	k3 = calloc(d, sizeof(double));
	k4 = calloc(d, sizeof(double));

	int i,j;

	for (i=0; i<d; i++) {
		ys[i] = y0[i];
	}

	for (i=0; i<solsize; i++) {
		solution.ts[i] = t;

		f(t,ys,k1);

		for (j=0;j<d; j++) {
			tmp[j] = ys[j]+h/2*k1[j];
		}
		f(t+h/2, tmp, k2);

		for (j=0;j<d; j++) {
			tmp[j] = ys[j]+h/2*k2[j];
		}
		f(t+h/2, tmp, k3);

		for (j=0;j<d; j++) {
			tmp[j] = ys[j]+h*k3[j];
		}
		f(t+h, tmp, k4);

		for (j=0; j<d; j++) {

			*(solution.yss + i*d + j) = ys[j];
			ys[j] += h/6*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
		}

		t += h;
	}

	free(ys);
	free(dy);
	free(tmp);
	free(k1);
	free(k2);
	free(k3);
	free(k4);

	return solution;
}


/*int main() {
	double* yss;
	double ys[2];

	double h = 0.01;
	double t0 = 0;
	double tf = 50;

	ys[0] = 10;
	ys[1] = 0;
	struct Solution solution = eulersolve(&shm, h, t0, tf, ys, 2);

	int i,j,n;
	n = (int)(tf-t0) / h;
	for (i=0; i<n; i++) {
		printf("%f ", solution.ts[i]);
		for (j=0; j<2; j++) {
			printf("%f ", *(solution.yss + i*2 + j));
		}
		printf("\n");
	}

	free(solution.ts);
	free(solution.yss);
}*/