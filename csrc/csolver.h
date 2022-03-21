#ifndef CSOLVER_H
#define CSOLVER_H

typedef void (*func)(double, double*, double*);

struct Solution {
	double* ts;
	double* yss;
};

struct Solution eulersolve(func f, double h, double t0, double tf, double* y0, int d);
struct Solution rk4solve(func f, double h, double t0, double tf, double* y0, int d);



#endif