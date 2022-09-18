#include <stdio.h> 
#include <math.h>

double a1 =  0.254829592;
double a2 = -0.284496736;
double a3 =  1.421413741;
double a4 = -1.453152027;
double a5 =  1.061405429;
double p  =  0.3275911;

double SQRT_2 = 1.41421356237;

// Error Function
double erf(double x) {
	int sign = 1;
	if (x < 0) {
		sign = -1;
	}
	x = fabs(x);

	double t = 1.0 / (1.0 + p*x);

	double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

	return sign * y;
}

// Cumulative Distribution Function
double cdf(double x) {
	double y = (erf(x / SQRT_2) + 1.0) / 2.0;
	return y;
}

// D1 Black Scholes
double D1(double S, double K, double T, double r, double sigma) {
	double d1 = (log(S/K) + (r + pow(sigma, 2)/2)) / (sigma * sqrt(T));
	return d1;
}

// D2 Black Scholes
double D2(double d1, double sigma, double T) {
	double d2 = d1 - sigma * sqrt(T);
	return d2;

}

// Black Scholes Call
double BS_CALL(double S, double K, double T, double r, double sigma) {
	double d1 = D1(S, K, T, r, sigma);
	double d2 = D2(d1, sigma, T);
	double C = S * cdf(d1) - K * exp(-r*T) * cdf(d2);

	return C;
}

// Black Scholes Put
double BS_PUT(double S, double K, double T, double r, double sigma) {
	double d1 = D1(S, K, T, r, sigma);
	double d2 = D2(d1, sigma, T);
	double C = K * exp(-r*T) * cdf(-d2) - S * cdf(-d1);
	
	return C;
}

// Black Scholes Call Delta
double delta_CALL(double S, double K, double T, double r, double sigma) {
	double d1 = D1(S, K, T, r, sigma);
	double delta = cdf(d1);
	return delta;
}

// Black Scholes Put Delta
double delta_PUT(double S, double K, double T, double r, double sigma) {
	double d1 = D1(S, K, T, r, sigma);
	double delta = cdf(d1) - 1;
	return delta;
}

double main(void) {

	double params[4];

	for (int i=0; i<=4; i++) {
		scanf("%lf", &params[i]);
	}

	double S = params[0];
	double K = params[1];
	double T = params[2];
	double r = params[3];
	double sigma = params[4];

	double call = BS_CALL(S, K, T, r, sigma);
	double put = BS_PUT(S, K,T, r, sigma);

	printf("CALL: %lf\n", call);
	printf("PUT:  %lf\n", put);	

	return 0;
}
