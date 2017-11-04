#include "misc.h"

#ifdef _WIN32
double log1p(double x) {
	return log(1+x);
}
#endif

int
NumExtrema(Vector in, int min, int max) {
	int num = 2;
	for (int s=min+1; s<=max-1; s++) {
		while (in[s] > in[s-1] && s<max) {
			s++;
			if (in[s] < in[s-1]) {
				num++;
			}
		}
		while (in[s] < in[s-1] && s<max) {
			s++;
			if (in[s] > in[s-1]) {
				num++;
			}
		}
	}
	return num;
}
Vector
PlaceExtrema(Vector in, int min, int max, int num) {
	Vector place(1,num);
	place[1] = min;
	place[num] = max;
	num = 2;
	for (int s=min+1; s<=max-1; s++) {
		while (in[s] > in[s-1] && s<max) {
			s++;
			if (in[s] < in[s-1]) {
				place[num] = s-1;
				num++;
			}
		}
		while (in[s] < in[s-1] && s<max) {
			s++;
			if (in[s] > in[s-1]) {
				place[num] = s-1;
				num++;
			}
		}
	}
	return place;
}
Vector
ValueExtrema(Vector in, int min, int max, int num) {
	Vector value(1,num);
	value[1] = in[min];
	value[num] = in[max];
	num = 2;
	for (int s=min+1; s<=max-1; s++) {
		while (in[s] > in[s-1] && s<max) {
			s++;
			if (in[s] < in[s-1]) {
				value[num] = in[s-1];
				num++;
			}
		}
		while (in[s] < in[s-1] && s<max) {
			s++;
			if (in[s] > in[s-1]) {
				value[num] = in[s-1];
				num++;
			}
		}
	}
	return value;
}
