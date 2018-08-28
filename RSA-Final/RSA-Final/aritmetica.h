#pragma once
#include<math.h>
#include <algorithm> 
#include <iostream>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include<time.h>
#include <NTL/ZZ.h>
#include <sstream> 
#define M_PI 3.14159265358979323846
using namespace NTL;
using namespace std;

class aritmetica
{
public:
	/**
	Funcion de Transformada de Fourier en CUDA (Recursiva)
	data: Sera un vector que represente el polinomio del numero
	n : Sera una raiz primitiva por defecto usaremos -1 ya que es valida como raiz primitiva
	tam: tamaño del vector data, por razones de formula el vector debera ser de un tamaño que cumpla que sea un 2^n.
	*/
	ZZ Blum(long n);
	ZZ aleatorioBits(long long i);
	ZZ generaPrimo(long long bits);
	ZZ powM(ZZ a, ZZ m, ZZ modulo,int NumB);
	void four(double* data, unsigned long nn);
	void fourI(double* data, unsigned long nn);
	double* Mult(double* X, double* Y, unsigned long tam);
	void MultComple(double *X, double *Xi, double *Y, double *Yi, double* Resp, double* Respi);
	void ConjuComple(double* Resp, int tam);
	ZZ MultiFourier(ZZ X, ZZ Y, int NumbitsRSA);
};