
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <iostream>
#include<math.h>
#include <fstream>
#include <Windows.h>
#include "RSA.h"
double performancecounter_diff(LARGE_INTEGER *a, LARGE_INTEGER *b)
{
	LARGE_INTEGER freq;
	QueryPerformanceFrequency(&freq);
	return (double)(a->QuadPart - b->QuadPart) / (double)freq.QuadPart;
}
int main()
{

	RSA a(1024);
	string s = a.cifrarMensaje("Mi mama Me Mi r 0876543 LOL");
	cout << "Mensaje Encriptado: " << s;
	cout << a.descifrarMensaje(s) << endl;
	system("Pause");
	return 0;
}
