#include "aritmetica.h"
#include "aritmetica.h"

__global__
void FFTKernel(double* data, unsigned long istep, unsigned long m, unsigned long mmax, double wr, double wi, unsigned long tam){
	int i = (threadIdx.x + blockDim.x * blockIdx.x)*(istep)+(m);
	if (i <= (tam))
	{
		unsigned long j = i + (mmax);
		double tempr = (wr)* data[j - 1] - (wi)* data[j];
		double tempi = (wr)* data[j] + (wi)* data[j - 1];

		data[j - 1] = (data[i - 1] - tempr);
		data[j] = (data[i] - tempi);
		data[i - 1] = (data[i - 1] + tempr);
		data[i] = (data[i] + tempi);
	}
}
void fourpara(double* data, unsigned long* istep, unsigned long* m, unsigned long* mmax, double* wr, double* wi, unsigned long* tam)
{
	int size = (*tam) * sizeof(double);
	double *d_data;
	cudaMalloc((void **)&d_data, size);
	cudaMemcpy(d_data, data, size, cudaMemcpyHostToDevice);
	FFTKernel << < ceil(*tam/1024), 1024 >> > (d_data, *istep, *m, *mmax, *wr, *wi, *tam);
	cudaMemcpy(data, d_data, size, cudaMemcpyDeviceToHost);
	cudaFree(d_data);
}
void aritmetica::four(double* data, unsigned long nn)
{
	unsigned long n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta;

	n = nn << 1;
	j = 1;
	for (i = 1; i<n; i += 2) {
		if (j>i) {
			swap(data[j - 1], data[i - 1]);
			swap(data[j], data[i]);
		}
		m = nn;
		while (m >= 2 && j>m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	};
	mmax = 2;
	while (n>mmax) {
		istep = mmax << 1;
		theta = -(2 * M_PI / mmax);
		wtemp = sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m < mmax; m += 2) {
			fourpara(data, &istep, &m, &mmax, &wr, &wi, &n);
			wtemp = wr;
			wr += wr*wpr - wi*wpi;
			wi += wi*wpr + wtemp*wpi;
		}
		mmax = istep;
	}
}

__global__
void IFFTKernel(double* data, unsigned long istep, unsigned long m, unsigned long mmax, double wr, double wi, unsigned long tam){
	int i = (threadIdx.x + blockDim.x * blockIdx.x)*(istep)+(m);
	if (i <= (tam))
	{
		unsigned long j = i + (mmax);
		double tempr = (wr)* data[j - 1] - (wi)* data[j];
		double tempi = (wr)* data[j] + (wi)* data[j - 1];

		data[j - 1] = (data[i - 1] - tempr) / 2;
		data[j] = (data[i] - tempi) / 2;
		data[i - 1] = (data[i - 1] + tempr) / 2;
		data[i] = (data[i] + tempi) / 2;
	}
}
void fourIpara(double* data, unsigned long* istep, unsigned long* m, unsigned long* mmax, double* wr, double* wi, unsigned long* tam)
{
	int size = (*tam) * sizeof(double);
	double *d_data;
	cudaMalloc((void **)&d_data, size);
	cudaMemcpy(d_data, data, size, cudaMemcpyHostToDevice);
	IFFTKernel << <ceil(*tam / 1024), 1024 >> > (d_data, *istep, *m, *mmax, *wr, *wi, *tam);
	cudaMemcpy(data, d_data, size, cudaMemcpyDeviceToHost);
	cudaFree(d_data);
}
void aritmetica::fourI(double* data, unsigned long nn)
{
	unsigned long n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta;

	n = nn << 1;
	j = 1;
	for (i = 1; i<n; i += 2) {
		if (j>i) {
			swap(data[j - 1], data[i - 1]);
			swap(data[j], data[i]);
		}
		m = nn;
		while (m >= 2 && j>m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	};
	mmax = 2;
	while (n>mmax) {
		istep = mmax << 1;
		theta = -(2 * M_PI / mmax);
		wtemp = sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m < mmax; m += 2) {
			fourIpara(data, &istep, &m, &mmax, &wr, &wi, &n);
			wtemp = wr;
			wr += wr*wpr - wi*wpi;
			wi += wi*wpr + wtemp*wpi;
		}
		mmax = istep;
	}
}
void aritmetica::MultComple(double *X, double *Xi, double *Y, double *Yi, double* Resp, double* Respi)
{
	if (*Xi == 0 && *Yi == 0)
		*Resp = (*X) * (*Y);
	else
	{
		*Resp = ((*X) * (*Y)) - ((*Xi) * (*Yi));
		*Respi = ((*X) * (*Yi)) + ((*Xi) * (*Y));
	}
}
void aritmetica::ConjuComple(double* Resp, int tam)
{
	for (int i = 1; i < tam; i += 2)
	{
		if (Resp[i] != 0)
			Resp[i] = -Resp[i];
	}
}
double* aritmetica::Mult(double* X, double* Y, unsigned long tam)
{
	double* Resp = new double[tam];
	memset(Resp, 0, tam*sizeof(double));
	four(X, tam / 2);
	four(Y, tam / 2);
	for (unsigned long i = 0; i < tam; i += 2)
	{
		MultComple(X + i, X + i + 1, Y + i, Y + i + 1, Resp + i, Resp + i + 1);
	}
	ConjuComple(Resp, tam);
	fourI(Resp, tam / 2);
	return Resp;
}
ZZ aritmetica::powM(ZZ a, ZZ m, ZZ modulo,int NumB)
{

	ZZ respuesta;
	respuesta = 1;
	ZZ x;
	x = a;
	while (m != 0)
	{

		if ((m & 1) == 1)
		{
			respuesta = (MultiFourier(respuesta,x,NumB)) % modulo;
			//cout<<"respuesta_ "<<respuesta<<" x: "<<x<<endl;


		}
		x = (MultiFourier(x, x, NumB)) % modulo;
		m >>= 1;
		// cout<<"m: "<<m<<" x: "<<x<<" respuesta: "<<respuesta<<endl;
		//if(mod(x,modulo)==1) break;

	}
	//cout<<endl;
	return respuesta;
}
ZZ aritmetica::Blum(long n)
{
    ZZ N,semilla, p, q, bits,x, res, temp;

    p = to_ZZ("7171153257");
    q =to_ZZ("5");

    N = p * q;
    clock_t t;
    t=clock();
    semilla =t;
    x = semilla%N;
    res = 0, bits = 0;

    #pragma omp parallel for
    for(int i=n; i>0; i--)
    {
        x = powM(x, to_ZZ(2), N, n);
        bits = x-((x>>1)<<1);
        power(temp,to_ZZ(2),(i-1));
        res += bits*temp;
    }
    return res;
}
ZZ aritmetica::aleatorioBits(long long i)
{
    ZZ d =Blum(i);
    // cout<<d<<endl;
    if((d &1)==0)
        return d+1;
    else
        return d;
}
/*ZZ aritmetica::generaPrimo(long long bits)
{
    ZZ n =aleatorioBits(bits);
    while(MillerWitness((n),to_ZZ(80))==0)
        n=aleatorioBits(bits);
    return n;
}*/
ZZ aritmetica::generaPrimo(long long bits)
{
    ZZ n;
   
   do {
        n=aleatorioBits(bits);
        //cout<<n<<endl<<endl;
        //cout<<ProbPrime(n)<<endl<<endl;
    } while(ProbPrime(n)==0);
    return n;
}
ZZ aritmetica::MultiFourier(ZZ X, ZZ Y, int NumbitsRSA)
{
	stringstream convertX;
	convertX << X;
	string SX = convertX.str();
	stringstream convertY;
	convertY << Y;
	string SY = convertY.str();
	int xt = SX.size() - 1;
	int yt = SY.size() - 1;
	double * XVec = new double[NumbitsRSA * 2];
	double * YVec = new double[NumbitsRSA * 2];
	memset(XVec, 0, NumbitsRSA * 2 * sizeof(double));
	memset(YVec, 0, NumbitsRSA * 2 * sizeof(double));
	for (int i = xt, j = 0; i >= 0; i--, j += 2)
	{
		char nume = SX.at(i);
		XVec[j] = atof(&nume);
	}

	for (int i = yt, j = 0; i >= 0; i--, j += 2)
	{
		char nume = SY.at(i);
		YVec[j] = atof(&nume);
	}
	SX.clear();
	SY.clear();
	double* Mu = Mult(XVec, YVec, NumbitsRSA * 2);
	delete(XVec);
	delete(YVec);
	ZZ Respuesta = to_ZZ(0);
	ZZ diez = to_ZZ(1);
	for (int i = 0; Mu[i] >= 1; i += 2)
	{
		string X = to_string(Mu[i]);
		Respuesta += to_ZZ(X.c_str()) * diez;
		diez *= to_ZZ(10);
	}
	return Respuesta;
}




