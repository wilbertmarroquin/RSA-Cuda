#include "RSA.h"


RSA::RSA(int bits)
{
    datos=" abcdefghijklmnopqrstuvwxyz ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
	NumbitsRSA = bits;
    do
    {
        p=GenPrime_ZZ(bits);
        //p=Funciones.generaPrimo(bits);
        q=GenPrime_ZZ(bits);
        //q=Funciones.generaPrimo(bits);
    }
    while(p==q);

    n=p*q;
    on=(p-1)*(q-1);
    //e=Funciones.generaPrimo(bits);
    e=GenPrime_ZZ(bits);
    do {
        e=Funciones.aleatorioBits(bits);
    }
    while((GCD(e,on))!=1);


    d=InvMod(e,on);
}
ZZ  RSA::Encriptar(string msj)
{
    ZZ j;
    j = datos.find(msj);
    j = PowerMod(j, e, n);
    return j;
}
char RSA::Desencriptar(ZZ j)
{
    ZZ r = PowerMod(j, d, n);
    int s = to_int(r);
    cout << s<<endl;
    return datos[s];
}
template <class T>
string RSA::convert(T num,bool val)
{
    stringstream convert;
    convert<<num;
    return convert.str();
}
ZZ RSA::convertS(string Text)
{
    stringstream convert(Text);
    ZZ num;
    if ( !(convert >> num) )
        return  to_ZZ(-1);
    return num;
}
string RSA::ceros(string t,ZZ number)
{
    while((to_ZZ(t.size()))<lenNumber(number))
    {
        t="0"+t;
    }
    return t;
}
string RSA::cerosFin(string t,ZZ number)
{
    do
    {
        t+="0";
    }
    while(to_ZZ(t.length())%number!=0);
    return t;
}
string RSA::ceros(string t)
{
    return ceros(t,n);
}
ZZ RSA::lenNumber(ZZ dat)
{
    ZZ d;
    d=0;
    do
    {
        d++;
        dat=dat/10;
    }
    while(dat!=0);
    return d;
}
string  RSA::tratarMensaje(string msj)
{
    string total;
    for(int i=0; i<msj.size(); i++)
    {
        int dat= datos.find(msj[i]);
        //cout<<"mi dat "<<dat<<endl;
        string aux;
        aux=convert<int>(dat);

        while((aux.length())<lenNumber(to_ZZ(datos.length()-1)))
        {
            aux="0"+aux;
            cout<<aux<<endl;
        }
        //cout<<"aux: "<<aux<<endl;
        total+=aux;
    }
    return total;
}
void RSA::proc_paralell(long &i,ZZ & val,string&d,string&total)
{
        string t=d.substr(i,to_long(val));
        ZZ temp(INIT_VAL, t.c_str());
        cout<<"t vale "<<temp<<endl;

        t=convert(PowerMod(temp,e,n));
        cout<<"mi t vale "<<t<<endl;
        if(t.size()<lenNumber(n))
        {
            total+=ceros((t));
            cout<<ceros(t)<<endl;
        }
        else
        {
            total+=t;
            cout<<t<<endl;
        }
}
string RSA::cifrarMensaje(string datos)
{
    string d=tratarMensaje(datos);

    ZZ val=lenNumber((n))-1;
    string total;
    string aux;


    if((to_ZZ(d.length())%val)!=0)
        d=cerosFin(d,val);

    for(long i=0,j=1; i<to_long(d.length()); i+=to_long(val))
    {
        if(j<thread::hardware_concurrency())
        {
            async(launch::async,&RSA::proc_paralell,this,ref(i),ref(val),ref(d),ref(total)).get();j++;
        }
        else
        {
            proc_paralell(i,val,d,total);j--;
        }
    }
    return total;
}
string Postceros(string t, ZZ number)
{
	while ((to_ZZ(t.size()))<(number))
	{
		t = "0" + t;
		//    cout<<"CENTRO DE CDEROS "<<t<<endl;
	}
	//cout<<"t resultante "<<t<<endl;
	return t;
}
void RSA::dec_paralell(long &i, ZZ & val, string&total, string&datos)
{
	string t = datos.substr(i, to_long(val));
	ZZ temp(INIT_VAL, t.c_str());
	ZZ aux = PowerMod(temp, d, n);
	if (lenNumber(aux)<val)
	{
		t = Postceros(convert(aux), val - 1);
	}
	total += t;
}
string RSA::descifrarMensaje(string datos)
{
	ZZ val = lenNumber(n);
	string total, t;
	ZZ aux;

	for (long i = 0, j = 1; i<to_long(datos.length()); i += to_long(val))
	{
		if (j<std::thread::hardware_concurrency())
		{
			async(launch::async, &RSA::dec_paralell, this, ref(i), ref(val), ref(total), ref(datos)).get(); j++;
		}
		else
		{
			dec_paralell(i, val, total, datos); j--;
		}
	}
	val--;
	for (long i = 0; i<to_long(total.length()); i += to_long(val))
	{
		string t = datos.substr(i, to_long(val));

		if (t.size()<val)
		{
			t = Postceros(t, val - 1);
		}


	}

	string res;
	for (int i = 0; i<total.size(); i += to_int(lenNumber(to_ZZ(this->datos.size() - 1))))
	{
		int dat = to_int(convertS(total.substr(i, to_int(lenNumber(to_ZZ(this->datos.size() - 1))))));
		res += this->datos[dat];
	}
	return res;
}
RSA::~RSA()
{
}
