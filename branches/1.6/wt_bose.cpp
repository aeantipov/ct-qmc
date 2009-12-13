typedef double n_type; //accuracy
#include "headers.cpp"


ifstream Sw_str("Sw.dat");
ofstream St_str("St.dat");

double beta=value("beta");

int n_harm;

double w_max;

double * Sw;

double Sw_expect(int k)
{
	if (k<n_harm) return Sw[k];
   else {double r=(n_harm-1.)/k; return Sw[n_harm-1]*r*r;}
;}

double St(double t)
{
	double s=Sw_expect(0);
   int k=0;
   for (double w=2*Pi/beta; w<1000; w+=2*Pi/beta)
   {
   	s+=2*Sw_expect(k)*cos(w*t);
      k++;
   ;}
   return s/beta;
;}

void main(int argc, char **argv)
{
	n_harm=atoi(argv[1]);
	Sw=new double [n_harm];
	w_max=(n_harm-0.5)*2*Pi/beta;
   int k=0;
   for (double w=0; w<w_max; w+=2*Pi/beta)
   {
   double r; Sw_str>>r;
   	Sw_str>>Sw[k];
      k++;
   ;}
   for (double t=0; t<beta+1e-5; t+=beta/128)
   	St_str<<t<<"    "<<St(t)<<"\n";
;}
