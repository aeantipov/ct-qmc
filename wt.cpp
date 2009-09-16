int const n_zone=2, n_part=1;


typedef double n_type; //accuracy
#include "headers.cpp"

n_type beta=value("beta");
n_type mu1_Gt, mu2_Gt;    complex * GM00_st;



n_type Gw0i_Gt(n_type w)
{
	return  -0.5*w/(w*w+mu1_Gt*mu1_Gt)-0.5*w/(w*w+mu2_Gt*mu2_Gt);
;}

n_type Gw0r_Gt(n_type w)
{
	return -0.5*mu1_Gt/(w*w+mu1_Gt*mu1_Gt)-0.5*mu2_Gt/(w*w+mu2_Gt*mu2_Gt);
;}

n_type Fermi_dist_Gt(n_type mu, n_type t)
{
	return exp(mu*(beta-t))/(1+exp(mu*beta));
;}

n_type Gtinf_Gt(n_type t)
{
   return 0.5*Fermi_dist_Gt(mu1_Gt,t)+0.5*Fermi_dist_Gt(mu2_Gt,t);
;}

n_type Gt_single(n_type t)
{
	n_type s=0;
   int wn=0;
   for (n_type w=Pi/beta; w<2*Pi*wn_max/beta; w+=2*Pi/beta)
   {
   	s+=2*( (imag(GM00_st[wn])-Gw0i_Gt(w))*sin(w*t)+
             (real(GM00_st[wn])-Gw0r_Gt(w))*cos(w*t))/beta;
      wn++;
   ;}
   return s-Gtinf_Gt(t);
;}

void Gt_write()
{
   ofstream Gt_str("wt.dat");
   n_type w=(2*wn_max-1)*Pi/beta;

   {
   	n_type re=real(GM00_st[wn_max-1]), im=imag(GM00_st[wn_max-1]);
      n_type sq=re*re+im*im, eps=1e-10-w*(im+w*sq)*(1+4*w*(im+w*sq));
      if (eps>0)
      {
      	mu1_Gt=(-w*re+sqrt(eps))/(im+2*w*sq);
      	mu2_Gt=(-w*re-sqrt(eps))/(im+2*w*sq);
      ;}
      else
      {
      	mu1_Gt=0; mu2_Gt=0;
      ;}
   ;}


   for (n_type t=0; t<beta+1e-5; t+=beta/256)
   {
   	Gt_str<<t;
        Gt_str<<"    "<<Gt_single(t);
      Gt_str<<"   "<<0<<"\n";
   ;}
;}



int main (int argc, char *argv[])
{
	GM00_st=new complex [wn_max+1];
   ifstream inp(argv[1]);
   for (int wn=0; wn<wn_max; wn++)
   	{n_type x,y=0;  inp>>x;  inp>>y; GM00_st[wn]=complex(x,y);}
   Gt_write();
	return 1;
;}
