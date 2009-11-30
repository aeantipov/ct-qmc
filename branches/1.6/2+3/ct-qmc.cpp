int const n_zone=2, n_part=2;


typedef double n_type; //accuracy
#include "..//headers.cpp"
			//standard headers are included
			//some starting things are declared





//================ Action ========================

n_type U=value("U"), J=value("J"), beta=value("beta");




Matrix & g0w(int z, int wn)
//the unperturbed Green function in a Fourier space;
//z - zone, wn - number of Matsubara frequency
{
   n_type w=(2*wn+1)*Pi/beta;
   static Matrix Hw(0);
    Hw=I*w;
   //Hw=0.5*I*(w+sqrt(w*w+1)); 
   //Hw=Hw-Delta[z][wn]; static int f=0; //if (f<50) {cout<<w<<"  "<<Hw.x[0][0]-I*w<<"\n";f++;}
   Inverse(Hw);
 	return Hw;
;}



void W(point & r1, point_ & r1_, point & r2, point_ & r2_, n_type &u, n_type & a1, n_type &a2)
//"perturbation generator"
//generates time points and zone numbers, so that U=\bar{u}
{
double tau=beta*rnd();  
r1.t=tau; r1_.t=tau; r2.t=tau; r2_.t=tau;
  
  
int n_terms=5, term=rnd(n_terms);
if (term==0)
{
   int i=rnd(2);
   r1.i=i; r1_.i=i; r2.i=i; r2_.i=i;
   int z=rnd(2);
   r1.z=z; r1_.z=z; r2.z=1-z; r2_.z=1-z;   
   u=n_terms*beta*U*2.;
;}
if (term==1)
{
   int i=rnd(2);
   r1.i=i; r1_.i=i; r2.i=1-i; r2_.i=1-i;
   int z=rnd(2);
   r1.z=z; r1_.z=z; r2.z=1-z; r2_.z=1-z;   
   u=n_terms*beta*(U-2*J)*2.;
;}
if (term==2)
{
   int i=rnd(2);
   r1.i=i; r1_.i=i; r2.i=1-i; r2_.i=1-i;
   int z=rnd(2);
   r1.z=z; r1_.z=z; r2.z=z; r2_.z=z;   
   u=n_terms*beta*(U-3*J)*2.;
;}
if (term==3)
{
   int i=rnd(2);
   r1.i=i; r1_.i=1-i; r2.i=1-i; r2_.i=i;
   int z=rnd(2);
   r1.z=z; r1_.z=z; r2.z=1-z; r2_.z=1-z;   
   u=-n_terms*beta*J*2.;
;}
if (term==4) 
{
   int i=rnd(2);
   r1.i=1-i; r1_.i=i; r2.i=1-i; r2_.i=i;
   int z=rnd(2);
   r1.z=z; r1_.z=z; r2.z=1-z; r2_.z=1-z;   
   u=-n_terms*beta*J*2.;
;}

 define_alpha(r1,r1_,r2,r2_,u,a1,a2);//alpha's are defined by a standard recipe
   
;}


#ifdef use_mpi
#include "..//main_mpi.cpp"
#else


#include "..//main.cpp"
#endif

