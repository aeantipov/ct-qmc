#ifdef orb1
#define n_part_ 1
#elif defined orb2
#define n_part_ 2
#elif defined orb3
#define n_part_ 3
#elif defined orb4
#define n_part_ 4
#elif defined orb5
#define n_part_ 5
#elif defined orb6
#define n_part_ 6
#endif

#ifdef zone1
#define n_zone_ 1
#elif defined zone2
#define n_zone_ 2
#elif defined zone3
#define n_zone_ 3
#elif defined zone4
#define n_zone_ 4
#elif defined zone5
#define n_zone_ 5
#elif defined zone6
#define n_zone_ 6
#endif

typedef double n_type; //accuracy
int const n_zone=n_zone_; 
int const n_part=n_part_;

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
   Hw=Hw-Delta[z][wn];// static int f=0; //if (f<50) {cout<<w<<"  "<<Hw.x[0][0]-I*w<<"\n";f++;}
   Inverse(Hw);
   return Hw;
;}

#ifdef zone2

void W(point & r1, point_ & r1_, point & r2, point_ & r2_, n_type &u, n_type & a1, n_type &a2)
//"perturbation generator"
//generates time points and zone numbers, so that U=\bar{u}
{
double tau=beta*rnd();  
r1.t=tau; r1_.t=tau; r2.t=tau; r2_.t=tau;
  
int i=rnd(n_part);
int j=rnd(n_part);
while (j==i) j=rnd(n_part);
n_type C=n_part;

int n_terms=5, term=rnd(n_terms);
if (term==0)
{
   r1.i=i; r1_.i=i; r2.i=i; r2_.i=i;
   int z=rnd(2);
   r1.z=z; r1_.z=z; r2.z=1-z; r2_.z=1-z;   
   u=n_terms*beta*U*C;
;}
if (term==1) // nn
{
   r1.i=i; r1_.i=i; r2.i=j; r2_.i=j;
   int z=rnd(2);
   r1.z=z; r1_.z=z; r2.z=1-z; r2_.z=1-z;   
   u=n_terms*beta*(U-2*J)*C*(C-1);
;}
if (term==2) //nn
{
   r1.i=i; r1_.i=i; r2.i=j; r2_.i=j;
   int z=rnd(2);
   r1.z=z; r1_.z=z; r2.z=z; r2_.z=z;   
   u=n_terms*beta*(U-3*J)*C*(C-1);
;}
if (term==3) //pair hopping
{
   r1.i=i; r1_.i=j; r2.i=j; r2_.i=i;
   int z=rnd(2);
   r1.z=z; r1_.z=z; r2.z=1-z; r2_.z=1-z;   
   u=-n_terms*beta*J*C*(C-1);
;}
if (term==4) //pair hopping
{
   r1.i=j; r1_.i=i; r2.i=j; r2_.i=i;
   int z=rnd(2);
   r1.z=z; r1_.z=z; r2.z=1-z; r2_.z=1-z;   
   u=-n_terms*beta*J*C*(C-1);
;}

 define_alpha(r1,r1_,r2,r2_,u,a1,a2);//alpha's are defined by a standard recipe
   
;}

#elif defined zone6
#ifdef orb1
void W(point & r1, point_ & r1_, point & r2, point_ & r2_, n_type &u, n_type & a1, n_type &a2)
//"perturbation generator"
//generates time points and zone numbers, so that U=\bar{u}
{
double tau=beta*rnd();  
r1.t=tau; r1_.t=tau; r2.t=tau; r2_.t=tau;
r1.i=0; r1_.i=0; r2.i=0; r2_.i=0;

int i=rnd(n_zone/2);
int j=rnd(n_zone/2);
while (j==i) j=rnd(n_zone/2);

int n_terms=5;
int term=rnd(n_terms);
n_type C=n_zone/2;

switch (term)
	{
	case 0 : // nn
		{
//		   r1.i=i; r1_.i=i; r2.i=i; r2_.i=i;
		   int z=rnd(2);
		   r1.z=2*i+z; r1_.z=2*i+z; r2.z=2*i+1-z; r2_.z=2*i+1-z;   
		   u=n_terms*beta*U*C;
		   break;
		;}
		
	case 1 : // nn
		{
//		   r1.i=i; r1_.i=i; r2.i=j; r2_.i=j;
		   int z=rnd(2);
		   r1.z=2*i+z; r1_.z=2*i+z; r2.z=2*j+1-z; r2_.z=2*j+1-z;   
		   u=n_terms*beta*(U-2*J)*C*(C-1);
		   break;
		;}
	case 2 : //nn
		{
//		   r1.i=i; r1_.i=i; r2.i=j; r2_.i=j;
		   int z=rnd(2);
		   r1.z=2*i+z; r1_.z=2*i+z; r2.z=2*j+z; r2_.z=2*j+z;   
		   u=n_terms*beta*(U-3*J)*C*(C-1);
		   break;
		;}
	case 3 : //pair hopping 1
		{
//		   r1.i=i; r1_.i=j; r2.i=j; r2_.i=i;
		   int z=rnd(2);
		   r1.z=2*i+z; r1_.z=2*j+z; r2.z=2*j+1-z; r2_.z=2*i+1-z;   
		   u=-n_terms*beta*J*C*(C-1);
		   break;
		;}
	case 4 : //pair hopping 2
		{
//		   r1.i=j; r1_.i=i; r2.i=j; r2_.i=i;
		   int z=rnd(2);
		   r1.z=2*j+z; r1_.z=2*i+z; r2.z=2*j+1-z; r2_.z=2*i+1-z;   
		   u=-n_terms*beta*J*C*(C-1);
		   break;
		;}
	}
 define_alpha(r1,r1_,r2,r2_,u,a1,a2);//alpha's are defined by a standard recipe
   
;}
#endif
#endif

#ifdef use_mpi
#include "..//main_mpi.cpp"
#else


#include "..//main.cpp"
#endif

