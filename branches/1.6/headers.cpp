
#include "input.cpp"


int const N_max=5000;   //maximal perturbation order

int main_argc;  char * main_argv[256];

			//classes for operators
struct point {n_type t; int i; int z;   //creation operator position
point (){i=0; t=0; z=0;}
point (int z1, int i1, n_type t1) {z=z1;i=i1;t=t1;}
void operator = (point  & r) {t=r.t;i=r.i;z=r.z;}
;};

struct point_{n_type t; int i; int z;  //annihilation
point_ (){i=0; t=0; z=0;}
point_ (int z1, int i1, n_type t1) {z=z1;i=i1;t=t1;}
void operator = (point_ & r) {t=r.t;i=r.i;z=r.z;}
;};

ofstream au("G0.dat");



int WL_count=1;


double to0(n_type x)
{
	static n_type f=value ("numerical_zero_for_output");
	if (fabs(x)<f) return 0;
	return double (x);	
;}







int n_tau=int_value("N_tau"), WN_max=int_value("number_of_Matsubara_frequencies")+2;
      	//n_tau - auxiliary number; number of points in t-direction
      	//for spline creation


//small things
n_type Pi=4.*atan(n_type(1.0));
complex I(0,1);

void read_Delta();
void Import_Interaction (point &, point_ &, point &, point_ &, n_type &, n_type &, n_type &);
int wn_max=int_value("number_of_Matsubara_frequencies");


int VectorDimension=n_part;
#include "math_lib.cpp"
   	   //matrix manipulations, splines, random generator etc...
         //this file must be included exactly at this place, because it uses
         //the value of n_part

Matrix * Delta[n_zone];
         //this matrix contains values from the file "Delta.dat"
         //it can be used in g0w(z, wn)


void define_alpha(point & r1, point_ & r1_, point & r2, point_ & r2_, n_type &u, n_type & a1, n_type &a2)
{
	static n_type alpha=value("alpha");
	static n_type alpha_off=value("alpha_offdiagonal");
	static n_type n_aver=value("expected_occupancy");

	a1=alpha_off*(rnd()-0.5);a2=alpha_off*(rnd()-0.5);

	if (r1.i==r1_.i && r2.i==r2_.i && u>0)
	{
		n_type x,y;
		if (n_aver<0.5) {x=alpha+2*n_aver-1; y=1-alpha;} else {x=alpha; y=2*n_aver-alpha;} 		

		if (rnd()<0.5) {a1+=x; a2+=y;} else {a1+=y; a2+=x;} 
   			
	;}

	if (r1.i==r1_.i && r2.i==r2_.i && u<0)

   			{a1+=n_aver; a2+=n_aver;}
;}


