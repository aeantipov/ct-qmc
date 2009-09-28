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
//int const n_zone=2; int n_part=1;
#include "..//headers.cpp"
			//standard headers are included
			//some starting things are declared

//================ Action ========================

n_type U=value("U"), beta=value("beta"), mu=value("mu"), delta=value("Delta");
n_type alpha=value("alpha");



Matrix & g0w(int z, int wn)
{
	n_type w=(2*wn+1)*Pi/beta;
   static Matrix Hw(0), d(0);
   if (z==0 || z==1) d.x[0][0]=-delta/2.; else d.x[0][0]=delta/2.;
   Hw=mu+I*w; Hw=Hw-d; Hw=Hw-Delta[z][wn];
   Inverse(Hw);
 	return Hw;
;}



void W(point & r1, point_ & r1_, point & r2, point_ & r2_, n_type &u, n_type & a1, n_type &a2)
{
	int i1=rnd(n_zone), i2;
   do i2=rnd(n_zone);while (i1==i2);

   r1.z=i1;r1_.z=i1; r2.z=i2;r2_.z=i2;

   r1.i=0; r1_.i=0; r2.i=0; r2_.i=0;
   n_type t=beta*rnd();

   r1.t=t; r1_.t=t; r2.t=t; r2_.t=t;

   u=n_zone*(n_zone-1)*beta*U/2;

//	if (rnd()<0.5) {a1=alpha-0.5; a2=1-alpha;} else {a1=1-alpha; a2=alpha-0.5;}
   define_alpha(r1,r1_,r2,r2_,u,a1,a2);//alpha's are defined by a standard recipe
;}


#ifdef mpi
#include "..//main_mpi.cpp"
#else
#include "..//main.cpp"
#endif
