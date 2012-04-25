//========== SORRY. THIS IS JUST A (black) MESS =======
//I put here some general "library" things
//without any order








//===============  small things =======================

n_type num_zero=1e-6;


int INT_RANDOM=time(NULL);

int int_rnd()
{
//for (int j=1; j<3; j++)
{
	 int k=INT_RANDOM/127773;
	 INT_RANDOM=16807*(INT_RANDOM-k*127773)-2836*k;
	 if (INT_RANDOM<0) INT_RANDOM+=2147483647;
;}
	 return INT_RANDOM;
;}

int INT_RANDOM_AUX_VARIABLE=int_rnd()+int_rnd()+int_rnd()+int_rnd();  //randomising...

n_type rnd ()
{
	return int_rnd()/2147483647.0
;}

int rnd (int k)
{
	int d=2147483647%k, d1=2147483647-d, r=int_rnd();
	if (r>=d1) return rnd(k);
	return r/(d1/k)
;}



n_type sqr (n_type x)
{
        return x*x;
;}

typedef n_type None_Arg_function ();
typedef n_type One_Arg_function (n_type);
typedef n_type One_Arg_function_i1 (n_type , int );
typedef n_type One_Arg_function_i2 (n_type , int , int);
typedef n_type One_Arg_function_i3 (n_type , int , int, int);
typedef n_type None_Arg_Rfunction ();
typedef void One_Agr_subroutine (n_type t);
typedef void None_Agr_subroutine ();
typedef n_type *p_n_type;
typedef n_type *p_n_type;




//==========================  spline ===========================
class Spline
{
	int N;
	int flag;
	void def_y2(n_type, n_type);
	n_type x0, x1;
	n_type x_n (int n) {return x0+(n*(x1-x0))/N;} //almost not necessary to be equidistant...
	public:
	n_type *y, *y2;
	Spline () {N=0; y=NULL; y2=NULL;}
	Spline (int NN) {N=NN; y=new n_type[N+1]; y2=new n_type[N+1];}
	void re_init (int NN) {delete y; delete y2; N=NN; y=new n_type[N+1]; y2=new n_type[N+1];}

void define (One_Arg_function  f, n_type x00, n_type x11, n_type yp0, n_type ypn)
{
	flag=1;
	x0=x00; x1=x11;
	for (int i=0; i<N+1; i++) y[i]=f(x_n(i));
	def_y2(yp0, ypn)
;}

void define (One_Arg_function_i1 f, n_type x00, n_type x11, n_type yp0, n_type ypn, int t)
{
	flag=1;
	x0=x00; x1=x11;
	for (int i=0; i<N+1; i++) y[i]=f(x_n(i), t);
	def_y2(yp0, ypn)
;}

void define (One_Arg_function_i2 f, n_type x00, n_type x11, n_type yp0, n_type ypn, int t1, int t2)
{
	flag=1;
	x0=x00; x1=x11;
	for (int i=0; i<N+1; i++) y[i]=f(x_n(i), t1, t2);
	def_y2(yp0, ypn)
;}

void define (One_Arg_function_i3 f, n_type x00, n_type x11, n_type yp0, n_type ypn, int t1, int t2, int t3)
{
	flag=1;
	x0=x00; x1=x11;
	for (int i=0; i<N+1; i++) y[i]=f(x_n(i), t1, t2, t3);
	def_y2(yp0, ypn)
;}


void define_nat (One_Arg_function f, n_type x00, n_type x11)
{
	flag=0;
	x0=x00; x1=x11;
	for (int i=0; i<N+1; i++) y[i]=f(x_n(i));
	def_y2(0, 0)
;}

void define_nat (One_Arg_function_i1 f, n_type x00, n_type x11, int t)
{
	flag=0;
	x0=x00; x1=x11;
	for (int i=0; i<N+1; i++) y[i]=f(x_n(i), t);
	def_y2(0, 0)
;}

void define_nat (One_Arg_function_i2 f, n_type x00, n_type x11, int t1, int t2)
{
	flag=0;
	x0=x00; x1=x11;
	for (int i=0; i<N+1; i++) y[i]=f(x_n(i), t1, t2);
	def_y2(0, 0)
;}

void define_nat (One_Arg_function_i3 f, n_type x00, n_type x11, int t1, int t2, int t3)
{
	flag=0;
	x0=x00; x1=x11;
	for (int i=0; i<N+1; i++) y[i]=f(x_n(i), t1, t2, t3);
	def_y2(0, 0)
;}


	n_type val(n_type xc);
   n_type val_no_check(n_type xc);
};

void Spline::def_y2(n_type yp0, n_type ypn)
{
	int i,k;
	n_type p, qn, sig, un, *u;
	u=new n_type [N+1];
	if (flag==0 || yp0>0.999e30) {y2[0]=0; u[0]=0;}
		else
		{y2[0]=-0.5;
//		cout<<x1<<"   "<<x0<<"   "<<x_n(1)<<"   "<<x_n(0);
		u[0]=(3/(x_n(1)-x_n(0)))*((y[1]-y[0])/(x_n(1)-x_n(0))-yp0);}

	for (i=1; i<N; i++)
	{
		sig=(x_n(i)-x_n(i+1))/(x_n(i+1)-x_n(i-1));
		p=sig*y2[i-1]+2;
		y2[i]=(sig-1)/p;
		u[i]=(y[i+1]-y[i])/(x_n(i+1)-x_n(i))-(y[i]-y[i-1])/(x_n(i)-x_n(i-1));
		u[i]=(6*u[i]/(x_n(i+1)-x_n(i-1))-sig*u[i-1])/p;
	;}

	if (flag==0 || ypn>0.999e30) {qn=0; un=0;}
		else {qn=0.5; un=(3/(x_n(N)-x_n(N-1)))*(ypn-(y[N]-y[N-1])/(x_n(N)-x_n(N-1)));}

	y2[N]=(un-qn*u[N-1])/(qn*y2[N-1]+1);
	for (k=N-1; k>=0; k--)  y2[k]=y2[k]*y2[k+1]+u[k];
;}





n_type Spline::val(n_type x)
{
	n_type ksi=(x-x0)*N/(x1-x0), h=(x1-x0)/N;                    //the only place where the equdistance is used
	if (ksi<0 || ksi>=N)
	{
		n_type eps=1e-10;
		if (ksi<0) ksi+=eps;
		if (ksi>=N) ksi-=eps;
		if (ksi<0 || ksi>=N) cout<<"spline "<<ksi<<"\n"; return 0;
	;}
	int k=int (ksi);
	n_type a=(x_n(k+1)-x)/h, b=(x-x_n(k))/h;
	return a*y[k]+b*y[k+1]+(a*(a*a-1)*y2[k]+b*(b*b-1)*y2[k+1])*h*h/6;

;}

n_type Spline::val_no_check(n_type x)
{
	n_type ksi=(x-x0)*N/(x1-x0), h=(x1-x0)/N;                    //the only place where the equdistance is used #2

	int k=int (ksi);
	n_type a=(x_n(k+1)-x)/h, b=(x-x_n(k))/h;
	return a*y[k]+b*y[k+1]+(a*(a*a-1)*y2[k]+b*(b*b-1)*y2[k+1])*h*h/6;

;}


