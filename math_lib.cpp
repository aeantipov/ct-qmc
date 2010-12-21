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

n_type norm2(ComplexType & x)
{
	return x.real()*x.real()+x.imag()*x.imag();
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

//==================================================================================================

struct Vector
{
	ComplexType *x;//[VectorDimension];
	int to_del;
	Vector();
	Vector(ComplexType *);
	Vector (Vector &);
	~Vector();
	void  operator = (Vector &);
	void  operator += (Vector &);
	void  operator -= (Vector &);
	void  operator *= (ComplexType);
	void  operator /= (ComplexType);
};

Vector::Vector()
{
	x= new ComplexType[VectorDimension]; to_del=1;
;}

Vector::Vector(ComplexType * p)
{
	x=p; to_del=0;
;}

Vector::Vector(Vector & v)
{
	x= new ComplexType[VectorDimension]; to_del=1;
	for (int i=0; i<VectorDimension; i++) x[i]=v.x[i];
;}

Vector::~Vector()
{
	if (to_del==1) delete x;
;}


void Vector:: operator = (Vector & a)
{
	for (int i=0; i<VectorDimension; i++) x[i]=a.x[i];
;}

void Vector:: operator += (Vector & a)
{
	for (int i=0; i<VectorDimension; i++) x[i]+=a.x[i];
;}

void Vector:: operator -= (Vector & a)
{
	for (int i=0; i<VectorDimension; i++) x[i]-=a.x[i];
;}

void Vector:: operator *= (ComplexType a)
{
	for (int i=0; i<VectorDimension; i++) x[i]*=a;
;}

void Vector:: operator /= (ComplexType a)
{
	for (int i=0; i<VectorDimension; i++) x[i]/=a;
;}

Vector operator - (Vector & a)
{
	Vector z;
	for (int i=0; i<VectorDimension; i++) z.x[i]=-a.x[i];
	return z;
;}


Vector operator + (Vector & a, Vector & b)
{
	Vector z;
	for (int i=0; i<VectorDimension; i++) z.x[i]=a.x[i]+b.x[i];
	return z;
;}

Vector operator - (Vector & a, Vector & b)
{
	Vector z;
	for (int i=0; i<VectorDimension; i++) z.x[i]=a.x[i]-b.x[i];
	return z;
;}

Vector operator * (ComplexType a, Vector & b)
{
	Vector z;
	for (int i=0; i<VectorDimension; i++) z.x[i]=a*b.x[i];
	return z;
;}

Vector operator * (Vector & b, ComplexType a)
{
	Vector z;
	for (int i=0; i<VectorDimension; i++) z.x[i]=b.x[i]*a;
	return z;
;}

Vector operator / (Vector & b, ComplexType a)
{
	Vector z;
	for (int i=0; i<VectorDimension; i++) z.x[i]=b.x[i]/a;
	return z;
;}

ComplexType operator * (Vector & a, Vector & b)
{
	ComplexType z=0;
	for (int i=0; i<VectorDimension; i++) z+=a.x[i]*b.x[i];
	return z;
;}

ComplexType sqr (Vector &a)
{
	return a*a;
;}

Vector VectorNull()
{
	Vector v;
	for (int i=0; i<VectorDimension; i++) v.x[i]=0;
	return v;
;}

void print(Vector & a)
{
	for (int i=0; i<VectorDimension; i++) cout<<real(a.x[i])<<"   "<<imag(a.x[i])<<"  ";
	cout<<"\n";
;}

void print(ofstream &f, Vector & a)
{
	for (int i=0; i<VectorDimension; i++) f<<i<<"   "<<real(a.x[i])<<"   "<<imag(a.x[i])<<"\n";
;}




struct Matrix
{
	ComplexType **x;//[VectorDimension];                //x[i] указывает на строку i
	int flag; //f==0 => delete after the first usage
	Matrix ();
	Matrix (ComplexType );
   	Matrix (n_type x);
	Matrix (Matrix &);
	Matrix (Vector&, Vector &); //direct product
	Matrix (Vector &); //diagonal
	~Matrix();
	void new_memory(); //explicitly get new memory for x...

   	void  operator = (ComplexType );
	void  operator = (n_type );
	void  operator = (Matrix &);
	void  operator += (Matrix &);
	void  operator -= (Matrix &);
	void  operator *= (ComplexType);
	void  operator /= (ComplexType);

	Vector column(int );
	Vector row(int );
   	void trans();
};

Matrix:: Matrix ()
{
	x=new ComplexType *  [VectorDimension];
	for (int i=0; i<VectorDimension; i++) x[i]=new ComplexType[VectorDimension];
	flag=1;
;}

Matrix:: Matrix (ComplexType  r)
{
	x=new ComplexType *  [VectorDimension];
	for (int i=0; i<VectorDimension; i++) x[i]=new ComplexType[VectorDimension];
	for (int l=0; l<VectorDimension; l++)
	for (int j=0; j<VectorDimension; j++) if (l==j) x[l][j]=r; else x[l][j]=0;
	flag=1;
;}

Matrix:: Matrix (n_type  r)
{
	x=new ComplexType *  [VectorDimension];
	for (int i=0; i<VectorDimension; i++) x[i]=new ComplexType[VectorDimension];
	for (int l=0; l<VectorDimension; l++)
	for (int j=0; j<VectorDimension; j++) if (l==j) x[l][j]=ComplexType(r,0); else x[l][j]=0;
	flag=1;
;}


Matrix::Matrix (Matrix & m)
{
	x=new ComplexType *  [VectorDimension];
	for (int i=0; i<VectorDimension; i++) x[i]=new ComplexType[VectorDimension];
	for (int l=0; l<VectorDimension; l++)
	for (int j=0; j<VectorDimension; j++) x[l][j]=m.x[l][j];
	flag=1;
	if (m.flag==0) m.~Matrix();
;}

Matrix::Matrix (Vector &v1, Vector & v2)
{
	x=new ComplexType *  [VectorDimension];
	for (int i=0; i<VectorDimension; i++) x[i]=new ComplexType[VectorDimension];
	for (int l=0; l<VectorDimension; l++)
	for (int j=0; j<VectorDimension; j++) x[l][j]=v1.x[l]*v2.x[j];
	flag=1;
;}

Matrix::Matrix (Vector &v)
{
	x=new ComplexType *  [VectorDimension];
	for (int i=0; i<VectorDimension; i++) x[i]=new ComplexType[VectorDimension];
	for (int l=0; l<VectorDimension; l++)
	for (int j=0; j<VectorDimension; j++)
   	if (l!=j) x[l][j]=0; else x[l][j]=v.x[l];
	flag=1;
;}




Matrix:: ~Matrix ()
{
	for (int i=0; i<VectorDimension; i++) delete[] x[i];
   	delete[] x;
;}

	
void Matrix:: new_memory()
{
	static int f=0;
	if (f==0)
   {
   x=new ComplexType *  [VectorDimension];
   for (int i=0; i<VectorDimension; i++) x[i]=new ComplexType[VectorDimension];
   ;}
   f=1;
   for (int i=0; i<VectorDimension; i++)
   for (int j=0; j<VectorDimension; j++) x[i][j]=0;
	flag=1;
;}





void Matrix:: operator = (n_type r)
{
	for (int l=0; l<VectorDimension; l++)
	for (int j=0; j<VectorDimension; j++) if (l==j) x[l][j]=r; else x[l][j]=0;
;}

void Matrix:: operator = (ComplexType r)
{
	for (int l=0; l<VectorDimension; l++)
	for (int j=0; j<VectorDimension; j++) if (l==j) x[l][j]=r; else x[l][j]=0;
;}





void Matrix:: operator = (Matrix & a)
{
	for (int i=0; i<VectorDimension; i++)
	for (int j=0; j<VectorDimension; j++) x[i][j]=a.x[i][j];
	if (a.flag==0) a.~Matrix();
;}



void Matrix:: operator += (Matrix & a)
{
	for (int i=0; i<VectorDimension; i++)
	for (int j=0; j<VectorDimension; j++) x[i][j]+=a.x[i][j];
	if (a.flag==0) a.~Matrix();
;}

void Matrix:: operator -= (Matrix & a)
{
	for (int i=0; i<VectorDimension; i++)
	for (int j=0; j<VectorDimension; j++) x[i][j]-=a.x[i][j];
	if (a.flag==0) a.~Matrix();
;}

void Matrix:: operator *= (ComplexType  a)
{
	for (int i=0; i<VectorDimension; i++)
	for (int j=0; j<VectorDimension; j++) x[i][j]*=a;
;}

void Matrix:: operator /= (ComplexType a)
{
	for (int i=0; i<VectorDimension; i++)
	for (int j=0; j<VectorDimension; j++) x[i][j]/=a;
;}


Matrix & operator - (Matrix & a)
{
	Matrix * z=new Matrix(); (*z).flag=0;
	for (int i=0; i<VectorDimension; i++)
	for (int j=0; j<VectorDimension; j++) (*z).x[i][j]=-a.x[i][j];
	if (a.flag==0) a.~Matrix();
	return *z;
;}



Matrix & operator + (Matrix & a, Matrix & b)
{
	Matrix * z=new Matrix(); (*z).flag=0;
	for (int i=0; i<VectorDimension; i++)
	for (int j=0; j<VectorDimension; j++) (*z).x[i][j]=a.x[i][j]+b.x[i][j];
	if (a.flag==0) a.~Matrix();
	if (b.flag==0) b.~Matrix();
	return *z;
;}


Matrix & operator - (Matrix & a, Matrix & b)
{
	Matrix * z=new Matrix(); (*z).flag=0;
	for (int i=0; i<VectorDimension; i++)
	for (int j=0; j<VectorDimension; j++) (*z).x[i][j]=a.x[i][j]-b.x[i][j];
	if (a.flag==0) a.~Matrix();
	if (b.flag==0) b.~Matrix();
	return *z;
;}



Matrix & operator * (Matrix & a, ComplexType b)
{
	Matrix * z=new Matrix(); (*z).flag=0;
	for (int i=0; i<VectorDimension; i++)
	for (int j=0; j<VectorDimension; j++) (*z).x[i][j]=a.x[i][j]*b;
	if (a.flag==0) a.~Matrix();
	return * z;
;}

Matrix & operator * (ComplexType b, Matrix & a)
{
	Matrix * z=new Matrix(); (*z).flag=0;
	for (int i=0; i<VectorDimension; i++)
	for (int j=0; j<VectorDimension; j++) (*z).x[i][j]=a.x[i][j]*b;
	if (a.flag==0) a.~Matrix();
	return * z;
;}

Matrix & operator / (Matrix & a, ComplexType b)
{
	Matrix * z=new Matrix(); (*z).flag=0;
	for (int i=0; i<VectorDimension; i++)
	for (int j=0; j<VectorDimension; j++) (*z).x[i][j]=a.x[i][j]/b;
	if (a.flag==0) a.~Matrix();
	return * z;
;}

Vector & operator * (Matrix & a, Vector & b)
{
	Vector * z=new Vector;
	for (int i=0; i<VectorDimension; i++)
	{
		(*z).x[i]=0;
		for (int j=0; j<VectorDimension; j++) (*z).x[i]+=a.x[i][j]*b.x[j];
	;}
	if (a.flag==0) a.~Matrix();
	return * z;
;}


Vector & operator * (Vector & b, Matrix & a)
{
	Vector *z=new Vector;
	for (int i=0; i<VectorDimension; i++)
	{
		(*z).x[i]=0;
		for (int j=0; j<VectorDimension; j++) (*z).x[i]+=b.x[j]*a.x[j][i];
	;}
	if (a.flag==0) a.~Matrix();
	return *z ;
;}


Matrix &  operator * (Matrix & a, Matrix & b)
{
	Matrix * z=new Matrix(); (*z).flag=0;
	for (int i=0; i<VectorDimension; i++)
	for (int j=0; j<VectorDimension; j++)
	{
		ComplexType w=0;
		for (int k=0; k<VectorDimension; k++) w+=a.x[i][k]*b.x[k][j];
		(*z).x[i][j]=w;
	;}
	if (a.flag==0) a.~Matrix();
	if (b.flag==0) b.~Matrix();
	return *z;
;}

Matrix operator %  (Matrix & a, Matrix & b) //elemant-by-element multiplication
{
	Matrix * z=new Matrix(); (*z).flag=0;
	for (int i=0; i<VectorDimension; i++)
	for (int j=0; j<VectorDimension; j++)
	{
		(*z).x[i][j]=a.x[i][j]*b.x[i][j];
	;}
	if (a.flag==0) a.~Matrix();
	if (b.flag==0) b.~Matrix();
	return *z;
;}



Vector Matrix:: row(int i)
{
	Vector v; for (int j=0; j<VectorDimension; j++) v.x[j]=x[i][j];
	return v;
;}

Vector Matrix:: column(int i)
{
	Vector v; for (int j=0; j<VectorDimension; j++) v.x[j]=x[j][i];
	return v;
;}

void Matrix:: trans()
{
	for (int i=0; i<VectorDimension; i++)
   {
   x[i][i]=conj(x[i][i]);
   for (int j=i+1; j<VectorDimension; j++)
   	{ComplexType a=x[i][j];
      x[i][j]=conj(x[j][i]);
      x[j][i]=conj(a);}
   ;}
;}



void SetUnity (Matrix & a)
{
	for (int i=0; i<VectorDimension; i++)
	for (int j=0; j<VectorDimension; j++)
		if (i==j) a.x[i][j]=1; else a.x[i][j]=0;
;}

void SetNull (Matrix & a)
{
	for (int i=0; i<VectorDimension; i++)
	for (int j=0; j<VectorDimension; j++)
		a.x[i][j]=0;
;}

Matrix & Conj (Matrix & a)
{
	Matrix * z=new Matrix(); (*z).flag=0;
	for (int i=0; i<VectorDimension; i++)
	for (int j=0; j<VectorDimension; j++)
	(*z).x[i][j]=conj(a.x[j][i]);
	if (a.flag==0) a.~Matrix();
	return *z;
;}



void print (Matrix & a)
{
	for (int i=0; i<VectorDimension; i++)
	{
		for (int j=0; j<VectorDimension; j++) cout<<real(a.x[i][j])<<"   "<<imag(a.x[i][j])<<"  ";
		cout<<"\n";
	;}
;}

void print (Matrix & a, ofstream & st)
{
	for (int i=0; i<VectorDimension; i++)
	{
		for (int j=0; j<VectorDimension; j++) st<<real(a.x[i][j])<<"   "<<imag(a.x[i][j])<<"  ";
		st<<"\n";
	;}
;}

ComplexType Sp(Matrix & M)
{
	ComplexType s=0;
	for (int i=0; i<VectorDimension; i++) s+=M.x[i][i];
	if (M.flag==0) M.~Matrix();
   return s;
;}



n_type norm2(Matrix & M)
{
	n_type s=0;
	for (int i=0; i<VectorDimension; i++)
	for (int j=0; j<VectorDimension; j++)
		s+=real(M.x[i][j]*conj(M.x[i][j]));
	if (M.flag==0) M.~Matrix();
	return s;
;}



void Inverse(Matrix & a, int improve=0) //Gauss with partial pivoting
//if flag!=0 => tries to improve accuracy
{
	Matrix r(1), ao=a;
	for (int i=0; i<VectorDimension; i++)
	{
		//pivoting
      {
		n_type fmax=-1; int jmax;
		for (int j=i; j<VectorDimension; j++)
			if ( fmax<abs(a.x[j][i]) )
         {jmax=j; fmax=abs(a.x[j][i]);}
		ComplexType * aux;
		aux=a.x[i]; a.x[i]=a.x[jmax]; a.x[jmax]=aux;
		aux=r.x[i]; r.x[i]=r.x[jmax]; r.x[jmax]=aux;
		if (norm2(a.x[i][i])<=0) {cout<<"LInverse!";SetNull(a); return;}
      ;}
		//main body
		for (int j=0; j<VectorDimension; j++)
		if (j!=i)
		{
			ComplexType f=a.x[j][i]/a.x[i][i];
			{for (int l=i; l<VectorDimension; l++) a.x[j][l]-=a.x[i][l]*f;}
			{for (int l=0; l<VectorDimension; l++) r.x[j][l]-=r.x[i][l]*f;}
		;}
		else
		{
			ComplexType f=n_type(1)/a.x[i][i];
			{for (int l=i; l<VectorDimension; l++) a.x[j][l]*=f;}
			{for (int l=0; l<VectorDimension; l++) r.x[j][l]*=f;}
		;}
	;}
	if (improve==0) a=r;
   else {Matrix d(0); d=ao*r; a=2*r-r*d;}
;}

/*

Matrix & operator / (Matrix & A, Matrix & B)
{
	//Matrix * z=new Matrix();
	Matrix C=B; Inverse(C);
	return A*C;
;}

Matrix & operator / (ComplexType &A, Matrix & B)
{
	Matrix C=B; Inverse(C);
	return A*C;
;}

Matrix operator / (n_type A, Matrix & B)
{
	Matrix C=B; Inverse(C);
	return A*C;
;}
*/

ComplexType det(Matrix &a)
{
	Matrix b=a;
	ComplexType s=1;
	for (int i=0; i<VectorDimension-1; i++)
	{
		//pivoting
      {
		n_type fmax=0; int jmax=i;
		for (int j=i; j<VectorDimension; j++)
		if ( fmax<abs(b.x[j][i]) ) {jmax=j; fmax=abs(b.x[j][i]);}
		if (jmax!=i)
		{
			ComplexType * aux;
			aux=b.x[i]; b.x[i]=b.x[jmax]; b.x[jmax]=aux;
			s=-s;
		;}
      ;}
		//main body
		if (norm2(b.x[i][i])<=0) return 0;
		for (int j=i+1; j<VectorDimension; j++)
		{
		ComplexType f=b.x[j][i]/b.x[i][i];
		for (int k=i; k<VectorDimension; k++) b.x[j][k]-=f*b.x[i][k];
		;}
	;}


   if (norm2(b.x[VectorDimension-1][VectorDimension-1])<=0) return 0;
	{for (int i=0; i<VectorDimension; i++) s*=b.x[i][i];}
   return s;
;}



