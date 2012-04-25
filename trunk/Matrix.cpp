#include "Matrix.h"
#include <fstream>

n_type norm2(ComplexType & x)
{
    return x.real()*x.real()+x.imag()*x.imag();
}

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



void Inverse(Matrix & a, int improve) //Gauss with partial pivoting
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



