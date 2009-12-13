//Import of the interaction parameters, Delta, etc. from files

Matrix & Rotate (int z)                              //no z-dependence for a while....
{
	static Matrix A; static int f=0;
   if (f==0)
   {
    int i=int_value("rotate_basis");
    if (i!=0)
    {
      ifstream rot("rotate.dat");
   	for (int i=0; i<n_part; i++)
   	for (int j=0; j<n_part; j++)
   		{double x; rot>>x; A.x[i][j]=x;}
    ;}
    else A=n_type(1);
    f=1;
   ;}
   return A;
;}


void read_Delta()
{
   int read_delta=int_value("read_delta");

   if (read_delta<2)
   for (int z=0; z<n_zone; z++)
   {
		ifstream delta("Delta.dat");
   	Delta[z]=new Matrix [wn_max];
      for (int wn=0; wn<wn_max; wn++)
      {
      	(Delta[z][wn]).new_memory();
         for (int i=0; i<n_part; i++)
//         for (int j=0; j<n_part; j++)
         {
           	n_type er=0, ei=0;
            if (read_delta==1) {delta >>er; delta >>ei;}

            (Delta[z][wn]).x[i][i]=ComplexType(er,ei);

         ;}
      ;}
   ;}
   else
   {
			ifstream delta("Delta.dat");
         {
         	for (int z=0; z<n_zone; z++)
            {
            	Delta[z]=new Matrix [wn_max];
               for (int wn=0; wn<wn_max; wn++)  (Delta[z][wn]).new_memory();
            ;}
         ;}
         for (int wn=0; wn<wn_max; wn++)
{
         for (int z=0; z<n_zone; z++)
         for (int i=0; i<n_part; i++)
         for (int j=0; j<n_part; j++)
         {
         	n_type er=0, ei=0; delta>>er; delta>>ei;
            (Delta[z][wn]).x[i][j]=ComplexType(er,ei);          //cout<<(Delta[z][wn]).x[i][j]<<"  ";
         ;} //cout<<"\n";
;}

   ;}

   //rotate Delta

   {
   	for (int z=0; z<n_zone; z++)
      {
         Matrix l=Rotate(z), r(l); r.trans();
      	for (int wn=0; wn<wn_max; wn++)
      		Delta[z][wn]=r*Delta[z][wn]*l;//from p- to r- basis 
      ;}
   ;}

;}




int number_of_fields=int_value("number_of_fields");
int * interaction_indices [6];
n_type * interaction_U;  n_type interaction_U_norm=0;


void set_interaction()
{
	{static int f=0; if (f==1) return; f=1;}
   {for (int i=0; i<6; i++) interaction_indices[i]=new int [number_of_fields];}
   interaction_U=new n_type [number_of_fields];
	ifstream interaction("interaction.dat");
   for (int j=0; j<number_of_fields; j++)
   {
   	for (int i=0; i<6; i++) interaction>>interaction_indices[i][j];
      interaction>>interaction_U[j]; interaction_U_norm+=fabs(interaction_U[j]);
   ;}
;}

void Import_Interaction (point & r1, point_ & r1_, point & r2, point_ & r2_, n_type &u, n_type & a1, n_type &a2)
{
		static n_type alpha=value("alpha");
static n_type alpha_off=value("alpha_offdiagonal");
	set_interaction();
	   n_type t=beta*rnd();
   r1.t=t; r1_.t=t; r2.t=t; r2_.t=t;


	n_type r=rnd(); int j=-1;
   do {j++;r-=fabs(interaction_U[j])/interaction_U_norm;} while(r>0 && j<number_of_fields-1);
   //j is the number of field component...

   if (rnd()<0.5)
   {
   	r1.z=interaction_indices[0][j];
   	r1_.z=interaction_indices[0][j];
   	r2.z=interaction_indices[3][j];
   	r2_.z=interaction_indices[3][j];

   	r1.i=interaction_indices[1][j];
   	r1_.i=interaction_indices[2][j];
   	r2.i=interaction_indices[4][j];
   	r2_.i=interaction_indices[5][j];
   ;}
   else //1 <-> 2 symmetrization
   {
   	r2.z=interaction_indices[0][j];
   	r2_.z=interaction_indices[0][j];
   	r1.z=interaction_indices[3][j];
   	r1_.z=interaction_indices[3][j];

   	r2.i=interaction_indices[1][j];
   	r2_.i=interaction_indices[2][j];
   	r1.i=interaction_indices[4][j];
   	r1_.i=interaction_indices[5][j];
   ;}

   u=beta*interaction_U_norm;
   if (interaction_U[j]<0) u=-u;

   if (r1.i==r1_.i && r2.i==r2_.i)
   {//diagonal
   	if (u<0) {a1=-alpha; a2=-alpha;}
      else
      	{if (rnd()<0.5) {a1=alpha; a2=1-alpha;} else  {a1=1-alpha; a2=alpha;};}
   ;}
   else
   {//off-diagonal
   	a1=alpha_off; a2=alpha_off;
   ;}

;}


