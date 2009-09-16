 //================= separate subroutines for 4-point correlators =======================
complex **** GMchi_matrix[n_zone];
void set_GMchi();
void ini_GMchi();
                            

complex *** chi2_array [n_zone];
complex ******* chi4_array [n_zone][n_zone];               //w1,w2,W ->       (w1 w1+W; w2+W w2)
complex *********** chi6_array [n_zone][n_zone][n_zone];   //w0, w0_, w1, w1_, w2 
void ini_chi2_array();
void ini_chi4_array();
void ini_chi6_array();

complex chi2_extract(int z, int z_, int w, int w_, int n, int n_, double s)
{
	if (z!=z_ || w!=w_) return 0;
   return chi2_array[z][w][n][n_]/s;
;}

complex chi4_irreducible(int z1, int z2, int w1, int w2_, int W, int n1, int n1_, int n2, int n2_, double s)
{
	 if (W<0) {static int wm=2*(int_value("number_of_Matsubara_frequencies_for_Gamma4")/2);  return conj(chi4_irreducible(z1,z2,wm-1-w1,wm-1-w2_,-W,n1,n1_,n2,n2_,s));}
    int w1_=w1+W, w2=w2_+W;
    complex C=chi4_array[z1][z2][w1][w2_][W][n1][n1_][n2][n2_]/s;
    C-=chi2_extract(z1, z1, w1, w1_, n1, n1_,s)*chi2_extract(z2, z2, w2, w2_, n2, n2_,s);
    C+=chi2_extract(z1, z2, w1, w2_, n1, n2_,s)*chi2_extract(z2, z1, w2, w1_, n2, n1_,s);
    return C;
;}

complex chi4_extract(int z1, int z1_, int z2, int z2_, int w1, int w1_, int w2, int w2_, int n1, int n1_, int n2, int n2_, double s)
{
	if (w1+w2-w1_-w2_!=0) return 0;
   if (z1==z1_ && z2==z2_) return chi4_irreducible(z1, z2, w1, w2_, w1_-w1, n1, n1_, n2, n2_,s);
   if (z1==z2_ && z2==z1_) return -chi4_extract(z1, z2_, z2, z1_, w1, w2_, w2, w1_, n1, n2_, n2, n1_,s);
   return 0;

;}



complex chi6_irreducible(int z0,int z1,int z2,int w0, int w0_, int w1, int w1_, int w2, int n0, int n0_, int n1, int n1_, int n2, int n2_, double s)
{
static int wm=2*(int_value("number_of_Matsubara_frequencies_for_Gamma4")/2);
if (w0<0) {return conj(chi6_irreducible(z0, z1, z2, wm-1-w0, wm-1-w0_, wm-1-w1, wm-1-w1_, wm-1-w2, n0, n0_, n1, n1_, n2, n2_, s));}

int w2_=w2+w1+w0-w1_-w0_;
	complex C=chi6_array[z0][z1][z2][w0][w0_][w1][w1_][w2][n0][n0_][n1][n1_][n2][n2_]/s;
   // 9 terms 2-4
   			 C-=chi2_extract(z0,z0,w0,w0_,n0,n0_,s)*chi4_extract(z1,z1,z2,z2,w1,w1_,w2,w2_,n1,n1_,n2,n2_,s);
             C+=chi2_extract(z0,z1,w0,w1_,n0,n1_,s)*chi4_extract(z1,z0,z2,z2,w1,w0_,w2,w2_,n1,n0_,n2,n2_,s);
             C+=chi2_extract(z0,z2,w0,w2_,n0,n2_,s)*chi4_extract(z1,z1,z2,z0,w1,w1_,w2,w0_,n1,n1_,n2,n0_,s);

   			 C-=chi2_extract(z1,z1,w1,w1_,n1,n1_,s)*chi4_extract(z0,z0,z2,z2,w0,w0_,w2,w2_,n0,n0_,n2,n2_,s);
             C+=chi2_extract(z1,z0,w1,w0_,n1,n0_,s)*chi4_extract(z0,z1,z2,z2,w0,w1_,w2,w2_,n0,n1_,n2,n2_,s);
             C+=chi2_extract(z1,z2,w1,w2_,n1,n2_,s)*chi4_extract(z0,z0,z2,z1,w0,w0_,w2,w1_,n0,n0_,n2,n1_,s);



    			 C-=chi2_extract(z2,z2,w2,w2_,n2,n2_,s)*chi4_extract(z0,z0,z1,z1,w0,w0_,w1,w1_,n0,n0_,n1,n1_,s);
             C+=chi2_extract(z2,z0,w2,w0_,n2,n0_,s)*chi4_extract(z0,z2,z1,z1,w0,w2_,w1,w1_,n0,n2_,n1,n1_,s);
             C+=chi2_extract(z2,z1,w2,w1_,n2,n1_,s)*chi4_extract(z0,z0,z1,z2,w0,w0_,w1,w2_,n0,n0_,n1,n2_,s);


   //6 terms 2-2-2
   C-=chi2_extract(z0,z0,w0,w0_,n0,n0_,s)*chi2_extract(z1,z1,w1,w1_,n1,n1_,s)*chi2_extract(z2,z2,w2,w2_,n2,n2_,s);

   C+=chi2_extract(z0,z0,w0,w0_,n0,n0_,s)*chi2_extract(z1,z2,w1,w2_,n1,n2_,s)*chi2_extract(z2,z1,w2,w1_,n2,n1_,s);
   C+=chi2_extract(z0,z1,w0,w1_,n0,n1_,s)*chi2_extract(z1,z0,w1,w0_,n1,n0_,s)*chi2_extract(z2,z2,w2,w2_,n2,n2_,s);
   C+=chi2_extract(z0,z2,w0,w2_,n0,n2_,s)*chi2_extract(z1,z1,w1,w1_,n1,n1_,s)*chi2_extract(z2,z0,w2,w0_,n2,n0_,s);

	C-=chi2_extract(z0,z2,w0,w2_,n0,n2_,s)*chi2_extract(z1,z0,w1,w0_,n1,n1_,s)*chi2_extract(z2,z1,w2,w1_,n2,n1_,s);
   C-=chi2_extract(z0,z1,w0,w1_,n0,n1_,s)*chi2_extract(z1,z2,w1,w2_,n1,n2_,s)*chi2_extract(z2,z0,w2,w0_,n2,n0_,s);

   return C;
;}

void write_chi4(double);
void write_chi6(double);

//===============   main body  ================

void chi(int write_flag=0)
 {
 	static int todo=int_value("calculate_Gamma4"); if (todo==0) return;
   static int todo6=int_value("calculate_Gamma6");
   static int wc_max=int_value("number_of_Matsubara_frequencies_for_Gamma4");

	static double s=1e-100;
   static int count=0;



  // static int f=1;
   ini_GMchi(); ini_chi2_array(); ini_chi4_array(); if (todo6>0) ini_chi6_array();

// ======================== write to the files =============================
 	if (write_flag==1)
   {
   write_chi4(s);
   if (todo6!=0) write_chi6(s);

   return;
   ;}


// =======================   calcultaion procedure ==================
count++; if (count%(n_part*wn_max*(1+31*todo6))!=0) return;   //to mantain N^2 complexity for chi4...
set_GMchi();
s+=curr_sgn/prev_weight;



{
		for (int z=0; z<n_zone; z++)
      for (int w=0; w<2*wc_max; w++)
      for (int n1=0; n1<n_part; n1++)
      for (int n2=0; n2<n_part; n2++)
      	chi2_array[z][w][n1][n2]+=GMchi_matrix[z][w][w][n1][n2]*(curr_sgn/prev_weight);

;}

{
      for (int z1=0; z1<n_zone; z1++)
      for (int z2=0; z2<n_zone; z2++)
      for (int n1=0; n1<n_part; n1++)
      for (int n1_=0; n1_<n_part; n1_++)
      for (int n2=0; n2<n_part; n2++)
      for (int n2_=0; n2_<n_part; n2_++)
      for (int w1=0; w1<2*wc_max/2; w1++)
      for (int w2_=0; w2_<2*wc_max/2; w2_++)
      for (int w1_=w1; w1_<2*wc_max/2; w1_++)
//      for (int w2=0; w2<2*wc_max/2; w2++) if (w1-w1_+w2-w2_==0 && w1_-w1>=0 && w1_-w1<2*(wc_max/2))
      {  int w2=w1_+w2_-w1;
      	if (w2>=0 && w2<2*wc_max/2)
         {
      	int W=w1_-w1; double f=(curr_sgn/prev_weight);
      	chi4_array[z1][z2][w1][w2_][W][n1][n1_][n2][n2_]+=
              	GMchi_matrix[z1][w1][w1_][n1][n1_]*GMchi_matrix[z2][w2][w2_][n2][n2_]*f;
         if (z1==z2)
         chi4_array[z1][z2][w1][w2_][W][n1][n1_][n2][n2_]-=
         		GMchi_matrix[z1][w1][w2_][n1][n2_]*GMchi_matrix[z2][w2][w1_][n2][n1_]*f;
         ;}
      ;}

;}
if (todo6==0) return;

{
	for (int z0=0; z0<n_zone; z0++)
   for (int z1=0; z1<n_zone; z1++)
   for (int z2=0; z2<n_zone; z2++)
   for (int n0_=0; n0_<n_part; n0_++)
   for (int n0=0; n0<n_part; n0++)
   for (int n1_=0; n1_<n_part; n1_++)
   for (int n1=0; n1<n_part; n1++)
   for (int n2_=0; n2_<n_part; n2_++)
   for (int n2=0; n2<n_part; n2++)
   for (int w0=0; w0<2*wc_max/2; w0++)
   for (int w0_=0; w0_<2*wc_max/2; w0_++)
   for (int w1=0; w1<2*wc_max/2; w1++)
   for (int w1_=0; w1_<2*wc_max/2; w1_++)
   for (int w2=0; w2<2*wc_max/2; w2++)
   {  int w2_=w0+w1+w2-w0_-w1_;
     	if (w2_>=0 && w2_<2*wc_max/2)
      {
   	double f=(curr_sgn/prev_weight); complex C=0.;
                   C+=GMchi_matrix[z0][w0][w0_][n0][n0_]*GMchi_matrix[z1][w1][w1_][n1][n1_]*GMchi_matrix[z2][w2][w2_][n2][n2_]*f;
      if (z1==z2)  C-=GMchi_matrix[z0][w0][w0_][n0][n0_]*GMchi_matrix[z1][w1][w2_][n1][n2_]*GMchi_matrix[z2][w2][w1_][n2][n1_]*f;
      if (z1==z0)  C-=GMchi_matrix[z0][w0][w1_][n0][n1_]*GMchi_matrix[z1][w1][w0_][n1][n0_]*GMchi_matrix[z2][w2][w2_][n2][n2_]*f;
      if (z2==z0)  C-=GMchi_matrix[z0][w0][w2_][n0][n2_]*GMchi_matrix[z1][w1][w1_][n1][n1_]*GMchi_matrix[z2][w2][w0_][n2][n0_]*f;
      if (z0==z1 && z0==z2) C+=GMchi_matrix[z0][w0][w2_][n0][n2_]*GMchi_matrix[z1][w1][w0_][n1][n0_]*GMchi_matrix[z2][w2][w1_][n2][n1_]*f;
      if (z0==z1 && z0==z2) C+=GMchi_matrix[z0][w0][w1_][n0][n1_]*GMchi_matrix[z1][w1][w2_][n1][n2_]*GMchi_matrix[z2][w2][w0_][n2][n0_]*f;


      chi6_array[z0][z1][z2][w0][w0_][w1][w1_][w2][n0][n0_][n1][n1_][n2][n2_]+=C;
      ;}
   ;}
;}

;}





















void ini_GMchi()
{
   static int f=1;
   static int wc_max=int_value("number_of_Matsubara_frequencies_for_Gamma4");
   if (f==1)
   {
   for (int z=0; z<n_zone; z++)
   {
   GMchi_matrix[z]=new complex *** [2*wc_max];
   	for (int w1=0; w1<2*wc_max; w1++)
   	{
      	GMchi_matrix[z][w1]=new complex ** [2*wc_max];
         for (int w2=0; w2<2*wc_max; w2++)
         {
         	GMchi_matrix[z][w1][w2]=new complex * [n_part];
            for (int n1=0; n1<n_part; n1++)
            	GMchi_matrix[z][w1][w2][n1]=new complex [n_part];
         ;}
   	;}
   ;}
   ;}
   f=0;
;}








void set_GMchi()
 {
 	static int wc_max=int_value("number_of_Matsubara_frequencies_for_Gamma4");
   static complex *** MR [n_zone];
   static int K_max[n_zone]; {for (int z=0; z<n_zone; z++) K_max[z]=-1;}

   static int f=0; if (f==0)
   {
   	for (int z=0; z<n_zone; z++)
      {
      	MR[z]=new complex ** [2*wc_max];
         for (int w=0; w<2*wc_max; w++)
         {
         	MR[z][w]=new complex * [n_part];
            for (int n=0; n<n_part; n++)
            	MR[z][w][n]=new complex [1];
         ;}
      ;}
   	f=1;
   ;}

   {
   for (int z=0; z<n_zone; z++)
   if (Nc_z[z]>K_max[z])
   {
      for (int w=0; w<2*wc_max; w++)
      for (int n=0; n<n_part; n++)
      {
      	delete MR[z][w][n];
         MR[z][w][n]=new complex [Nc_z[z]];
      ;}
   K_max[z]=Nc_z[z];
   ;}
   ;}

 	for (int z=0; z<n_zone; z++)
   {
   	{
   	for (int i=0; i<Nc_z[z]; i++)
      for (int w=0; w<wc_max; w++)
      for (int n=0;n<n_part;n++)
      {
      	int W2=w-wc_max/2;
      	MR[z][w][n][i]=0;
         static Matrix grot(0); if (W2>=0) grot=(*Grot(z,W2,-1)); else grot=(*Grot(z,-W2-1,-1));
         int kk=-1;
         for (int j=0;j<Nc_z[z];j++)
         {
         	next(z,kk);
            complex EW; if (W2>=0) EW=ewt[kk][W2]*grot.x[n][p_[kk].i]; else EW=conj(ewt[kk][-W2-1]*grot.x[n][p_[kk].i]);
	        	MR[z][w][n][i]+=EW*M[z][j][i];//ewt[kk][w]
         ;}
      ;}
      ;}

      {
      for (int w=0; w<wc_max; w++)
      for (int n=0;n<n_part;n++)
      for (int w2=0; w2<wc_max; w2++)
      for (int n2=0;n2<n_part;n2++)
      {
      	int W2=w2-wc_max/2;
      	GMchi_matrix[z][w][w2][n][n2]=0;
         static Matrix grot(0); if (W2>=0) grot=(*Grot(z,W2,1)); else grot=(*Grot(z,-W2-1,1));
         if (w==w2 && W2>=0) GMchi_matrix[z][w][w2][n][n2]=(*Grot(z,W2,0)).x[n][n2];
         if (w==w2 && W2<0) GMchi_matrix[z][w][w2][n][n2]=conj((*Grot(z,-W2-1,0)).x[n][n2]);
         int kk=-1;
         for (int j=0;j<Nc_z[z];j++)
         {
         	next(z,kk);
            complex EW; if (W2>=0) EW=ewt_[kk][W2]*grot.x[p[kk].i][n2]; else EW=conj(ewt_[kk][-W2-1]*grot.x[p[kk].i][n2]);
	        	GMchi_matrix[z][w][w2][n][n2]-=MR[z][w][n][j]*EW/beta;//ewt_[kk][w2]
         ;}
      ;}
      ;}




   ;}
;}

//================== Write to the files ==============================

void write_chi4(double s)
{
	double S4=0;
   ofstream chi_f("Gamma4.dat");
   chi_f<<"Re         Im               z1 z2          w1' w1 w2' w2        n1' n1 n2' n2\n";
   static int wc_max=int_value("number_of_Matsubara_frequencies_for_Gamma4");

      for (int z1=0; z1<n_zone; z1++)
      for (int z2=n_zone-1; z2>=0; z2--)
      for (int n1=0; n1<n_part; n1++)
      for (int n1_=0; n1_<n_part; n1_++)
      for (int n2=0; n2<n_part; n2++)
      for (int n2_=0; n2_<n_part; n2_++)
      {

//      for (int W=0; W<2*wc_max/2; W++)
      {
      for (int w1=0; w1<2*wc_max/2; w1++)
      for (int w1_=0; w1_<2*wc_max/2; w1_++)
      for (int w2=0; w2<2*wc_max/2; w2++)
      for (int w2_=0; w2_<2*wc_max/2; w2_++)   if (w1-w1_+w2-w2_==0 && abs(w1_-w1)<2*(wc_max/2))
      {
      	int W=w1_-w1;
      	complex z=chi4_irreducible(z1, z2, w1, w2_, W,n1,n1_,n2,n2_,s);
         z/=chi2_array[z1][w1][n1][n1]/s;                       //true only for a diagonal basis!!!!!
         z/=chi2_array[z1][w1+W][n1_][n1_]/s;
         z/=chi2_array[z2][w2_+W][n2][n2]/s;
         z/=chi2_array[z2][w2_][n2_][n2_]/s;
         z*=beta;

	   static n_type acc=value("chi4_numerical_zero");
	   if (abs(z)>acc)
	   {		
         	chi_f<<real(z)<<" "<<imag(z)<<"        ";
         	chi_f<<z1<<" "<<z2<<"              "<<w1-wc_max/2<<" "<<w1_-wc_max/2<<" "<<w2-wc_max/2<<" "<<w2_-wc_max/2<<"            "<<n1<<" "<<n1_<<" "<<n2<<" "<<n2_<<"\n";
         ;} 
         S4+=real(z)*real(z)+imag(z)*imag(z);
      ;}

      ;}
      ;}
chi_f<<"0    0";
//          debug<<S4<<"  ";
;}

void write_chi6(double s)
{ double S6=0;
	ofstream chi_f6("Gamma6.dat");
   chi_f6<<"Re        Im      z1   z2   z3   w0'  w0   w1'  w1   w2'  w2     n0'   n0   n1'   n1   n2'   n2\n";
   static int wc_max=int_value("number_of_Matsubara_frequencies_for_Gamma4");
   for (int z0=0; z0<n_zone; z0++)
   for (int z1=0; z1<n_zone; z1++)
   for (int z2=0; z2<n_zone; z2++)
   for (int n0_=0; n0_<n_part; n0_++)
   for (int n0=0; n0<n_part; n0++)
   for (int n1_=0; n1_<n_part; n1_++)
   for (int n1=0; n1<n_part; n1++)
   for (int n2_=0; n2_<n_part; n2_++)
   for (int n2=0; n2<n_part; n2++)

   for (int w0=0; w0<2*wc_max/2; w0++)
   for (int w0_=0; w0_<2*wc_max/2; w0_++)
   for (int w1=0; w1<2*wc_max/2; w1++)
   for (int w1_=0; w1_<2*wc_max/2; w1_++)
   for (int w2=0; w2<2*wc_max/2; w2++)
   for (int w2_=0; w2_<2*wc_max/2; w2_++)  if (w0-w0_+w1-w1_+w2-w2_==0)
   {
      	int WM=wc_max/2;
      	complex z=chi6_irreducible(z0, z1, z2, w0, w0_, w1, w1_, w2, n0, n0_, n1, n1_, n2, n2_, s);

         z/=chi2_array[z0][w0][n0][n0]/s;        	               //true only for a diagonal basis!!!!!
         z/=chi2_array[z0][w0_][n0_][n0_]/s;
         z/=chi2_array[z1][w1][n1][n1]/s;
         z/=chi2_array[z1][w1_][n1_][n1_]/s;
         z/=chi2_array[z2][w2][n2][n2]/s;
         z/=chi2_array[z2][w2_][n2_][n2_]/s;
         z*=beta*beta;

         static n_type acc=value("chi6_numerical_zero");
         if (abs(z)>acc)
         {
           chi_f6<<real(z)<<" "<<imag(z)<<" "<<z0<<" "<<z1<<" "<<z2<<" "<<w0-WM<<" "<<w0_-WM<<" "<<w1-WM<<" "<<w1_-WM<<" "<<w2-WM<<" "<<w2_-WM<<" ";
           chi_f6<<n0<<" "<<n0_<<" "<<n1<<" "<<n1_<<" "<<n2<<" "<<n2_<<"\n";
         ;}
         S6+=real(z)*real(z)+imag(z)*imag(z);
   ;}


   chi_f6<<"0    0";
    //debug<<S6<<"\n"<<flush;
;}


//================== Initialisation of arrays ========================






void ini_chi2_array()
{
	static int f=0; if (f==1) return;
   int wc_max=int_value("number_of_Matsubara_frequencies_for_Gamma4");
   for (int z=0; z<n_zone; z++)
   {
   chi2_array[z]=new complex ** [2*wc_max];
   	for (int w1=0; w1<2*wc_max; w1++)
   	{
         chi2_array[z][w1]=new complex * [n_part];
            for (int n=0; n<n_part; n++)
            {
            	chi2_array[z][w1][n]= new complex [n_part];
               for (int n2=0; n2<n_part; n2++)  chi2_array[z][w1][n][n2]=0;
            ;}
   	;}
   ;}
   f=1;
;}



void ini_chi4_array()
{
	static int f=0; if (f==1) return;
   int wc_max=int_value("number_of_Matsubara_frequencies_for_Gamma4");
   for (int z1=0; z1<n_zone; z1++)
   for (int z2=0; z2<n_zone; z2++)
   {
   	chi4_array[z1][z2]=new complex ******[2*wc_max/2];
      for (int w1=0; w1<2*wc_max/2; w1++)
      {
      	chi4_array[z1][z2][w1]=new complex *****[2*wc_max/2];
         for (int w2_=0; w2_<2*wc_max/2; w2_++)
         {
	         chi4_array[z1][z2][w1][w2_]=new complex ****[2*wc_max/2];
            for (int W=0; W<2*wc_max/2; W++)
            {
            	chi4_array[z1][z2][w1][w2_][W]=new complex ***[n_part];
               for (int n1=0; n1<n_part; n1++)
               {
                  chi4_array[z1][z2][w1][w2_][W][n1]=new complex **[n_part];
               	for (int n1_=0; n1_<n_part; n1_++)
                  {
                  	chi4_array[z1][z2][w1][w2_][W][n1][n1_]=new complex * [n_part];
                     for (int n2=0; n2<n_part; n2++)
                     {
                     	chi4_array[z1][z2][w1][w2_][W][n1][n1_][n2]=new complex [n_part];
               			for (int n2_=0; n2_<n_part; n2_++)  chi4_array[z1][z2][w1][w2_][W][n1][n1_][n2][n2_]=0
                     ;}
                  ;}
               ;}
            ;}
         ;}
      ;}
   ;}


   f=1;
;}


void ini_chi6_array()
{
	static int f=0; if (f==1) return;
   int wc_max=int_value("number_of_Matsubara_frequencies_for_Gamma4");
   for (int z1=0; z1<n_zone; z1++)
   for (int z2=0; z2<n_zone; z2++)
   for (int z3=0; z3<n_zone; z3++)
   {
   	chi6_array[z1][z2][z3]=new complex **********[2*wc_max/2];
      for (int w0=0; w0<2*wc_max/2; w0++)
      {
        chi6_array[z1][z2][z3][w0]=new complex *********[2*wc_max/2];
        for (int w0_=0; w0_<2*wc_max/2; w0_++)
        {
      	chi6_array[z1][z2][z3][w0][w0_]=new complex ********[2*wc_max/2];
         for (int w1=0; w1<2*wc_max/2; w1++)
         {
	         chi6_array[z1][z2][z3][w0][w0_][w1]=new complex *******[2*wc_max/2];
            for (int w1_=0; w1_<2*wc_max/2; w1_++)
            {
            	chi6_array[z1][z2][z3][w0][w0_][w1][w1_]=new complex ******[2*wc_max/2];
               for (int w2=0; w2<2*wc_max/2; w2++)
               {
                 chi6_array[z1][z2][z3][w0][w0_][w1][w1_][w2]=new complex *****[n_part];
                 for (int n1=0; n1<n_part; n1++)
                 {
                   chi6_array[z1][z2][z3][w0][w0_][w1][w1_][w2][n1]=new complex ****[n_part];
               	 for (int n1_=0; n1_<n_part; n1_++)
                   {
                  	chi6_array[z1][z2][z3][w0][w0_][w1][w1_][w2][n1][n1_]=new complex *** [n_part];
                     for (int n2=0; n2<n_part; n2++)
                     {
                     	chi6_array[z1][z2][z3][w0][w0_][w1][w1_][w2][n1][n1_][n2]=new complex **[n_part];
                        for (int n2_=0; n2_<n_part; n2_++)
                        {
                        	chi6_array[z1][z2][z3][w0][w0_][w1][w1_][w2][n1][n1_][n2][n2_]=new complex *[n_part];
                           for (int n3=0; n3<n_part; n3++)
                           {
                           	chi6_array[z1][z2][z3][w0][w0_][w1][w1_][w2][n1][n1_][n2][n2_][n3]=new complex [n_part];
                        		for (int n3_=0; n3_<n_part; n3_++) 	chi6_array[z1][z2][z3][w0][w0_][w1][w1_][w2][n1][n1_][n2][n2_][n3][n3_]=0;
                           ;}
                        ;}
                     ;}
                   ;}
                 ;}
               ;}
            ;}
         ;}
        ;}
      ;}
   ;}


   f=1;
;}


