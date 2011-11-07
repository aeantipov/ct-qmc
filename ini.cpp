
extern default_dict internal_input;
#ifdef use_mpi
extern CTQMC_WORLD CTQMC;
#endif

void Ini()
{

#ifdef use_mpi
    CTQMC.getStream()<<"\n       This is CT-QMC code version 1.6\n       Check www.ct-qmc.org for latest updates.\n";
    CTQMC.getStream()<<"\n       Using MPI Version.\n";
#else 
    cout<<"\n       This is CT-QMC code version 1.6\n       Check www.ct-qmc.org for latest updates.\n";
    cout<<"\n       Using single-processor version.\n";
#endif

   try{
   for (int ll=0;ll<N_max;ll++){point p0(0,0,0); p[ll]=p0; point_ p0_(0,0,0); p_[ll]=p0_;}
	mu1_Gt=new n_type * [n_zone]; {for (int i=0; i<n_zone; i++) mu1_Gt[i]=new n_type [n_part];}
	mu2_Gt=new n_type * [n_zone]; {for (int i=0; i<n_zone; i++) mu2_Gt[i]=new n_type [n_part];}   

   {
   for (int z=0; z<n_zone; z++)
   {
   	NumberOfPoints[z]=0;
   	M[z]=new n_type * [N_max];
      for (int j=0; j<N_max; j++)
      {
      	M[z][j]=new n_type [N_max];
      ;}                                               
   GM_matrix[z]=new ComplexType **[wn_max];
   Gtotal_st[z]=new ComplexType **[wn_max];
   GM00_st[z]=new ComplexType *[wn_max];
   for (int w=0; w<wn_max; w++)
   {
      	GM_matrix[z][w]=new ComplexType *[n_part];
         Gtotal_st[z][w]=new ComplexType *[n_part];
         GM00_st[z][w]=&dataGM00[z*wn_max*n_part + w*n_part];
         for (int n1=0; n1<n_part; n1++)
         {
         	GM_matrix[z][w][n1]=new ComplexType [n_part];
            Gtotal_st[z][w][n1]=&dataGtotal[z*wn_max*n_part*n_part + w*n_part*n_part + n1*n_part];
            for (int n2=0; n2<n_part; n2++) 
         ;}
   ;}

   ;}
   ;}
   }
   catch (std::bad_alloc&)
   {
     #ifdef use_mpi
     cout << "Process " << CTQMC.getRank() << " couldn't allocate memory" << endl;
	 #else
     cout << "Couldn't allocate memory" << endl;
	 #endif
   };

   read_Delta();

   cluster_ini();

   {for (int i=0; i<N_max; i++) {ewt[i]=new ComplexType [wn_max]; ewt_[i]=new ComplexType [wn_max];};}


   {for (int i=0; i<N_max; i++) hist[i]=0;}


	{for (int i=0; i<N_max; i++) WL_stream>>WL_weight[i];}

WL_max=abs(int(n_part*U*beta*value("WL_factor"))*2)+2;

   from_scratch();

   #ifdef use_mpi
   DEBUG(CTQMC.getRank() << " has made it from scratch. ");
   #else 
   cout<<"+"<<flush;
   #endif
if (WL_flag==0)
{
	int n=WL_max*WL_max; if (n>int_value("Initial_steps_without_WL")) n=int_value("Initial_steps_without_WL");
	DEBUG("Doing " << n << " steps " << endl);
   {
   	for (int i=0; i<n; i++) {
	DEBUG("step : " << i << " Nc :  " << Nc << endl);
	DEBUG("NumberOfPoints[] : ");
	for (int z=0;z<n_zone;z++) DEBUG(NumberOfPoints[z] << " "); 
	steps();
	DEBUG(endl);
	DEBUG("After steps: Nc = " << Nc << " NumberOfPoints[]:  ")
	for (int z=0;z<n_zone;z++) DEBUG(NumberOfPoints[z] << " "); 
	DEBUG(endl);
/*	try {steps();}
	catch (std::bad_alloc&)
	{ cout << "AAAAAAAAAAA " << endl;
	  exit(0);
	}
	catch (std::bad_exception)
	{ cout << "Some unknown shit happened" << endl;
	  exit(0);
	};*/
	if (i%1000==0) 
   	{
   #ifdef use_mpi
 	  CTQMC.getStream()<<">"<<flush;
   #else 
 	  cout<<">"<<flush;
   #endif
	};}}

   #ifdef use_mpi
	DEBUG(CTQMC.getRank() << " has done steps" << endl);
	#endif
   from_scratch(); 
   STREAM <<"+"<<flush;
;}

if  (WL_flag==0 || WL_flag==2)  {for (int i=0; i<N_max; i++) WL_weight[i]=0;}
if  (WL_flag==1 || WL_flag==2)
{

WL_factor=value ("WL_initial_factor");
int lmax=int_value("WL_circles_number");
n_type WL_Tolerance=value("WL_Tolerance_factor");
ofstream WL_str("WL_hist.dat");
//ofstream deb("debug.dat");
for (int l=0; l<lmax; l++)
{

int flag; {for (int i=0; i<N_max; i++) hist[i]=0;}
do
{
   steps();
   flag =1;
   if (l==0) {for (int i=0; i<WL_max/2; i++) if (hist[i]<100) flag=0;}
   else
   {
   	n_type hmax=0, hmin=1e6;
      for (int i=0; i<WL_max/2; i++) {if (hist[i]<hmin) hmin=hist[i];if (hist[i]>hmax) hmax=hist[i];}
      if (hmin<WL_Tolerance*hmax) flag=0;
   ;}
;}
while (flag==0); 

#ifdef use_mpi
CTQMC.getStream()<<"+"<<flush;
#else
cout<<"+"<<flush;
#endif
for (int i=0; i<WL_max/2; i++) WL_str<<hist[i]<<"       "; WL_str<<"\n"<<flush;
if (l!=0) WL_factor/=2;    from_scratch();
;}
#ifdef use_mpi
CTQMC.getStream()<<"WL factors are determined up to "<<WL_max/2<<flush;
#else
cout<<"WL factors are determined up to "<<WL_max/2<<flush;
#endif
;}

   weight_sum=0; weight_sum2=0;
	sgn_sum=0;
	Nc_sum=0;
	DispG0_sum=0;
   GLOBAL_MC_COUNT=0, ACCEPT_MC_COUNT=0;

   WL_factor=0;
{for (int i=0; i<N_max; i++) hist[i]=0;}
//{for (int i=0;i<100000;i++) MC_step();}


   int max_j1=2, max_j;
   {n_type max_f=-100; {for (int i=2; i<WL_max/2; i++) if (WL_weight[i]>max_f) {max_f=WL_weight[i]; max_j1=i;};}
   n_type Delta_WL=value("Additional_Delta_WL");
   max_j=max_j1; for (int i=max_j1; i<WL_max/2; i++) if (WL_weight[i]>max_f-Delta_WL) max_j=i;}


   n_type maxWL=WL_weight[max_j];
if (!(WL_flag==1 || WL_flag==2)) {max_j=N_max-1; maxWL=0;}
	{for (int i=0; i<max_j; i++) {WL_out<<int(100*(WL_weight[i]-maxWL)-1)/100.0<<"\n"; WL_weight[i]-=maxWL;};}
   {for (int i=max_j; i<N_max; i++) {WL_out<<0<<"\n";WL_weight[i]=0;}} WL_out<<flush;


    WL_flag=0;


   from_scratch();


      {     for (int z=0; z<n_zone; z++)
         for (int wn=0; wn<wn_max; wn++)
         for (int j=0; j<n_part; j++)
         {
         	GM00_st[z][wn][j]=0;
            for (int j2=0; j2<n_part; j2++) Gtotal_st[z][wn][j][j2]=0;
         ;}
   ;}

      {
   	for (n_type t=0; t<beta;t+=beta/(4*n_tau)-1e-5)
      {
         au<<double(t)<<"  ";
			for (int z=0; z<n_zone; z++)
			{
           	point p0(z,0,0);
         	for (int j=0; j<n_part; j++)
				{
            	point_ p0_(z,j,beta-t);
      			au<<-double(G0(p0,p0_))<<"  ";
            ;}
         ;}
         au<<"\n"<<flush;
      ;}
   ;}
;}
