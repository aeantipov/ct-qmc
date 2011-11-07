//=================================== weight =========================================


ComplexType G_ww(int , int);
//ofstream corr("corr.dat");

n_type weight()
{
	static ComplexType * g0[n_zone];
   static int f=0;
   if (f==0)
   {
      for (int z=0; z<n_zone; z++)
      {
      	g0[z]=new ComplexType[n_part];
	   	for (int j=0; j<n_part; j++)
      	{
         	Matrix l=Rotate(z), r(l); r.trans();
      		g0[z][j]=(l*g0w(z,0).x[j][j]*r).x[j][j];
      	;}
      ;}
      f=1;
   ;}

  
	n_type x=0;

   for (int z=0; z<n_zone; z++)
   for (int j=0; j<n_part; j++)
   for (int j2=0; j2<n_part; j2++)
   {
   	ComplexType y=GM_matrix[z][0][j][j2]; if (j==j2) y-=g0[z][j];//GM00[z][w][j]-g0[z][w][j];
   	x+=sqr(real(y))+sqr(imag(y));
   ;}


static int RR=int_value("calculate_Gamma4");
if (RR!=0)
   {
   	for (int z1=0; z1<n_zone; z1++) 
	for (int z2=0; z2<n_zone; z2++)
	for (int n1=0; n1<n_part; n1++)
	for (int n2=0; n2<n_part; n2++)
	{
   		ComplexType y=GM_matrix[z1][0][n1][n1]*GM_matrix[z2][0][n2][n2];
   		x+=sqr(real(y))+sqr(imag(y));
	;}
   ;}



static int NN=int_value("calculate_nn");
if (NN!=0 && RR==0)
{
	static int n1=int_value("nn_number1");
	static int n2=int_value("nn_number2");
	static int z1=int_value("nn_zone1");
	static int z2=int_value("nn_zone2");
	
	ComplexType y=GM_matrix[z1][0][n1][n1]*GM_matrix[z2][0][n2][n2];
	x+=sqr(real(y))+sqr(imag(y));	
;}


   return sqrt(x+n_part*n_zone)*exp(-WL_weight[N_WL]);
;}

//===================================== MC ===========================================
// trial step itself;
//explicitly modifies p, p_, u, a

void update_statistics()
{
if (isnan(curr_sgn/prev_weight)){cout<<"!";return;}
	
	GLOBAL_MC_COUNT++; hist[N_WL]++;
WL_weight[N_WL]+=WL_factor;

for (int z=0; z<n_zone; z++)
for (int wn=0; wn<wn_max; wn++)  
for (int j=0; j<n_part; j++)
{
   GM00_st[z][wn][j]+=(curr_sgn/prev_weight)*GM_interpolate(z,wn,j,j);
   for (int j2=0; j2<n_part; j2++)
   Gtotal_st[z][wn][j][j2]+=(curr_sgn/prev_weight)*GM_interpolate(z,wn,j,j2);
;}

weight_sum+=curr_sgn/prev_weight; weight_sum2+=1/prev_weight;
sgn_sum+=curr_sgn;
Nc_sum+=(Nc/2)*curr_sgn/prev_weight;
DispG0_sum+=sqr(real(GM_matrix[0][0][0][0])/prev_weight);
;}



void MC_single_step()
{
		if (rnd(2)==1)
   { //try add

DEBUG(" Add step ");
   	if ( (Nc>N_max-3) || (WL_flag==1 && 2*N_WL>WL_max-3) ) return;

      n_type uc;
      keep_state_all(0,0);
			W(p[Nc], p_[Nc], p[Nc+1], p_[Nc+1], uc, a[Nc], a[Nc+1]);
      	u[Nc]=uc; u[Nc+1]=1; wl[Nc]=WL_count; wl[Nc+1]=1;
//	#ifdef use_mpi
//	CTQMC.getStream() << " Before add_two " << flush;
//	#endif
      	uc*=-add_two(0);
//	#ifdef use_mpi
//	CTQMC.getStream() << " After add_two " << flush;
//	#endif
      	uc/=Nc/2;
         N_WL+=WL_count; WL_count=1;
      	n_type w=weight(); uc*=w/prev_weight; prev_weight=w;
         if (u[Nc-2]>0) curr_sgn*=-1;
		rnd();//for safety
      if (fabs(uc)>rnd()) {DEBUG("Accepted");add_two(1);ACCEPT_MC_COUNT++;} //accept
		else {DEBUG("Rejected"); add_two(-1);keep_state_all(1);}
   ;}
   else
   { //try remove
DEBUG(" Starting remove step ");


   	if (Nc<2) return;
      keep_state_all(0,0); //keep_state(0, p[K].z, p[K+1].z);
	DEBUG(endl);
	DEBUG("Here: Nc = " << Nc << " NumberOfPoints[]:  ")
	for (int z=0;z<n_zone;z++) DEBUG(NumberOfPoints[z] << " "); 
	DEBUG(endl);
      	int K=2*rnd(Nc/2); n_type uk=u[K];
			n_type uc=1/uk;  WL_count=wl[K];
      	uc*=Nc/2;
         N_WL-=WL_count; WL_count=1;
      	uc*=-remove_two(0,K); //if (uc<0) cout<<"!";
	DEBUG(endl);
	DEBUG("Here 2: Nc = " << Nc << " NumberOfPoints[]:  ")
	for (int z=0;z<n_zone;z++) DEBUG(NumberOfPoints[z] << " "); 
	DEBUG(endl);
      	n_type w=weight(); uc*=w/prev_weight;
         if (uk>0)   curr_sgn*=-1;prev_weight=w;
		rnd();//for safety
      if (fabs(uc)>rnd()) {DEBUG("Accepted"); 
	  
	  
	DEBUG(endl);
	DEBUG("00000000000: Nc = " << Nc << " NumberOfPoints[]:  ")
	for (int z=0;z<n_zone;z++) DEBUG(NumberOfPoints[z] << " "); 
	DEBUG(endl);
	  
	  
	  remove_two(1,K);ACCEPT_MC_COUNT++;}//accept
		else {DEBUG("Rejected"); remove_two(-1,K);keep_state_all(1);}
   ;}

;}


void MC_step(int n=1)
{
update_statistics();
if (Nc==0) curr_sgn=1;    //to escape accumulation of errors...
if (WL_flag!=0 && WL_weight[N_WL]>1e-10) {for (int i=0; i<N_max; i++) WL_weight[i]-=WL_factor;}

if (n==1) {MC_single_step(); 
return;}

	if (rnd(2)==1)
   { //try add
   	DEBUG("Add step");
   	if (Nc>N_max-1-2*n) return;
      if (WL_flag==1 && 2*N_WL>WL_max-1-2*n) return;
      n_type uc_=1, uc;
      keep_state_all(0,0);          //keep_state(0,p[Nc].z, p[Nc+1].z);
      for (int n_=0; n_<n; n_++)
      {
			W(p[Nc], p_[Nc], p[Nc+1], p_[Nc+1], uc, a[Nc], a[Nc+1]);
      	u[Nc]=uc; u[Nc+1]=1; wl[Nc]=WL_count; wl[Nc+1]=1;
      	uc*=-add()*add();
      	uc/=Nc/2;
         N_WL+=WL_count; WL_count=1;
      	n_type w=weight(); uc*=w/prev_weight; prev_weight=w;  uc_*=uc;
         if (u[Nc-2]>0) curr_sgn*=-1;
   	;}
		rnd();//for safety
      if (fabs(uc_)>rnd()) {ACCEPT_MC_COUNT++;} //accept
		else {keep_state_all(1);}
   ;}
   else
   { //try remove
   	DEBUG("Remove step");
   	if (Nc<2*n) return;
      keep_state_all(0,0); //keep_state(0, p[K].z, p[K+1].z);
      n_type uc_=1;
      for (int n_=0; n_<n; n_++)
      {
      	int K=2*rnd(Nc/2);
			n_type uc=1/u[K];  WL_count=wl[K];
      	uc*=Nc/2;
         N_WL-=WL_count; WL_count=1;
      	uc*=-remove(K+1); uc*=remove(K); //if (uc<0) cout<<"!";
      	n_type w=weight(); uc*=w/prev_weight;   uc_*=uc;
         if (u[Nc]>0)   curr_sgn*=-1;prev_weight=w;
      ;}
		rnd();//for safety
      if (fabs(uc_)>rnd()) {ACCEPT_MC_COUNT++;}//accept
		else {keep_state_all(1);}
   ;}

;}


n_type from_scratch();

void single_global_move(int ** z_new, int ** i_new, int num)
{
		static int gen_n=int_value("W_group_generators");
    	static int *** gen;
      static int f=0; if (f==0)
      {
      	gen=new int ** [2*gen_n];
      	{for (int n=0; n<2*gen_n; n++) {gen[n]=new int * [n_zone];for (int z=0; z<n_zone; z++) gen[n][z]=new int [n_part];};}
         for (int n=0; n<gen_n; n++)
         for (int z=0; z<n_zone; z++)
         for (int i=0; i<n_part; i++)
         {
         	gen[n][z][i]=int_value("W_group_generators",1+i+z*n_part+n_zone*n_part*n); //cout<<gen[n][z][i]<<"  "<<flush;
            int z1=gen[n][z][i]/n_part, i1=gen[n][z][i]%n_part;// cout<<z1<<i1<<" ";
            gen[n+gen_n][z1][i1]=z*n_part+i; //inverse element: used to deliver a detailed balance in Global Updates 
         ;}
      	f=1;
      ;}

      for (int z=0; z<n_zone; z++)
      for (int i=0; i<n_part; i++)
      {
      	int z1=z_new[z][i], i1=i_new[z][i];
         z_new[z][i]=gen[num][z1][i1]/n_part;
         i_new[z][i]=gen[num][z1][i1]%n_part;
      ;}

;}

void do_global_move(int ** z_new, int ** i_new)
{
   static int gen_n=int_value("W_group_generators");

    if (gen_n==0)
    {
      for (int z=0;z<n_zone; z++)
   	{
   		int f;
         do
   			{z_new[z][0]=rnd(n_zone); f=1; for (int i=0; i<z; i++) if (z_new[i][0]==z_new[z][0]) f=0;}
   		while(f==0);
         for (int i=0; i<n_part; i++) {z_new[z][i]=z_new[z][0]; i_new[z][i]=i;}
      ;}
    return;              //default: zones are identical; particles aren't
    ;}

    

    for (int z=0;z<n_zone; z++)
    for (int i=0; i<n_part; i++)
    	{z_new[z][i]=z; i_new[z][i]=i;}
    for (int j=0; j<rnd(n_part*n_zone);j++) single_global_move(z_new, i_new, rnd(2*gen_n));
;}

void GlobalMove()
{
	static int f=int_value("Use_Global_moves");
   n_type w0=from_scratch(); if (f==0) return;
   keep_state_all(0);
   int * z_new[n_zone], * i_new[n_zone];
   for (int z=0; z<n_zone;z++) {z_new[z]=new int [n_part]; i_new[z]=new int [n_part];}

   do_global_move(z_new, i_new);

   for(int j=0; j<Nc; j++)
   {
      int z=(p[j]).z, i=(p[j]).i; (p[j]).z=z_new[z][i]; (p[j]).i=i_new[z][i];
      z=(p_[j]).z; i=(p_[j]).i; (p_[j]).z=z_new[z][i]; (p_[j]).i=i_new[z][i];
   ;}

   n_type w1=from_scratch();
   if (w1>w0*rnd()) return;//accept
   else keep_state_all(1);
;}




