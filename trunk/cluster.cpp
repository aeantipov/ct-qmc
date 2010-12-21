//============== random walk in the cluster =================
// not implemented yet


int cluster_size=int_value("cluster_size");
int cluster_curr_size;
int cluster_nz [n_zone];

n_type ** cluster_M[n_zone], ** cluster_M0[n_zone];
ComplexType ** cluster_gM[n_zone], ** cluster_Mg[n_zone];
ComplexType ** cluster_gM0[n_zone], ** cluster_M0g[n_zone];
int * cluster_state; int * cluster_state_z[n_zone];
int * cluster_z1, * cluster_z2, * cluster_k1, * cluster_k2;
int * cluster_KM1, *cluster_KM2, *cluster_KN;
n_type * cluster_u; int * cluster_wl;
ComplexType ** cluster_Gmatr[n_zone];



void cluster_print_M(int z)  //debug
{
	cout<<"\nclusterM z="<<z<<"\n";
	for (int i=0; i<cluster_nz[z]; i++)
   {
	   for (int j=0; j<cluster_nz[z]; j++)
      	cout<<cluster_M[z][i][j]<<"  "<<flush;
      cout<<"\n";
   ;}
;}

void cluster_print_M0(int z)  //debug
{
	cout<<"\nclusterM0 z="<<z<<"\n";
	for (int i=0; i<cluster_nz[z]; i++)
   {
	   for (int j=0; j<cluster_nz[z]; j++)
      	cout<<cluster_M0[z][i][j]<<"  "<<flush;
      cout<<"\n";
   ;}
;}


void cluster_ini() //memory alloc.
{
	{static int f=0; if (f==1) return; f=1;}
	for (int z=0; z<n_zone; z++)
   {
   	cluster_state_z[z]=new int [2*cluster_size];
   	cluster_M[z] =new n_type * [2*cluster_size];
      cluster_M0[z]=new n_type * [2*cluster_size];
      for (int k=0; k<2*cluster_size; k++)
      {
	      cluster_M[z][k]=new n_type [2*cluster_size];
   	   cluster_M0[z][k]=new n_type [2*cluster_size];
      ;}
      cluster_Mg[z]=new ComplexType * [2*cluster_size];
      cluster_M0g[z]=new ComplexType * [2*cluster_size];
      for (int k1=0; k1<2*cluster_size; k1++)
      {
      	cluster_Mg[z][k1]=new ComplexType [n_part];
         cluster_M0g[z][k1]=new ComplexType [n_part];
      ;}
      cluster_gM[z]=new ComplexType * [n_part];
      cluster_gM0[z]=new ComplexType * [n_part];
      for (int k2=0; k2<n_part; k2++)
      {
      	cluster_gM[z][k2]=new ComplexType [2*cluster_size];
         cluster_gM0[z][k2]=new ComplexType [2*cluster_size];
      ;}
   	cluster_Gmatr[z]=new ComplexType * [n_part];
 	   for (int q=0; q<n_part; q++)
   	cluster_Gmatr[z][q]=new ComplexType[n_part];

   ;}
   cluster_state=new int [cluster_size];
   cluster_z1=new int [cluster_size];
   cluster_z2=new int [cluster_size];
   cluster_k1=new int [cluster_size];
   cluster_k2=new int [cluster_size];
   cluster_u= new n_type [cluster_size];
   cluster_wl=new int [cluster_size];
   cluster_KM1=new int [cluster_size];
   cluster_KM2=new int [cluster_size];
   cluster_KN= new int [cluster_size];

;}

n_type remove_single(int z, int k)
{
	cluster_state_z[z][k]=0;
	n_type y=cluster_M[z][k][k];
	static n_type * Ml=new n_type [2*cluster_size];
	static n_type * Mr=new n_type [2*cluster_size];
   for (int i=0; i<cluster_nz[z]; i++)
   {
   	Ml[i]=-cluster_M[z][i][k];
      Mr[i]=cluster_M[z][k][i]/y;
   ;}

//change of G
	{
   for (int q=0;  q<n_part;  q++)
   for (int q2=0; q2<n_part; q2++)
   	GM_matrix[z][0][q][q2]+=cluster_gM[z][q][k]*cluster_Mg[z][k][q2]/(beta*y);
   ;}

//change of GM, MG
   {
   	for (int q=0; q<n_part; q++)
      for (int i=0; i<cluster_nz[z]; i++)
      if (cluster_state_z[z][i]!=0)
      {
      	cluster_gM[z][q][i]-=cluster_gM[z][q][k]*cluster_M[z][k][i]/y;
         cluster_Mg[z][i][q]-=cluster_Mg[z][k][q]*cluster_M[z][i][k]/y;
      ;}
   ;}

//change of M
   {
   for (int i=0; i<cluster_nz[z]; i++)
   for (int j=0; j<cluster_nz[z]; j++)
   	if (cluster_state_z[z][i]!=0 && cluster_state_z[z][j]!=0)
      cluster_M[z][i][j]+=Ml[i]*Mr[j];
   ;}


	return y;
;}

n_type cluster_weight()
{

	for (int z=0; z<n_zone; z++)
   for (int i=0; i<cluster_nz[z]; i++)
   {
   	for (int k=0; k<cluster_nz[z]; k++) cluster_M[z][i][k]=cluster_M0[z][i][k];
      for (int q=0; q<n_part; q++)
      {
      	cluster_Mg[z][i][q]=cluster_M0g[z][i][q];
         cluster_gM[z][q][i]=cluster_gM0[z][q][i];
      ;}
      cluster_state_z[z][i]=1;
   ;}

	n_type w=1; int wl=0;
   for (int i=0; i<cluster_size; i++)
 	if (cluster_state[i]==0)
   {
      	w*=remove_single(cluster_z1[i], cluster_k1[i]);
         w*=remove_single(cluster_z2[i], cluster_k2[i]);
         w/=cluster_u[i];
         wl+=cluster_wl[i];
   ;}
   N_WL-=wl; w*=weight(); N_WL+=wl;

   {
   	for (int z=0; z<n_zone; z++)
   	for (int q1=0;q1<n_part; q1++)
      for (int q2=0;q2<n_part; q2++)
   	GM_matrix[z][0][q1][q2]=cluster_Gmatr[z][q1][q2];
   ;}
	return w;
;}



void cluster_define()
{
   {
   	for (int z=0; z<n_zone; z++)
   	for (int q1=0;q1<n_part; q1++)
      for (int q2=0;q2<n_part; q2++)
   	cluster_Gmatr[z][q1][q2]=GM_matrix[z][0][q1][q2];
   ;}

   {for (int z=0; z<n_zone; z++) cluster_nz[z]=0;}
   for (int i=0; i<cluster_size; i++)
   {
   	cluster_u[i]=u[2*cluster_KN[i]];
      cluster_wl[i]=wl[2*cluster_KN[i]];                            //??????  wl[Nc]

      int z=p[2*cluster_KN[i]].z;
      cluster_z1[i]=z;
      cluster_k1[i]=cluster_nz[z];

      int K=-1;{int kk=-1; for (int j=0; kk!=2*cluster_KN[i] && j<NumberOfPoints[z]; j++) {next(z,kk);K++;}}
      cluster_KM1[i]=K;
      cluster_nz[z]++;


      z=p[2*cluster_KN[i]+1].z;
      cluster_z2[i]=z;
      cluster_k2[i]=cluster_nz[z];

      K=-1;{int kk=-1; for (int j=0; kk!=2*cluster_KN[i]+1 && j<NumberOfPoints[z]; j++) {next(z,kk);K++;}}
      cluster_KM2[i]=K;   
      cluster_nz[z]++;
   ;}



{
for (int z=0; z<n_zone; z++)
{
	static int * kk_array=new int [N_max];
   {int kk=-1; for (int j=0; j<Nc; j++) {next(z,kk); kk_array[j]=kk;};}

   for (int i=0; i<cluster_size; i++)
   {
   	for (int j=0; j<cluster_size; j++)
   	{
   		if (z==cluster_z1[i] && z==cluster_z1[j])
         	cluster_M0[z][cluster_k1[i]][cluster_k1[j]]=M[z][cluster_KM1[i]][cluster_KM1[j]];
         if (z==cluster_z2[i] && z==cluster_z2[j])
      		cluster_M0[z][cluster_k2[i]][cluster_k2[j]]=M[z][cluster_KM2[i]][cluster_KM2[j]];
   		if (z==cluster_z1[i] && z==cluster_z2[j])
         	cluster_M0[z][cluster_k1[i]][cluster_k2[j]]=M[z][cluster_KM1[i]][cluster_KM2[j]];
   		if (z==cluster_z2[i] && z==cluster_z1[j])
         	cluster_M0[z][cluster_k2[i]][cluster_k1[j]]=M[z][cluster_KM2[i]][cluster_KM1[j]];
   	;}
      for (int q=0; q<n_part; q++)
      {
      	if (z==cluster_z1[i])
         {
      		cluster_gM0[z][q][cluster_k1[i]]=GM(z,0,cluster_KM1[i],q,kk_array);
            cluster_M0g[z][cluster_k1[i]][q]=MG(z,0,cluster_KM1[i],q,kk_array);
         ;}
         if (z==cluster_z2[i])
         {
      		cluster_gM0[z][q][cluster_k2[i]]=GM(z,0,cluster_KM2[i],q,kk_array);
            cluster_M0g[z][cluster_k2[i]][q]=MG(z,0,cluster_KM2[i],q,kk_array);
         ;}
      ;}
   ;}

;}
;}
;}

n_type factorial (int i)
{
	n_type s=1;
	for (int j=1; j<=i; j++) s*=j;
   return s;
;}



n_type cluster_factor()
{
   n_type w=1; int n=0;
   for (int i=0; i<cluster_size; i++)
   	if (cluster_state[i]==0) {w*=Nc/2-n; n++;}
   w*=factorial(n)*factorial(cluster_size-n);
   return w;
;}


int cluster_walks()
{

	cluster_define();
//   {for (int k=0; k<cluster_size; k++)  if (rnd()<0.1) cluster_state[k]=1; else cluster_state[k]=0; return 1;}



   static int num=int_value("number_of_trials_in_cluster");
	static n_type *  ww=new n_type [num+1];
   static int ** states;
   {static int f=0; if (f==0)
   {states= new int * [cluster_size]; for (int i=0; i<cluster_size; i++) states[i]=new int [num+1];f=1;};}



   n_type norm=0;
   for (int i=0; i<num+1; i++)
   {
   	ww[i]=cluster_factor();
      n_type clw=cluster_weight();
      ww[i]*=fabs(clw);

      norm+=ww[i];
      {for (int k=0; k<cluster_size; k++) {states[k][i]=cluster_state[k];}}


		for (;i<num;)
		{

      	{for (int k=0; k<cluster_size; k++)  cluster_state[k]=rnd(2);}

         int f=0;
         for (int l=0; l<i+1; l++)
         {
         	int f1=0;
         	for (int k=0; k<cluster_size; k++) if (cluster_state[k]!=states[k][l]) f1=1;
            if (f1==0) f=1;
         ;}

         if (f==0) break;
		;}

   ;}
   n_type s=0, r=rnd();
   int j=-1;
   do
   {
   	j++;
   	s+=ww[j]/norm;
   ;}
   while (j<num && s<r);


   for (int k=0; k<cluster_size; k++) cluster_state[k]=states[k][j];



	if (j==0) return 0; else return 1;
;}


void cluster_step()
{

   update_statistics();

	int n_new=rnd(cluster_size+1);
   if (Nc/2+n_new-cluster_size<0) return; 			 //incorrect cluster - reject

   {for (int i=0; i<n_new; i++)                     //points 2b added
   {
   		cluster_KN[i]=Nc/2;
         cluster_state[i]=0;
			W(p[Nc], p_[Nc], p[Nc+1], p_[Nc+1], u[Nc], a[Nc], a[Nc+1]);
         if (u[Nc]>0) curr_sgn*=-1;
      	u[Nc+1]=1; N_WL+=WL_count; wl[Nc]=WL_count; wl[Nc+1]=1;
      	add();add(); WL_count=1;
   ;}
   ;}

   {for (int i=0; i<cluster_size-n_new; i++)        	//points 2b deleted
   	{
      	int c;
         do
         {  c=0;
         	int l=0; if (Nc/2-n_new>1) l=rnd(Nc/2-n_new);
   			cluster_KN[i+n_new]=l;  cluster_state[i+n_new]=1;
         	for (int j=0; j<i; j++)
            	if (cluster_KN[i+n_new]==cluster_KN[j+n_new]) c=1;
         ;}
         while (c==1);
      ;}
   ;}


int ac=cluster_walks();


   ACCEPT_MC_COUNT+=ac;
   //ordering

   {
	   for (int i=0; i<cluster_size-1; i++)
      for (int j=i+1; j<cluster_size; j++)
      {
      	if (cluster_KN[i]<cluster_KN[j])
         {
            int k=cluster_KN[i], s=cluster_state[i];
            cluster_KN[i]=cluster_KN[j]; cluster_state[i]=cluster_state[j];
            cluster_KN[j]=k; cluster_state[j]=s;
         ;}
      ;}
   ;}

   {
	   for (int i=0; i<cluster_size; i++)
      if (cluster_state[i]==0) {int k=2*cluster_KN[i]; N_WL-=wl[k]; if (u[k]>0) curr_sgn*=-1; remove(k+1);remove(k);}
      prev_weight=weight();
   ;}
;}





