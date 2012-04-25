//====================================================================================


void keep_state_all(int f, int keep_ewt=1) //f==0 => keep, f==1 => recover
{

   static int curr_sgn0, N_WL0;
   static n_type prev_weight0;
	static int Nc0, NumberOfPoints0[n_zone];
	static point p0[N_max]; static point_ p0_[N_max];
   static n_type u0[N_max], a0[N_max];
   static int wl0[N_max];
   static ComplexType * ewt0[N_max], *ewt0_[N_max];
   static ComplexType * ewt00[N_max], *ewt00_[N_max];  //no memory 2b allocated for these arrays...
   static int ewt_kept;
   static n_type **M0[n_zone];
   static ComplexType *** GM_matrix0[n_zone];




   static int m=0;
   if (m==0)
   {
	   for (int z=0; z<n_zone; z++)
	   {
   		NumberOfPoints0[z]=0;
   		M0[z]=new n_type * [N_max];
      	for (int j=0; j<N_max/3; j++) //Why?
      		M0[z][j]=new n_type [N_max];

   	GM_matrix0[z]=new ComplexType **[wn_max];
   	for (int w=0; w<wn_max; w++)
   	{
   	  GM_matrix0[z][w]=new ComplexType *[n_part];
         for (int n1=0; n1<n_part; n1++)
         	GM_matrix0[z][w][n1]=new ComplexType [n_part];
   	  ;}

   	;}
      for (int i=0; i<N_max; i++){ewt0[i]=new ComplexType [wn_max]; ewt0_[i]=new ComplexType [wn_max];}
      m=1;
   ;}


   if (f==0)
   {

   	Nc0=Nc; N_WL0=N_WL; {for (int i=0; i<n_zone; i++) NumberOfPoints0[i]=NumberOfPoints[i];} curr_sgn0=curr_sgn; prev_weight0=prev_weight;

      {
      for (int i=0; i<Nc; i++) {p0[i]=p[i]; p0_[i]=p_[i]; a0[i]=a[i]; u0[i]=u[i]; wl0[i]=wl[i];}
       														 // memcpy(ewt0[i],ewt[i],wn_max*sizeof(ComplexType));   memcpy(ewt0_[i],ewt_[i],wn_max*sizeof(ComplexType));
      ewt_kept=0;
      if (keep_ewt==1)
         {
         	for (int i=0; i<Nc; i++)
         	for (int wn=0; wn<wn_max; wn++)
         		{ewt0[i][wn]=ewt[i][wn];ewt0_[i][wn]=ewt_[i][wn];}
            ewt_kept=1;
         ;}
         else
         for (int i=0; i<Nc; i++)
         {
         	ewt00[i]=ewt[i];ewt00_[i]=ewt_[i];
         ;}

      ;}

      {
         for(int z=0; z<n_zone; z++)
      	for(int i=0; i<NumberOfPoints[z]; i++)
         	memcpy(M0[z][i],M[z][i],NumberOfPoints[z]*sizeof(n_type));
      	//for(int j=0; j<NumberOfPoints[z]; j++) M0[z][i][j]=M[z][i][j];
      ;}

	   {
    		for(int z=0; z<n_zone; z++)
    		for (int nn=0; nn<n_part; nn++)
    		for (int nn2=0; nn2<n_part; nn2++)
     		for (int wn=0; wn<wn_max; wn++)
             	GM_matrix0[z][wn][nn][nn2]=GM_matrix[z][wn][nn][nn2];
      ;}

   ;}

   if (f==1)
   {
   	Nc=Nc0; N_WL=N_WL0; {for (int i=0; i<n_zone; i++) NumberOfPoints[i]=NumberOfPoints0[i];} curr_sgn=curr_sgn0; prev_weight=prev_weight0;
      {for (int i=0; i<Nc; i++) {p[i]=p0[i]; p_[i]=p0_[i]; a[i]=a0[i]; u[i]=u0[i]; wl[i]=wl0[i];}

      if (ewt_kept==1)
         {
         	for (int i=0; i<Nc; i++)
         	for (int wn=0; wn<wn_max; wn++)
         	{ewt[i][wn]=ewt0[i][wn];ewt_[i][wn]=ewt0_[i][wn];}
         ;}
         else
         for (int i=0; i<Nc; i++)
         {
         	ewt[i]=ewt00[i];ewt_[i]=ewt00_[i];
         ;}
      ;}
      {
      	for(int z=0; z<n_zone; z++)
      	for(int i=0; i<NumberOfPoints[z]; i++)
         	memcpy(M[z][i],M0[z][i],NumberOfPoints[z]*sizeof(n_type));
//      	for(int j=0; j<NumberOfPoints[z]; j++) M[z][i][j]=M0[z][i][j];
      ;}

      {
      	for(int z=0; z<n_zone; z++)
    		for (int nn=0; nn<n_part; nn++)
        	for (int nn2=0; nn2<n_part; nn2++)
    		for (int wn=0; wn<wn_max; wn++)
                  	GM_matrix[z][wn][nn][nn2]=GM_matrix0[z][wn][nn][nn2];

    	;}


   ;}


;}

