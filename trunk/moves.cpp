//================== Matsubara frequencies 2b used===================
bool use_wn(int wn)
{
	static bool * r=new bool[wn_max+1];
   static int f=0;
   if (f==0)
   {
      int k=int_value("sparse_Matsubaras");
		int k0=int_value("sparse_Matsubaras_start_from");
      for (int i=0; i<=wn_max; i++)
      {
      	r[i]=true;
         if (k>1 && i>0 && (k*(i/k)-k)>=k0 && (k*(i/k)+2*k)<wn_max && (i%k)!=0) r[i]=false;
         //if (r[i]) cout<<i<<" ";
      ;}
   	f=1;
   ;}
   return r[wn];
;}

ComplexType GM_interpolate(int z, int wn, int q1, int q2)
{
	if (use_wn(wn)) return GM_matrix[z][wn][q1][q2];
   static int k=int_value("sparse_Matsubaras");
   int w1=k*(wn/k); n_type x=(n_type(wn%k))/k;                      //  cout<<wn<<">"<<w1<<" "<<flush;
   ComplexType y0=GM_matrix[z][w1-k][q1][q2], y1=GM_matrix[z][w1][q1][q2], y2=GM_matrix[z][w1+k][q1][q2],y3=GM_matrix[z][w1+2*k][q1][q2];
   return y1-(2.*y0+3.*y1-6.*y2+y3)*x/6.+(y0-2.*y1+y2)*x*x/2.-(y0-3.*y1+3.*y2-y3)*x*x*x/6.;
;}



//==================== simple aux. functions =========================
void next(int z, int & prev_) //gives the next point of a given zone
{
	do
   {
   	prev_++;
      if ((p[prev_]).z==z) return;
   ;}
   while(prev_<Nc-1);
   prev_++; //no point => return Nc;
;}

void prev(int z, int & next_) //gives the previous point of a given zone
{
	do
   {
   	next_--;
      if ((p[next_]).z==z) return;
   ;}
   while(next_>0);
   next_--; //no point => return -1
;}

//========================= subroutines for Green function calcs =====================

ComplexType MG(int z, int wn, int k, int q, int* kk)     //k: index in M, q: index in observable-space
{
   ComplexType s=0;
	ComplexType ** x=(*Grot(z,wn,1)).x; 
	for (int j=0; j<NumberOfPoints[z]; j++)
       	s+=M[z][k][j]*(x[p[kk[j]].i][q]*ewt_[kk[j]][wn]);
	return s;
;}



ComplexType GM(int z, int wn, int k, int q, int* kk)
{
   ComplexType s=0;
	ComplexType **x=(*Grot(z,wn,-1)).x;
	for (int j=0; j<NumberOfPoints[z]; j++)
			s+=M[z][j][k]*(x[q][p_[kk[j]].i]*ewt[kk[j]][wn]);
	return s;
;}




//=============================  add-remove single point =============================
//low-level operations
//only these functions "constructively" modify global variables M, GM00 etc.
//(formally, these variables are modified also by keep_state() and from_scratch())

void add_point()  //just adds a point; M 2b recalculated separately
{
			n_type t=p[Nc].t, t_=p_[Nc].t;
         if (sqr(t-t_)>1e-13)
         for (int wn=0; wn<wn_max; wn++)
         {
        		n_type w=Pi*(2*wn+1)/beta;
            n_type a=w*t, a_=w*t_;
        		ewt[Nc][wn]=ComplexType(cos(a), sin(a));//exp(I*w*p[Nc].t);
            ewt_[Nc][wn]=ComplexType(cos(a_), -sin(a_));//exp(-I*w*p_[Nc].t);
         ;}
         else
         for (int wn=0; wn<wn_max; wn++)
         {
          	n_type w=Pi*(2*wn+1)/beta;
            n_type a=w*t, ca=cos(a), sa=sin(a);                       //try more optimization with sin(x+y) formula!
        		ewt[Nc][wn]=ComplexType(ca, sa);//exp(I*w*p[Nc].t);
            ewt_[Nc][wn]=ComplexType(ca,-sa);//exp(-I*w*p_[Nc].t);
         ;}


		  int z=p[Nc].z;Nc++;

        if (z==n_zone) return;
        NumberOfPoints[z]++;
        for (int i=0; i<NumberOfPoints[z]; i++)
        {
        		M[z][NumberOfPoints[z]-1][i]=0; M[z][i][NumberOfPoints[z]-1]=0;
        ;}

;}

void remove_point(int K) //same thing
{

		  int z=p[K].z;
        int K0=-1, k=-1;
        do {next(z,K0);k++;} while (K0!=K);

        Nc--;
        point pK=p[K]; point_ pK_=p_[K];n_type uK=u[K], aK=a[K]; int wlK=wl[K];
        ComplexType *ewtK=ewt[K], *ewtK_=ewt_[K];
        for (int i=K; i<Nc; i++)
        		{p[i]=p[i+1];p_[i]=p_[i+1];u[i]=u[i+1];a[i]=a[i+1];wl[i]=wl[i+1];ewt[i]=ewt[i+1];ewt_[i]=ewt_[i+1];}
        p[Nc]=pK; p_[Nc]=pK_;u[Nc]=uK;a[Nc]=aK; wl[Nc]=wlK; ewt[Nc]=ewtK; ewt_[Nc]=ewtK_;

        if (z==n_zone) return;
        NumberOfPoints[z]--;
        {for (int i=k; i<NumberOfPoints[z]; i++)
         for (int l=0; l<k; l++)
         {
                M[z][i][l]=M[z][i+1][l];
                M[z][l][i]=M[z][l][i+1];
         ;}
        ;}

        {for (int i=k; i<NumberOfPoints[z]; i++)
         for (int l=k; l<NumberOfPoints[z]; l++)
         {
                M[z][i][l]=M[z][i+1][l+1];
         ;}
        ;}
;}

n_type modif_remove_1(int K)
{
        int k=-1, K_=-1, z=(p[K]).z;  
        do {next(z,K_); k++;} while (K_!=K); //k is determined

   	  n_type y=M[z][k][k];                             //y is determined

		  static ComplexType * gm=new ComplexType[n_part], * mg=new ComplexType[n_part];
        n_type d=beta*M[z][k][k];
        {
			static int * kk_array= new int [N_max];
			{int kk=-1; for (int j=0; j<NumberOfPoints[z]; j++) {next(z,kk); kk_array[j]=kk;};}
	        for (int q=0; q<n_part; q++) {gm[q]=GM(z,0,k,q, kk_array); mg[q]=MG(z,0,k,q, kk_array);}
        ;}

        for (int q=0; q<n_part; q++)
        for (int q2=0; q2<n_part; q2++)
              		GM_matrix[z][0][q][q2]+=gm[q]*mg[q2]/d;  //G[wn=0] is determined
        return y;
;}






void modif_remove_2(int K)
{
        int k=-1, K_=-1, z=(p[K]).z;  
        do {next(z,K_); k++;} while (K_!=K); //k is determined

   	  n_type y=M[z][k][k];                             //y is determined

        static ComplexType ** gm, ** mg;//[WN_max][n_part];
        static int gf=0; if (gf==0)
        {
        		gm=new ComplexType *[wn_max];mg=new ComplexType *[wn_max];
            for (int i=0; i<wn_max; i++)
            {
            	gm[i]=new ComplexType [n_part];
               mg[i]=new ComplexType [n_part];
            ;}
        gf=1;
        ;}

        {
			static int * kk_array= new int [N_max];
			{int kk=-1; for (int j=0; j<NumberOfPoints[z]; j++) {next(z,kk); kk_array[j]=kk;};}
//t_counter-=clock();
           for (int wn=1; wn<wn_max; wn++) if(use_wn(wn))
	        for (int q=0; q<n_part; q++) {gm[wn][q]=GM(z,wn,k,q, kk_array); mg[wn][q]=MG(z,wn,k,q, kk_array);}
//t_counter+=clock();
        ;}

		  //fast update of G
        n_type d=beta*M[z][k][k];
        {
        for (int wn=1; wn<wn_max; wn++) if(use_wn(wn))
        for (int q=0; q<n_part; q++)
        for (int q2=0; q2<n_part; q2++)
           		GM_matrix[z][wn][q][q2]+=gm[wn][q]*mg[wn][q2]/d;

        ;}              //rest of Matsubara's is determinted



  	     n_type Ml[N_max], Mr[N_max];
        for (int i=0; i<NumberOfPoints[z]; i++)
        {
            Ml[i]=-M[z][i][k];
            Mr[i]=M[z][k][i]/y;
        ;}

        {for (int i=0; i<NumberOfPoints[z]; i++)
         {n_type ml=Ml[i]; int ncz=NumberOfPoints[z];
         	#pragma vector always
         	for (int j=0; j<ncz; j++)  M[z][i][j]+=ml*Mr[j];
         ;}
        ;}

;}


n_type modif_remove(int K)
{
	n_type r=modif_remove_1(K);  modif_remove_2(K); return r;
;}


n_type modif_add_1(int f2)
{
        int z=(p[Nc-1]).z;

   	  n_type y=G0(p[Nc-1], p_[Nc-1]);
  	     n_type Ml[N_max], Mr[N_max];
        for (int i=0; i<NumberOfPoints[z]; i++)  {Ml[i]=0; Mr[i]=0;}

        n_type r[N_max], l[N_max];
        int k2=-1, K_=-1; do {next(z,K_); k2++; l[k2]=G0(p[K_], p_[Nc-1]); r[k2]=G0(p[Nc-1], p_[K_]);} while (K_<Nc);

        if (f2==1) {r[NumberOfPoints[z]-1]-=a[Nc-1];l[NumberOfPoints[z]-1]-=a[Nc-1];y-=a[Nc-1];}


        {for (int i=0; i<NumberOfPoints[z]-1; i++)
        #pragma vector always
         for (int j=0; j<NumberOfPoints[z]-1; j++)
         {
         	Ml[i]+=M[z][i][j]*l[j];
            Mr[i]+=r[j]*M[z][j][i];
         ;}
         ;}

     		Ml[NumberOfPoints[z]-1]=-1;Mr[NumberOfPoints[z]-1]=-1;
        {for (int i=0; i<NumberOfPoints[z]-1; i++) y-=l[i]*Mr[i];}         //y is determined
        {for (int i=0; i<NumberOfPoints[z]; i++) Mr[i]/=y;}

        {for (int i=0; i<NumberOfPoints[z]; i++)
         {n_type ml=Ml[i]; int ncz=NumberOfPoints[z];
         	#pragma vector always
         	for (int j=0; j<ncz; j++) M[z][i][j]+=ml*Mr[j];
         ;}
        ;}                                                        //M is changed


        static ComplexType * gm=new ComplexType[n_part], * mg=new ComplexType[n_part];

        static int * kk_array= new int [N_max];
		  {int kk=-1; for (int j=0; j<NumberOfPoints[z]; j++) {next(z,kk); kk_array[j]=kk;};}
        for (int q=0; q<n_part; q++) {gm[q]=GM(z,0,NumberOfPoints[z]-1,q, kk_array); mg[q]=MG(z,0,NumberOfPoints[z]-1,q, kk_array);}



        n_type d=beta*M[z][NumberOfPoints[z]-1][NumberOfPoints[z]-1];
		  //fast update of G
        {
        for (int q=0; q<n_part; q++)
        for (int q2=0; q2<n_part; q2++)
              		GM_matrix[z][0][q][q2]-=gm[q]*mg[q2]/d;        //G[w=0]
        ;}


               return y;
;}




void modif_add_2(int z)
{

        static ComplexType ** gm, ** mg;//[WN_max][n_part];
        static int gf=0; if (gf==0)
        {
        		gm=new ComplexType *[wn_max];mg=new ComplexType *[wn_max];
            for (int i=0; i<wn_max; i++)
            {
            	gm[i]=new ComplexType [n_part];
               mg[i]=new ComplexType [n_part];
            ;}
        gf=1;
        ;}

        {static int * kk_array= new int [N_max];
			{int kk=-1; for (int j=0; j<NumberOfPoints[z]; j++) {next(z,kk); kk_array[j]=kk;};}

//t_counter-=clock();
           for (int wn=1; wn<wn_max; wn++)  if(use_wn(wn))
	        for (int q=0; q<n_part; q++) {gm[wn][q]=GM(z,wn,NumberOfPoints[z]-1,q, kk_array); mg[wn][q]=MG(z,wn,NumberOfPoints[z]-1,q, kk_array);}
//t_counter+=clock();
        ;}


        n_type d=beta*M[z][NumberOfPoints[z]-1][NumberOfPoints[z]-1];
		  //fast update of G
        {
        for (int wn=1; wn<wn_max; wn++) if(use_wn(wn))
        for (int q=0; q<n_part; q++)
        for (int q2=0; q2<n_part; q2++)
              		GM_matrix[z][wn][q][q2]-=gm[wn][q]*mg[wn][q2]/d;       //rest of Matsubara's
        ;}


;}


n_type modif_add(int f2)
{
	n_type r=modif_add_1(f2);
   int z=p[Nc-1].z;
   modif_add_2(z);
   return r;
;}

//============================== add-remove ==========================================
//higher-level operations

n_type add(int f2=1)
{
   add_point();
   n_type f=modif_add(f2);
   if (f<0) curr_sgn*=-1;

   return f;
;}

n_type add_two(int f)  //f==1=>finalize
{
static bool same_zone;   static int z1,z2;
if (f==0)
{
	same_zone=(p[Nc].z==p[Nc+1].z);      z1=p[Nc].z; z2=p[Nc+1].z;
   if (same_zone) return add()*add();
   else
   {
   	add_point();n_type f=modif_add_1(1);
      add_point();      f*=modif_add_1(1); 
      if (f<0) curr_sgn*=-1; return f;
   ;}
;}
else
{
	if (!same_zone) {if (f>0) {modif_add_2(z1);modif_add_2(z2);};}
	return 1;
;}
;}

n_type remove(int K)
{

	n_type f=modif_remove(K);
	remove_point(K);
   if (f<0) curr_sgn*=-1;

   return f;
;}

n_type remove_two(int f,int K)    //f==1=>finalize
{
static bool same_zone;
if (f==0)
{
	 same_zone=(p[K].z==p[K+1].z);
      if (same_zone) {n_type r=remove(K+1); return r*remove(K);}
    else{n_type f=modif_remove_1(K+1); f*=modif_remove_1(K); if (f<0) curr_sgn*=-1; return f;}
}
else
{
	if (!same_zone) {if (f>0) {modif_remove_2(K+1); modif_remove_2(K);}     remove_point(K+1); remove_point(K);}
	return 1;
;}

;}



