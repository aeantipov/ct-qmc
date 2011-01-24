//=================== effective g0 are calculated here ===========================
#ifdef use_mpi
extern CTQMC_WORLD CTQMC;
#endif

Matrix & Rotate (int); //the body is in input.cpp




Matrix & alpha_W(int z)
{
	int NNN=int(alphaW_number); if (NNN<10) NNN=100000000;
	static Matrix * aW=new Matrix [n_zone];
	static n_type n_aver=value("expected_occupancy");
	static int f=0;
   if (f==0)
   {
   	{for (int i0=0; i0<n_zone; i0++)
    		{aW[i0].new_memory(); aW[i0]=0;}
      ;}
      int t=time(NULL);
      #ifdef use_mpi
      CTQMC.getStream()<<"\nalpha_W data accumulation ..."<<flush;
      #else
      cout<<"\nalpha_W data accumulation ..."<<flush;
      #endif
      point r1, r2; point_ r1_, r2_; n_type u_disp=0, u, a1, a2;
      for (int i=0; i<NNN; i++)
      {
      	W(r1, r1_, r2, r2_, u, a1, a2);
         (aW[r1.z]).x[r1.i][r1_.i]+=ComplexType(a2*u/(beta*NNN));
         if (r1.i==r1_.i) (aW[r1.z]).x[r1.i][r1_.i]-=n_aver*u/(beta*NNN);
         (aW[r2.z]).x[r2.i][r2_.i]+=ComplexType(a1*u/(beta*NNN));
         if (r2.i==r2_.i) (aW[r2.z]).x[r2.i][r2_.i]-=n_aver*u/(beta*NNN);
         u_disp+=sqr(u/(beta*NNN));

         if (i==NNN/10)
         {
         	n_type s=0;
      		{for (int i0=0; i0<n_zone; i0++) s+=norm2(aW[i0]);}
         	if (s<u_disp*30)
            {
      		#ifdef use_mpi
            	CTQMC.getStream()<<"\nforced no correction to G0\n";
		#else
            	cout<<"\nforced no correction to G0\n";
		#endif
            	{for (int i0=0; i0<n_zone; i0++) aW[i0]=0;}
               break;
            ;}
         ;}
      ;}
      #ifdef use_mpi
      CTQMC.getStream()<<" done in "<<time(NULL)-t<<" sec.\n"<<flush;
      #else
      cout<<" done in "<<time(NULL)-t<<" sec.\n"<<flush;
      #endif
      f=1;

//      if (norm2(aW[0]-aW[1])<1e-4) {Matrix a=0.5*(aW[0]+aW[1]); aW[0]=a; aW[1]=a; cout<<"alpha_W symmetrized\n"<<flush;}
   ;}
   return aW[z];
;}


Matrix & G0w_(int z, int wn)
{
	static Matrix g;
	g=g0w(z,wn);
   	Inverse(g);
			g=g-alpha_W(z);
   	Inverse(g);
	return g;
;}



Matrix * G0w(int z, int wn)
{
	static Matrix ** G0wc [n_zone];
   static int k=0;
   if (k==0)
   {		
	{for (int zz=0; zz<n_zone; zz++) G0wc[zz]=new Matrix * [WN_max+1];}
      for (int zz=0; zz<n_zone; zz++)
  	for (int ww=0; ww<wn_max; ww++)
  	{
      	G0wc[zz][ww]=new Matrix(0); (*G0wc[zz][ww]).new_memory();
         *G0wc[zz][ww]=G0w_(zz,ww);
      ;}
   k=1;
   ;}
	return G0wc[z][wn];
;}

Matrix * Grot(int z, int wn, int f) //-1: RG 1: GR^+ 0: RGR^+
{
	static int rb=int_value("rotate_basis");   if (rb==0) return G0w(z,wn);

	static Matrix ** G0l  [n_zone];//[WN_max+1];
	static Matrix ** G0r  [n_zone];//[WN_max+1];
	static Matrix ** G0lr [n_zone];//[WN_max+1];
   static int k=0;
   if (k==0)
   {
      {for (int zz=0; zz<n_zone; zz++) {G0l[zz]=new Matrix * [WN_max+1];G0r[zz]=new Matrix * [WN_max+1];G0lr[zz]=new Matrix * [WN_max+1];};}
      for (int zz=0; zz<n_zone; zz++)
  	for (int ww=0; ww<wn_max; ww++)
  	{
      	G0l[zz][ww]=new Matrix(0);  (*G0l[zz][ww]).new_memory();
      	G0r[zz][ww]=new Matrix(0);  (*G0r[zz][ww]).new_memory();
      	G0lr[zz][ww]=new Matrix(0); (*G0lr[zz][ww]).new_memory();
         Matrix l=Rotate(zz), r(l); r.trans();
         *G0l[zz][ww]=l*G0w_(zz,ww);
         *G0r[zz][ww]=G0w_(zz,ww)*r;
         *G0lr[zz][ww]=l*G0w_(zz,ww)*r;
   ;}
   k=1;
   ;}
	if (f==-1) return G0l[z][wn];
	if (f==1)  return G0r[z][wn];
   return G0lr[z][wn];
;}


Vector * G_asimptotics_mu1 [n_zone], * G_asimptotics_mu2 [n_zone];
Matrix * G_asimptotics_m1 [n_zone], * G_asimptotics_m2 [n_zone], * G_asimptotics_m3 [n_zone];

void G_asimptotics_define()
{
static int k=0; if (k==1) return; k=1;
	for (int z=0; z<n_zone; z++)
   {
   G_asimptotics_mu1[z]=new Vector();
   G_asimptotics_mu2[z]=new Vector();
   G_asimptotics_m1[z]=new Matrix (0); (*G_asimptotics_m1[z]).new_memory();
   G_asimptotics_m2[z]=new Matrix (0); (*G_asimptotics_m2[z]).new_memory();
   G_asimptotics_m3[z]=new Matrix (0); (*G_asimptotics_m3[z]).new_memory();

   n_type w_max=(2*(wn_max-1)+1)*Pi/beta;
   Matrix Ginf(*G0w(z,(wn_max-1))); //Matrix N(0);
   //*G_asimptotics_m[z]=(0.25*beta*w_max*w_max)*Ginf;

   for (int i=0; i<n_part; i++)
    for (int j=0; j<n_part; j++)
    {
    	if (i!=j)
      {
      	n_type re1=real(Ginf.x[i][j]),im1=imag(Ginf.x[i][j]);

//      	if (fabs(re1)<1e-10)
//         {(*G_asimptotics_m1 [z]).x[i][j]=(0,0);(*G_asimptotics_m2 [z]).x[i][j]=(0,0);(*G_asimptotics_m3 [z]).x[i][j]=(0,0);}
//      	else
			if (re1*re1>im1*im1)
      	{
			(*G_asimptotics_m1 [z]).x[i][j]=-re1*w_max*w_max*(((im1*im1)/(re1*re1))+1);
	   	(*G_asimptotics_m2 [z]).x[i][j]=im1*w_max/re1;
         (*G_asimptotics_m3 [z]).x[i][j]=(0,0);
      	;}
         else
         {
			(*G_asimptotics_m1 [z]).x[i][j]=0;
	   	(*G_asimptotics_m2 [z]).x[i][j]=0;
         (*G_asimptotics_m3 [z]).x[i][j]=im1*w_max*w_max*w_max;
         ;}

//      	cout<<i<<"	"<<j<<"	"<<re1<<"	"<<im1<<"	"<<real((*G_asimptotics_m1 [z]).x[i][j])<<"	"<<real((*G_asimptotics_m2 [z]).x[i][j])<<"	"<<"\n"<<flush;
      ;}
      else
      {
      n_type re=real(Ginf.x[i][i]),im=imag(Ginf.x[i][i]), sq=re*re+im*im;
      if (-w_max*(im+w_max*sq)*(1+4*w_max*(im+w_max*sq))>-1e-12)
      {
			(*G_asimptotics_mu1 [z]).x[i]=
			(-w_max*re+sqrt(1e-12-w_max*(im+w_max*sq)*(1+4*w_max*(im+w_max*sq))))/(im+2*w_max*sq);

	   	(*G_asimptotics_mu2 [z]).x[i]=
			(-w_max*re-sqrt(1e-12-w_max*(im+w_max*sq)*(1+4*w_max*(im+w_max*sq))))/(im+2*w_max*sq);
      ;}
      else
      {
      	(*G_asimptotics_mu1 [z]).x[i]=0; (*G_asimptotics_mu2 [z]).x[i]=0;
      	static int h=0; if (h==0) {cout<<"\n!!! problems with high-frequency asymptotics of G - suspect unphysical bath  !!!\n";h=1;}
      ;}
      ;}
    ;}

	;}
   //int rr; cin>>rr;
;}

Matrix & G0_inf_w(int z, n_type w)
{
	G_asimptotics_define();
	//Matrix M((*G_asimptotics_m[z])*(4/(beta*w*w)));
	static Matrix M; M=0;
   for (int i=0; i<n_part; i++)
    for (int j=0; j<n_part; j++)
    {
     if (i!=j)
     {
      M.x[i][j]=(*G_asimptotics_m1[z]).x[i][j]/(I*w*(*G_asimptotics_m2[z]).x[i][j]-(w*w))
               +I*(*G_asimptotics_m3[z]).x[i][j]/(w*w*w);
     ;}
     else
     {
   	M.x[i][i]=n_type(0.5)/(I*w-(*G_asimptotics_mu1[z]).x[i])
  			+n_type(0.5)/(I*w-(*G_asimptotics_mu2[z]).x[i]);
     ;}
    ;}
   return M;
;}

n_type Fermi_dist(n_type mu, n_type t)
{
	if (mu*(beta-t)>200.) return 0.;
	if (mu*t>200.) return exp(mu*(t-beta));
	return 1./(exp(-mu*t)+exp(mu*(beta-t)));
;}

Matrix & G0_inf_t(int z, n_type t)
{
	G_asimptotics_define();
	static Matrix M; M=0;
   for (int i=0; i<n_part; i++)
    for (int j=0; j<n_part; j++)
    {
     if (i!=j)
     {
      if (fabs(real((*G_asimptotics_m2[z]).x[i][j]))<1e-10)
          {M.x[i][j]=real((*G_asimptotics_m1[z]).x[i][j])*(2*t-beta)/4.;
		  }
      else
      {
          M.x[i][j]=real((*G_asimptotics_m1[z]).x[i][j])*(0.5-Fermi_dist(-real((*G_asimptotics_m2[z]).x[i][j]),t))/real((*G_asimptotics_m2[z]).x[i][j]);
      ;}
      M.x[i][j]=M.x[i][j]+(*G_asimptotics_m3[z]).x[i][j]*((sqr(t-beta/2)-sqr(beta/2))/4);
     ;}
     else
     {
   	M.x[i][i]=0.5*Fermi_dist(real((*G_asimptotics_mu1[z]).x[i]),t)
      			+0.5*Fermi_dist(real((*G_asimptotics_mu2[z]).x[i]),t);
     ;}
    ;}
   return M;
;}


//==================
Matrix & G0t(int z, n_type t)
{
   Matrix G0(0);
//	for (n_type w=-(2*wn_max-1)*(Pi/beta); w<(2*wn_max+0.5)*(Pi/beta); w+=2*Pi/beta)
	for (int wn=0; wn<wn_max; wn++)
   {
   	n_type w=(2*wn+1)*Pi/beta;
   	G0+=2*(1/beta)*exp(I*w*t)*(*G0w(z,wn)-G0_inf_w(z,w)); //only real part is OK...
   ;}
   static Matrix res; res=G0+G0_inf_t(z,t);
   return res;
;}

n_type G0t_element(n_type t, int z, int i1, int i2)
{
	return real(G0t(z,t).x[i1][i2]);//very ineffective - do not care...
;}

n_type G0der_0(int z, int i1, int i2)
{
	n_type eps=1e-4*beta;
	return (G0t_element(eps,z,i1,i2)-G0t_element(0,z,i1,i2))/eps;
;}

n_type G0der_beta(int z, int i1, int i2)
{
	n_type eps=1e-4*beta;
	return (G0t_element(beta,z,i1,i2)-G0t_element(beta-eps,z,i1,i2))/eps;
;}

//===================================================================================

Spline *** G0_[n_zone];

void splines_def()
{
	for (int z=0;  z<n_zone;  z++)
   {
   	G0_[z]=new Spline ** [n_part];
   	for (int p1=0; p1<n_part; p1++)
      {
         G0_[z][p1]=new Spline * [n_part];
			for (int p2=0; p2<n_part; p2++)
   		{
   			G0_[z][p1][p2]=new Spline(n_tau);
      		(*G0_[z][p1][p2]).define(G0t_element, 0, beta,G0der_0(z,p1,p2),G0der_beta(z,p1,p2),z,p1,p2);
      	;}
      ;}
   ;}
;}

n_type G0(point & p, point_ & p_)
{
	{static int k=0; if (k==0) splines_def(); k=1;}
	if (p.z!=p_.z) return 0;
   n_type t=p_.t-p.t;
   if (t*t<1e-16) return (*G0_[p.z][p.i][p_.i]).val(0);
   if (t>0) return  (*G0_[p.z][p.i][p_.i]).val(t);
   else return  -(*G0_[p.z][p.i][p_.i]).val(t+beta);
;}



