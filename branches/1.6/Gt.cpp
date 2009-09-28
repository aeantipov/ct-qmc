n_type ** mu1_Gt, ** mu2_Gt;


n_type Gw0i_Gt(n_type w, int j, int z)
{
	return  -0.5*w/(w*w+mu1_Gt[z][j]*mu1_Gt[z][j])-0.5*w/(w*w+mu2_Gt[z][j]*mu2_Gt[z][j]);
;}

n_type Gw0r_Gt(n_type w, int j, int z)
{
	return -0.5*mu1_Gt[z][j]/(w*w+mu1_Gt[z][j]*mu1_Gt[z][j])-0.5*mu2_Gt[z][j]/(w*w+mu2_Gt[z][j]*mu2_Gt[z][j]);
;}

n_type Fermi_dist_Gt(n_type mu, n_type t)
{
	return exp(mu*(beta-t))/(1+exp(mu*beta));
;}                              

n_type Gtinf_Gt(n_type t, int j, int z)
{
   return 0.5*Fermi_dist_Gt(mu1_Gt[z][j],t)+0.5*Fermi_dist_Gt(mu2_Gt[z][j],t);
;}

n_type Gt_single(n_type t, int j, int z)
{
	n_type s=0;
   int wn=0;
   for (n_type w=Pi/beta; w<2*Pi*wn_max/beta; w+=2*Pi/beta)
   {
   	s+=2*( (imag(GM00_st[z][wn][j])/weight_sum-Gw0i_Gt(w,j,z))*sin(w*t)+
             (real(GM00_st[z][wn][j])/weight_sum-Gw0r_Gt(w,j,z))*cos(w*t))/beta;
      wn++;
   ;}
   return s-Gtinf_Gt(t,j,z);
;}

void Gt_write()
{
   ofstream Gt_str("Gt.dat");
   n_type w=(2*wn_max-1)*Pi/beta;
   for (int z=0; z<n_zone; z++)
   for (int j=0; j<n_part; j++)
   {
   	n_type re=real(GM00_st[z][wn_max-1][j])/weight_sum, im=imag(GM00_st[z][wn_max-1][j])/weight_sum;
      n_type sq=re*re+im*im, eps=1e-10-w*(im+w*sq)*(1+4*w*(im+w*sq));
      if (eps>0)
      {
      	mu1_Gt[z][j]=(-w*re+sqrt(eps))/(im+2*w*sq);
      	mu2_Gt[z][j]=(-w*re-sqrt(eps))/(im+2*w*sq);
      ;}
      else
      {
      	mu1_Gt[z][j]=0; mu2_Gt[z][j]=0;
      ;}
   ;}


   for (n_type t=0; t<beta+1e-5; t+=beta/256)
   {
   	Gt_str<<t;
      for (int z=0; z<n_zone; z++)
	   for (int j=0; j<n_part; j++)
        Gt_str<<"    "<<Gt_single(t, j, z);
      Gt_str<<"\n";
   ;}
;}

