//===============================  from scratch  =====================================

n_type from_scratch()
{
	int n=Nc; n_type w=1;
	Nc=0; for (int z=0; z<n_zone; z++) NumberOfPoints[z]=0;
   if (n%4==0) curr_sgn=1; else curr_sgn=-1;

   {for (int z=0; z<n_zone; z++)
    for (int wn=0; wn<wn_max; wn++)
    for (int nn=0; nn<n_part; nn++)
    {
    	for (int nn2=0; nn2<n_part; nn2++)
      {
         GM_matrix[z][wn][nn][nn2]=(*Grot(z,wn,0)).x[nn][nn2];
      ;}
    ;}
    ;}




   for (int k=0; k<n; k++)
   {
   	w*=add();
      if (Nc%2==0) {w/=(Nc/2);w*=u[Nc-2];}
      if (u[Nc-1]<0) curr_sgn=-curr_sgn;
   ;}
   prev_weight=weight();
	return fabs(w);
;}


