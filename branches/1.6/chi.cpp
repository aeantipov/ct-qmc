 //================= separate subroutines for 4-point correlators =======================
 complex **** GMchi_matrix[n_zone];


//ofstream chi_au("chi_au.dat");

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


void ini_GMchi()
{
   static int f=1;
   static int wc_max=int_value("number_of_Matsubara_frequencies_for_Gamma4");
   if (f==1)                                            									//initialization of all the arrays
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



 void chi(int write_flag=0)
 {
 	static int todo=int_value("calculate_Gamma4"); if (todo==0) return;
   static int wc_max=int_value("number_of_Matsubara_frequencies_for_Gamma4");

	static double s=1e-100;
   static int count=0;

 	static complex ******* chi_array [n_zone][n_zone];
   static complex *** G_chi [n_zone];

   static int f=1;
   ini_GMchi();
   if (f==1)                                            									//initialization of all the arrays
   {
   {
   for (int z=0; z<n_zone; z++)
   {
   G_chi[z]=new complex ** [2*wc_max];
   	for (int w1=0; w1<2*wc_max; w1++)
   	{
         G_chi[z][w1]=new complex * [n_part];
            for (int n=0; n<n_part; n++)
            {
            	G_chi[z][w1][n]= new complex [n_part];
               for (int n2=0; n2<n_part; n2++)  G_chi[z][w1][n][n2]=0;
            ;}
   	;}
   ;}
   ;}



   for (int z1=0; z1<n_zone; z1++)
   for (int z2=0; z2<n_zone; z2++)
   {
   	chi_array[z1][z2]=new complex ******[2*wc_max/2];
      for (int w1=0; w1<2*wc_max/2; w1++)
      {
      	chi_array[z1][z2][w1]=new complex *****[2*wc_max/2];
         for (int w2_=0; w2_<2*wc_max/2; w2_++)
         {
	         chi_array[z1][z2][w1][w2_]=new complex ****[2*wc_max/2];
            for (int W=0; W<2*wc_max/2; W++)
            {
            	chi_array[z1][z2][w1][w2_][W]=new complex ***[n_part];
               for (int n1=0; n1<n_part; n1++)
               {
                  chi_array[z1][z2][w1][w2_][W][n1]=new complex **[n_part];
               	for (int n1_=0; n1_<n_part; n1_++)
                  {
                  	chi_array[z1][z2][w1][w2_][W][n1][n1_]=new complex * [n_part];
                     for (int n2=0; n2<n_part; n2++)
                     {
                     	chi_array[z1][z2][w1][w2_][W][n1][n1_][n2]=new complex [n_part];
               			for (int n2_=0; n2_<n_part; n2_++)  chi_array[z1][z2][w1][w2_][W][n1][n1_][n2][n2_]=0
                     ;}
                  ;}
               ;}
            ;}
         ;}
      ;}
   ;}

   f=0;
   ;}
   int W2=wc_max/2;

 	if (write_flag==1)                                                                     //write to the files
   {
/*   ofstream chi_G("chi_G.dat");         //not so neseccary; can be commented
   for (int w=0; w<2*wc_max; w++)
   {
   	for (int z=0; z<n_zone; z++)
      for (int n1=0; n1<n_part; n1++)
      for (int n2=0; n2<n_part; n2++)
      	chi_G<<real(G_chi[z][w][n1][n2])/s<<"    "<<imag(G_chi[z][w][n1][n2])/s<<"   ";
      chi_G<<"\n";
   ;}
*/

   ofstream chi_f("Gamma4.dat");
   chi_f<<"Re         Im               z1 z2          w1' w1 w2' w2        n1' n1 n2' n2\n";

      for (int z1=0; z1<n_zone; z1++)
      for (int z2=n_zone-1; z2>=0; z2--)
      for (int n1=0; n1<n_part; n1++)
      for (int n1_=0; n1_<n_part; n1_++)
      for (int n2=0; n2<n_part; n2++)
      for (int n2_=0; n2_<n_part; n2_++)
      {
		//chi_f<<z1<<"  "<<z2<<"\n";
      for (int W=0; W<2*wc_max/2; W++)
      {
      for (int w1=0; w1<2*wc_max/2; w1++)
      for (int w2_=0; w2_<2*wc_max/2; w2_++)
      if ( (w1-W2)<wc_max/2  && (w1-W2)>=-wc_max/2  && (w1+W-W2)<wc_max/2  && (w1+W-W2)>=-wc_max/2)
      if ( (w2_-W2)<wc_max/2 && (w2_-W2)>=-wc_max/2 && (w2_+W-W2)<wc_max/2 && (w2_+W-W2)>=-wc_max/2)
      {
      	complex z=chi_array[z1][z2][w1][w2_][W][n1][n1_][n2][n2_]/s;
         if (W==0)
         	z-=(G_chi[z1][w1][n1][n1_]/s)*(G_chi[z2][w2_][n2][n2_]/s);
         if (z1==z2 && w1==w2_)
            z+=(G_chi[z1][w1][n1][n2_]/s)*(G_chi[z1][w1+W][n2][n1_]/s);

         z/=G_chi[z1][w1][n1][n1]/s;
         z/=G_chi[z1][w1+W][n1_][n1_]/s;
         z/=G_chi[z2][w2_+W][n2][n2]/s;
         z/=G_chi[z2][w2_][n2_][n2_]/s;

	   z*=beta;


         chi_f<<real(z)<<" "<<imag(z)<<"        ";
         chi_f<<z1<<" "<<z2<<"              "<<w1-W2<<" "<<w1+W-W2<<" "<<w2_+W-W2<<" "<<w2_-W2<<"            "<<n1<<" "<<n1_<<" "<<n2<<" "<<n2_<<"\n";

         if (W!=0)
         {
	         chi_f<<real(z)<<" "<<-imag(z)<<"        ";
         	chi_f<<z1<<" "<<z2<<"              "<<-1-(w1-W2)<<" "<<-1-(w1+W-W2)<<" "<<-1-(w2_+W-W2)<<" "<<-1-(w2_-W2)<<"            "<<n1<<" "<<n1_<<" "<<n2<<" "<<n2_<<"\n";
         ;}

      ;}
//      chi_f<<"\n";
      ;}
//      chi_f<<"\n";
      ;}
      return;
   ;}


count++; if (count%(n_part*wn_max)!=0) return;   //to mantain N^2 complexity...



set_GMchi();
s+=curr_sgn/prev_weight;
{
		for (int z=0; z<n_zone; z++)
      for (int w=0; w<2*wc_max; w++)
      for (int n1=0; n1<n_part; n1++)
      for (int n2=0; n2<n_part; n2++)
      	G_chi[z][w][n1][n2]+=GMchi_matrix[z][w][w][n1][n2]*(curr_sgn/prev_weight);

//      complex gggg=GMchi_matrix[0][0][0][0][0]*(curr_sgn/prev_weight);
//      chi_au<<real(gggg)<<"  "<<imag(gggg)<<"\n"

;}
   	for (int w1=0; w1<2*wc_max/2; w1++)
		if ((w1-W2)<wc_max/2  && (w1-W2)>=-wc_max/2)
       for (int w2_=0; w2_<2*wc_max/2; w2_++)
  		 if ((w2_-W2)<wc_max/2 && (w2_-W2)>=-wc_max/2)
        for (int W=0; W<2*wc_max/2; W++)
        if ((w1+W-W2)<wc_max/2  && (w1+W-W2)>=-wc_max/2 && (w2_+W-W2)<wc_max/2 && (w2_+W-W2)>=-wc_max/2)
{
      for (int z1=0; z1<n_zone; z1++)
      for (int z2=0; z2<n_zone; z2++)
      for (int n1=0; n1<n_part; n1++)
      for (int n1_=0; n1_<n_part; n1_++)
      for (int n2=0; n2<n_part; n2++)
      for (int n2_=0; n2_<n_part; n2_++)
      {
      	chi_array[z1][z2][w1][w2_][W][n1][n1_][n2][n2_]+=
              	GMchi_matrix[z1][w1][w1+W][n1][n1_]*GMchi_matrix[z2][w2_+W][w2_][n2][n2_]*(curr_sgn/prev_weight);
         if (z1==z2)
         chi_array[z1][z2][w1][w2_][W][n1][n1_][n2][n2_]-=
         		GMchi_matrix[z1][w1][w2_][n1][n2_]*GMchi_matrix[z2][w2_+W][w1+W][n2][n1_]*(curr_sgn/prev_weight);
      ;}

;}


;}



