 //optical susceptibility: n(0) n(tau)

 ComplexType * nn_w_data, * ss_w_data, * G_nn11, * G_nn22,* G_nn21,* G_nn12; n_type nn_sum=0;

// ofstream debug("debug.dat");

 void write_nn_files()
 {//debug<<nn_sum<<"\n"<<flush;
 static int wc_max=int_value("number_of_Matsubara_frequencies_for_Gamma4");
 if (fabs(nn_sum)<1e-8) return;
   static ComplexType * nn_w =new ComplexType [wc_max/2];
   static ComplexType * ss_w =new ComplexType [wc_max/2];
   {for (int WNN=0; WNN<wc_max/2; WNN++) {nn_w[WNN]=nn_w_data[WNN]/nn_sum;ss_w[WNN]=ss_w_data[WNN]/nn_sum;};}

 	ofstream nnw_str("nnw.dat"), nntau_str("nntau.dat");
   for (int wn=0; wn<wc_max/2; wn++)
   {
   	nnw_str<<wn;
      	nnw_str<<"   "<<real(nn_w[wn])<<" "<<imag(nn_w[wn])<<"   "<<real(ss_w[wn])<<" "<<imag(ss_w[wn]);
      nnw_str<<"\n";
   ;}
   for (n_type tau=0; tau<beta+1e-8; tau+=beta/n_tau)
   {
   	nntau_str<<tau;
      {
   		n_type nn=real(nn_w[0]);
      	{for (int wn=1; wn<wc_max/4; wn++) nn+=2*real(exp(-I*(2.*Pi*wn)*tau/beta)*nn_w[wn]); }     //still no correct asimptotics!!!!
         nntau_str<<"   "<<nn;
   		n_type ss=real(ss_w[0]);
      	{for (int wn=1; wn<wc_max/4; wn++) ss+=2*real(exp(-I*(2.*Pi*wn)*tau/beta)*ss_w[wn]); }     //still no correct asimptotics!!!!
         nntau_str<<"   "<<ss;

      ;}
      nntau_str<<"\n";
   ;}

   ofstream G_nn("Gnn.dat");
   {
   for (int wn=0; wn<2*wc_max/2; wn++)
      G_nn<<wn-wc_max/2<<"   "<<G_nn11[wn]/nn_sum<<"   "<<G_nn22[wn]/nn_sum<<"   "<<G_nn12[wn]/nn_sum<<"   "<<G_nn21[wn]/nn_sum<<"\n";
   ;}

;}


void update_nn()
 {
 	static int todo=int_value("calculate_nn"); if (todo==0) return;

static int wc_max=int_value("number_of_Matsubara_frequencies_for_Gamma4");
   ini_GMchi();
   static int count=0; count++; if (count%(wc_max/2)!=0) return;  //to mantain N^2 ComplexTypeity



   {
   	static int f=0;   if (f==0)
      {
      	nn_w_data=new ComplexType [wc_max]; ss_w_data=new ComplexType [wc_max]; f=1;
         G_nn11=new ComplexType [wc_max]; G_nn22=new ComplexType [wc_max]; G_nn12=new ComplexType [wc_max]; G_nn21=new ComplexType [wc_max];
         {for (int i=0; i<wc_max/2;i++) {nn_w_data[i]=0.;ss_w_data[i]=0.;G_nn11[i]=0;G_nn22[i]=0;G_nn12[i]=0;G_nn21[i]=0;};}



      ;}
   ;}

   static int n1=int_value("nn_number1");
   static int n2=int_value("nn_number2");
   static int z1=int_value("nn_zone1");
   static int z2=int_value("nn_zone2");

set_GMchi();      //gives Green function: GMchi_matrix


nn_sum+=curr_sgn/prev_weight;


for (int WNN=0;WNN<wc_max/2;WNN++)
{
   //G_11 G_22
   {
   	ComplexType s=0, s1=0, s2=0;
      if (WNN==0)
      {
      	s1=real(GMchi_matrix[z1][0][0][n1][n1]), s2=real(GMchi_matrix[z2][0][0][n2][n2]);
      ;}

   	for (int w1=0; w1<2*wc_max/2; w1++)
      for (int w2_=0; w2_<2*wc_max/2; w2_++)
      {
      	n_type runge=1;
         if (WNN==0&& (fabs(w1-wc_max/2+0.5)>3*wc_max/8 || fabs(w2_-wc_max/2+0.5)>3*wc_max/8))  runge=4;
      	s+=runge*(GMchi_matrix[z1][w1][w1+WNN][n1][n1]+0.*s1)*(GMchi_matrix[z2][w2_+WNN][w2_][n2][n2]+0.*s2);

      ;}
      nn_w_data[WNN]+=s*(curr_sgn/prev_weight)/beta;
   ;}
   //G_12 G_21
   if (z1==z2)
   {
   	ComplexType s12=0;
   	if (n1==n2) s12=real(GMchi_matrix[z1][0][0][n1][n1])/wc_max;
   	for (int w1=0; w1<2*wc_max/2; w1++)        //????? indices
      for (int w2_=0; w2_<2*wc_max/2; w2_++)
      	nn_w_data[WNN]+=(GMchi_matrix[z1][w1+WNN][w2_+WNN][n1][n2]+s12)*(GMchi_matrix[z2][w2_][w1][n2][n1]+s12)*(curr_sgn/prev_weight)/beta;
   ;}
   else
   {
   	ComplexType s1=real(GMchi_matrix[z1][0][0][n1][n1]), s2=real(GMchi_matrix[z2][0][0][n2][n2]);
   	for (int w1=0; w1<2*wc_max/2; w1++)
      for (int w2_=0; w2_<2*wc_max/2; w2_++)     //????? indices
      {
      	n_type runge=1;
        if (fabs(w1-wc_max/2+0.5)>3*wc_max/8 || fabs(w2_-wc_max/2+0.5)>3*wc_max/8)  runge=4;
         ss_w_data[WNN]+=runge*(GMchi_matrix[z1][w1+WNN][w2_+WNN][n1][n1])*(GMchi_matrix[z2][w2_][w1][n2][n2])*(curr_sgn/prev_weight)/beta;
      ;}
   ;}
;}





      for (int w=0;w<2*wc_max/2;w++)
      {
      	G_nn11[w]+=GMchi_matrix[z1][w][w][n1][n1]*(curr_sgn/prev_weight); G_nn22[w]+=GMchi_matrix[z1][w][w][n2][n2]*(curr_sgn/prev_weight);
         if (z1==z2)
         {G_nn21[w]+=GMchi_matrix[z1][w][w][n2][n1]*(curr_sgn/prev_weight); G_nn12[w]+=GMchi_matrix[z1][w][w][n1][n2]*(curr_sgn/prev_weight);}
      ;}




 ;}
