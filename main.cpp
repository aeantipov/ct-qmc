extern default_dict internal_input;
int alphaW_number=int_value("maximum_MC_steps");
n_type desired_accuracy=value("desired_accuracy");
#include "g_recalc.cpp"
			//effective Green function is calculated;

			// ============== this block ends up with functions =====================
			//void W(point & r1, point_ & r1_, point & r2, point_ & r2_, n_type &u, n_type & a1, n_type & a2);
			//n_type G0(point & p, point_ & p_)
			//Matrix G0w(int z, int wn)
         //Matrix g0w(int z, int wn)


int Nc=0, Nc_z[n_zone];
   		//modified by add_point, remove_point

int GLOBAL_MC_COUNT=0, ACCEPT_MC_COUNT=0, N_WL=0;  n_type prev_weight=0; int hist[N_max];
point p[N_max]; point_ p_[N_max]; n_type u[N_max], a[N_max]; int wl[N_max];
   		//modified by MC_step()

int curr_sgn=1;
n_type ** M [n_zone];
ComplexType * ewt[N_max], * ewt_[N_max];
ComplexType *** GM_matrix[n_zone];
		   //modified by low-level functions from moves.cpp, re-calculated by from_scratch()

			//!!!  ALL global variables are also stored-restored by keep_state() !!!

int WL_flag=int_value("use_Wang_Landau"), WL_max;
n_type WL_factor, WL_weight[N_max]; ifstream WL_stream("WL.dat"); ofstream WL_out("WL1.dat");
			//defined in ini and remains unchanged so far...

ComplexType ** GM00_st[n_zone]; ComplexType *** Gtotal_st [n_zone];//[WN_max][n_part];

n_type sgn_sum=0, weight_sum=0, weight_sum2=0, Nc_sum=0, DispG0_sum=0;
			//statistics

bool debugflag=false;
int t_counter=0;

#include "import.cpp"
			//import of interaction and self-energy from files
         //interaction.dat and Delta.dat

#include "keep_state.cpp"
			//stores-recovers the state; used for a fast roll back after rejected steps
			//contains a single procedure keep_state(...)

//#include "direct.cpp"
			//direct checks; this file is not used in a final build

#include "moves.cpp"
			//elementary moves; ends up with functions add() and remove(...)


#include "sampling.cpp"
			//the sampling procedure itself
			//function MC_step() should be called from outside

#include "cluster.cpp"
			//random walk in the cluster - further development of sampling

#include "scratch.cpp"
			//calculation from a scratch
			//contains a single procedure from_scratch()

#include "chi6.cpp"      //calculation of the three-point irreducible correlator
			//is not implemented yet


#include "Gt.cpp"
			//Green function in time domain

#include "nn.cpp"
			//optical susceptibility

void steps()
{

	static int enable_clusters=int_value("Enable_cluster_updates");
	static n_type num_clust=value("Part_of_cluster_steps");
	if (enable_clusters==0 || rnd()>num_clust)	MC_step(); else cluster_step();

;}


#include "ini.cpp"
			//Initialisation; memory for dynamical arrays; pre-sets of variables
			//and a number of initial MC steps ...



int main ()
{


   ofstream info_stream("info.dat"); info_stream<<n_zone<<"    "<<n_part<<"\n"<<flush;

   Ini();

   n_type s=0; int ii=0, ti=time(NULL);
   {ofstream st("stop.dat"); st<<2<<flush;}



//   ofstream history("history.dat");

   int scratch_period=int_value("scratch_period");
   int output_period=int_value("output_period");
   int output_period_time=int_value("output_period_time");

   int i_max=alphaW_number;
   int const N_autocorr=int_value("N_autocorr"),N_autocorr_max=500, N_autocorr_factor=20;//int(1+beta*beta*U*U/10);
   n_type autocorr[N_autocorr_max], autocorr_aux[N_autocorr_max];
   n_type autocorr_s[N_autocorr_max], autocorr_s_aux[N_autocorr_max];
   {for (int i=0; i<N_autocorr; i++) {autocorr[i]=0;autocorr_aux[i]=0;autocorr_s[i]=0;autocorr_s_aux[i]=0;};}



for (int i=0; i<i_max; i++)
{
 	steps();

 	chi(); update_nn();

if (i%N_autocorr_factor==0)
         {
         {for (int j=N_autocorr-1; j>0; j--) {autocorr_aux[j]=autocorr_aux[j-1];autocorr_s_aux[j]=autocorr_s_aux[j-1];};}
         autocorr_aux[0]=real(GM_matrix[0][0][0][0])*curr_sgn/prev_weight; autocorr_s_aux[0]=curr_sgn;
         if (i>N_autocorr*N_autocorr_factor)	{for (int j=0; j<N_autocorr; j++) {autocorr[j]+=autocorr_aux[j]*autocorr_aux[0];autocorr_s[j]+=autocorr_s_aux[j]*autocorr_s_aux[0];};}
         ;}

//if (i%scratch_period==0) from_scratch();
if (rnd(scratch_period)==0) {GlobalMove();}



static int curr_t=time(NULL)/output_period_time;
if ( (i%output_period==0 && i!=0 && curr_t!=time(NULL)/output_period_time) || (i%(100*output_period)==0 && i!=0))
//write to the files...
         {
             curr_t=time(NULL)/output_period_time;
//hist.dat
            {
            	ofstream ou_hist("hist.dat");
               ou_hist<<"#  "<<i<<"   "<<ACCEPT_MC_COUNT<<"   "<<weight_sum/weight_sum2<<"\n";
               for (int i=0; i<N_max; i++) ou_hist<<i<<"  "<<double(hist[i])<<"\n";
            ;}

//Green function (Gw.dat, Gt.dat, Gw_complete.dat, Sigma.dat)
            ComplexType ggg=0;
            {
            	ofstream ou_gg("Gw.dat"), ou_total("Gw_complete.dat");
               for (int wn=0; wn<wn_max; wn++)
               {
               ou_gg<<wn<<"  ";  ou_total<<wn<<"  ";
               for (int z=0; z<n_zone; z++)
               for (int j=0; j<n_part; j++)
               {
               	ou_gg<<double(real(GM00_st[z][wn][j])/weight_sum)<<"  "<<double(imag(GM00_st[z][wn][j])/weight_sum)<<"   ";
               	if (wn==0 && z==0) ggg+=GM00_st[z][wn][j]/(n_part*weight_sum);
                  for (int j2=0; j2<n_part; j2++)
                     ou_total<<double(real(Gtotal_st[z][wn][j][j2])/weight_sum)<<"  "<<double(imag(Gtotal_st[z][wn][j][j2])/weight_sum)<<"   ";
               ;}
               ou_gg<<"\n"; ou_total<<"\n";
               ;}         

               Gt_write();

               static int write_sigma=int_value("write_sigma");
               if (write_sigma==1)
               {
               ofstream Sigma_str("Sigma.dat");
               for (int wn=0; wn<wn_max; wn++)
               {
               Sigma_str<<wn<<"  ";
               for (int z=0; z<n_zone; z++)
               for (int j=0; j<n_part; j++)
               {
                  ComplexType g=GM00_st[z][wn][j]/weight_sum, g0=(*Grot(z,wn,0)).x[j][j];
               	Sigma_str<<double(real(1./g0-1./g))<<"  "<<double(imag(1./g0-1./g))<<"  ";
               ;}
               Sigma_str<<"\n";
               ;}
               ;}
            ;}


//autocorr.dat
				ofstream autocorr_str("autocorr.dat");
            n_type DispG0_sum00=sqr(real(GM00_st[0][0][0]))/i,BGRG=DispG0_sum00*(i-N_autocorr)/((1.*i)*N_autocorr_factor);
				if (i>N_autocorr*N_autocorr_factor && fabs(sgn_sum)>0)
            {for (int h=0; h<N_autocorr; h++)
            {
            autocorr_str<<double(h*N_autocorr_factor)<<"  ";
            autocorr_str<<double(N_autocorr_factor*(autocorr[h]-BGRG)*i/((i-N_autocorr*N_autocorr_factor)*weight_sum))<<"  ";
            autocorr_str<<N_autocorr_factor*autocorr_s[h]/((i-N_autocorr*N_autocorr_factor))-sqr(sgn_sum/i)<<"\n";
            ;}
            ;}
//chi
				chi(1);
//nn
				write_nn_files();

//screen output
            n_type tau=0; {for (int h=0; h<N_autocorr && ((autocorr[h]-BGRG)/(autocorr[0]-BGRG)>0); h++) tau+=N_autocorr_factor*(autocorr[h]-BGRG)/(autocorr[0]-BGRG);}
            {
              {
              static int jj=0;
              if (jj==0 && i>1e4 && sqrt((DispG0_sum-DispG0_sum00)*tau/(10*weight_sum*weight_sum))<desired_accuracy*beta)
            	{
                  cout<<"\n\nESTIMATIONS at "<<i<<" steps, after "<<time(NULL)-ti<<" sec.:\nautocorrelation time "<<int(tau)<<"\n";
         			cout<<"variance "<<double((DispG0_sum-DispG0_sum00)*i/(weight_sum*weight_sum))<<"\n";
               	cout<<"calc. time "<<int ((DispG0_sum-DispG0_sum00)*tau*(time(NULL)-ti)/sqr(weight_sum*desired_accuracy*beta))<<" sec.\n";
                  cout<<"Average pert. order "<<Nc_sum/weight_sum<<"\n"<<flush;
                  jj=1;
            	;}
              ;}
              ifstream st("stop.dat"); int k=0; st>>k;
              if (k!=2 || i>0.95*i_max || (i>1e5 && sqrt((DispG0_sum-DispG0_sum00)*tau/(weight_sum*weight_sum))<desired_accuracy*beta))
              {
            	      {
                     cout<<"\nRESULTS\n"<<i<<"   steps done in "<<time(NULL)-ti<<" sec.";
                     cout<<"\nautocorrelation time "<<int(tau)<<"\n";
          				cout<<"variance "<<double((DispG0_sum-DispG0_sum00)*i/(weight_sum*weight_sum))<<"\n";
         				cout<<"estim. error in G(w0) "<<double(sqrt((DispG0_sum-DispG0_sum00)*tau/(weight_sum*weight_sum)))<<"\n";//((DispG0_sum/weight_sum)*tau/i)
                     cout<<"Average pert. order "<<Nc_sum/weight_sum;

                     info_stream<<"\nRESULTS\n"<<i<<"   steps done in "<<time(NULL)-ti<<" sec.";
                     info_stream<<"\nautocorrelation time "<<int(tau)<<"\n";
         				info_stream<<"variance "<<double((DispG0_sum-DispG0_sum00)*i/(weight_sum*weight_sum))<<"\n";
         				info_stream<<"estim. error in G(w0) "<<double(sqrt((DispG0_sum-DispG0_sum00)*tau/(weight_sum*weight_sum)))<<"\n";//((DispG0_sum/weight_sum)*tau/i)
                     info_stream<<"Average pert. order "<<Nc_sum/weight_sum;
		      			;}
                     cout<<"\nthat's all"<<flush; //cout<<"   "<<(double)t_counter/ CLOCKS_PER_SEC;
                     return 1;
//            			sgn_sum/=(i+1); sgn_sum*=i_max; i=i_max+1
              ;}
            ;}
            
;}  //end for the "output" block
;}  //end for the main loop

return 1;
;}





