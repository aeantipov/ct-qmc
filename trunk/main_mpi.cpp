extern default_dict internal_input;

CTQMC_WORLD CTQMC;


int alphaW_number=int_value("maximum_MC_steps");
n_type desired_accuracy=value("desired_accuracy");
#include "g_recalc.cpp"
			//effective Green function is calculated;

			// ============== this block ends up with functions =====================
			//void W(point & r1, point_ & r1_, point & r2, point_ & r2_, n_type &u, n_type & a1, n_type & a2);
			//n_type G0(point & p, point_ & p_)
			//Matrix G0w(int z, int wn)
         //Matrix g0w(int z, int wn)


int Nc=0, NumberOfPoints[n_zone];
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

std::vector<ComplexType> dataGM00(n_zone*wn_max*n_part,0);
std::vector<ComplexType> dataGtotal(n_zone*wn_max*n_part*n_part,0);
std::vector<ComplexType> dataChi4;


ComplexType *** GM00_st = new ComplexType ** [n_zone]; 
ComplexType **** Gtotal_st = new ComplexType ***[n_zone];//[WN_max][n_part];

n_type sgn_sum=0, weight_sum=0, weight_sum2=0, Nc_sum=0, DispG0_sum=0;
			//statistics

bool debugflag=false;
int t_counter=0;

static const CTQMC_TAG GM00_TAG=CTQMC_TAG(1229);   // arbitrary number to define unique id for this datatype transmission
static const CTQMC_TAG GTOTAL_TAG=CTQMC_TAG(1745); // the same

/* include operations with measured data */
#include "MeasuredData.h"

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
	if (enable_clusters==0 || rnd()>num_clust)	MC_step(); 
	else { 
		cluster_step();
		 }

;}


#include "ini.cpp"
			//Initialisation; memory for dynamical arrays; pre-sets of variables
			//and a number of initial MC steps ...

void print_logo()
{   
   int WorkersAmount = CTQMC.getSize()-1; // optimization

   std::cout << "\n       This is CT-QMC code version 1.6.2\n       Check www.ct-qmc.org for latest updates.\n";
   std::cout << "\n       Using MPI Version.\n";
   std::cout << "       N of Markov walkers = " << CTQMC.getSize() << ".\n\n";
   std::cout << "Preparing alpha_W and G0 in background." << endl; 
   std::cout << "Preparing data containers and allocating memory." << endl; 
}

int worker_routines()
{
CTQMC.getStream()<<n_zone<<"    "<<n_part<<"\n"<<flush;
n_type s=0; int ii=0, ti=time(NULL);

int scratch_period=int_value("scratch_period");
int output_period=int_value("output_period");
int output_period_time=int_value("output_period_time");

/* Define autocorrelation function to check convergence */
int i_max=alphaW_number;
int const N_autocorr=int_value("N_autocorr"),N_autocorr_max=500, N_autocorr_factor=20;//int(1+beta*beta*U*U/10);
n_type autocorr[N_autocorr_max], autocorr_aux[N_autocorr_max];
n_type autocorr_s[N_autocorr_max], autocorr_s_aux[N_autocorr_max];
for (int i=0; i<N_autocorr; i++) {autocorr[i]=0;autocorr_aux[i]=0;autocorr_s[i]=0;autocorr_s_aux[i]=0;};

CTQMC.sync();
for (int i=0; i<i_max; i++){
    /* Print amount of steps done */
    if (!(i%output_period)) cout << "Walker " << CTQMC.getRank() << " : " << i << " steps done." << endl;
    /* Do the actual steps */
    steps();
    chi(); update_nn();

    /* Measure autocorrelation function to check convergence */
    if (i%N_autocorr_factor==0){
        for (int j=N_autocorr-1; j>0; j--)  {autocorr_aux[j]=autocorr_aux[j-1];autocorr_s_aux[j]=autocorr_s_aux[j-1];};
        autocorr_aux[0]=real(GM_matrix[0][0][0][0])*curr_sgn/prev_weight; autocorr_s_aux[0]=curr_sgn;
        if (i>N_autocorr*N_autocorr_factor)	{for (int j=0; j<N_autocorr; j++) {autocorr[j]+=autocorr_aux[j]*autocorr_aux[0];autocorr_s[j]+=autocorr_s_aux[j]*autocorr_s_aux[0];};}
        };

    /* Flush the Markov chain at random time define by scratch_period */
    //if (i%scratch_period==0) from_scratch();
    if (rnd(scratch_period)==0) {GlobalMove();} 

    /* Output data to files */
    static int curr_t=time(NULL)/output_period_time;
    if ( (i%output_period==0 && i!=0 && curr_t!=time(NULL)/output_period_time) || (i%(100*output_period)==0 && i!=0)){
    //if ( (i%output_period==0 && i!=0 ) || (i%(100*output_period)==0 && i!=0))
		stringstream s;
        curr_t=time(NULL)/output_period_time;

        /* Print history of current walker to hist#.dat */
        /* Omitting - useless in most cases 
        s << "hist" << CTQMC.getRank() << ".dat";
		ofstream ou_hist(s.str().c_str());
        s.str("");
        ou_hist<<"#  "<<i<<"   "<<ACCEPT_MC_COUNT<<"   "<<weight_sum/weight_sum2<<"\n";
        for (int j=0; j<N_max; i++) ou_hist<<"  "<<double(hist[i])<<"\n";
        ou_hist.close();
        */

	    /*================== Gw ================= */
			
		s << "Gw_proc" << CTQMC.getRank() << ".dat";
		ofstream ou_gg(s.str().c_str());
        s.str("");
		s << "Gw_complete_proc" << CTQMC.getRank() << ".dat";
		ofstream ou_total(s.str().c_str());
        s.str("");
        for (int wn=0; wn<wn_max; wn++){
            ou_gg<<wn<<"  ";  ou_total<<wn<<"  ";
            for (int z=0; z<n_zone; z++)
              for (int j=0; j<n_part; j++){
                ou_gg<<double(real(GM00_st[z][wn][j]))/weight_sum<<"  "<<double(imag(GM00_st[z][wn][j]))/weight_sum<<"   ";
                for (int j2=0; j2<n_part; j2++) 
                    ou_total<<double(real(Gtotal_st[z][wn][j][j2]))/weight_sum<<"  "<<double(imag(Gtotal_st[z][wn][j][j2]))/weight_sum<<"   ";
    		    };
     		  ou_gg<<"\n"; ou_total<<"\n";
     		  }         
        ou_gg.close();
		ou_total.close();

	    /* Omitted Gt writing for each walker - use green_parse utility */
        //Gt_write();

        //autocorr.dat
        s << "autocorr" << CTQMC.getRank() << ".dat";
		ofstream autocorr_str(s.str().c_str());
        s.str("");
        n_type DispG0_sum00=sqr(real(GM00_st[0][0][0]))/i,BGRG=DispG0_sum00*(i-N_autocorr)/((1.*i)*N_autocorr_factor);
		if (i>N_autocorr*N_autocorr_factor && fabs(sgn_sum)>0){
            for (int h=0; h<N_autocorr; h++){
                autocorr_str<<double(h*N_autocorr_factor)<<"  ";
                autocorr_str<<double(N_autocorr_factor*(autocorr[h]-BGRG)*i/((i-N_autocorr*N_autocorr_factor)*weight_sum))<<"  ";
                autocorr_str<<N_autocorr_factor*autocorr_s[h]/((i-N_autocorr*N_autocorr_factor))-sqr(sgn_sum/i)<<"\n";
                };
            };
//chi
	    chi(1);
//nn
		write_nn_files();

//screen output
        n_type tau=0; 
        for (int h=0; h<N_autocorr && ((autocorr[h]-BGRG)/(autocorr[0]-BGRG)>0); h++) tau+=N_autocorr_factor*(autocorr[h]-BGRG)/(autocorr[0]-BGRG);
        static int jj=0;
        if (jj==0 && i>1e4 && sqrt((DispG0_sum-DispG0_sum00)*tau/(10*weight_sum*weight_sum))<desired_accuracy*beta){
		    if (CTQMC.getRank() == 1){
                cout<<"\n\nESTIMATIONS at "<<i<<" steps, after "<<time(NULL)-ti<<" sec.:\nautocorrelation time "<<int(tau)<<"\n";
                cout<<"variance "<<double((DispG0_sum-DispG0_sum00)*i/(weight_sum*weight_sum))<<"\n";
              	cout<<"calc. time "<<int ((DispG0_sum-DispG0_sum00)*tau*(time(NULL)-ti)/sqr(weight_sum*desired_accuracy*beta))<<" sec.\n";
                cout<<"Average pert. order "<<Nc_sum/weight_sum<<"\n"<<flush;
		    };
		    CTQMC.getStream()<<"\n\nESTIMATIONS at "<<i<<" steps, after "<<time(NULL)-ti<<" sec.:\nautocorrelation time "<<int(tau)<<"\n";
            CTQMC.getStream()<<"variance "<<double((DispG0_sum-DispG0_sum00)*i/(weight_sum*weight_sum))<<"\n";
            CTQMC.getStream()<<"calc. time "<<int ((DispG0_sum-DispG0_sum00)*tau*(time(NULL)-ti)/sqr(weight_sum*desired_accuracy*beta))<<" sec.\n";
            CTQMC.getStream()<<"Average pert. order "<<Nc_sum/weight_sum<<endl;
            jj=1;
         	;}
/* Output condition */
        if (i>0.95*i_max || (i>1e5 && sqrt((DispG0_sum-DispG0_sum00)*tau/(weight_sum*weight_sum))<desired_accuracy*beta)){
	        CTQMC.getStream() << "Exiting at step " << i << endl;
		    if (CTQMC.getRank() == 0){
                cout<<"\nRESULTS\n"<<i<<"   steps done in "<<time(NULL)-ti<<" sec.";
                cout<<"\nautocorrelation time "<<int(tau)<<"\n";
          		cout<<"variance "<<double((DispG0_sum-DispG0_sum00)*i/(weight_sum*weight_sum))<<"\n";
         		cout<<"estim. error in G(w0) "<<double(sqrt((DispG0_sum-DispG0_sum00)*tau/(weight_sum*weight_sum)))<<"\n";//((DispG0_sum/weight_sum)*tau/i)
                cout<<"Average pert. order "<<Nc_sum/weight_sum;
		        };
      	    sgn_sum/=(i+1); sgn_sum*=i_max; i=i_max+1;
            CTQMC.getStream()<<"\nRESULTS\n"<<i<<"   steps done in "<<time(NULL)-ti<<" sec.";
            CTQMC.getStream()<<"\nautocorrelation time "<<int(tau)<<"\n";
         	CTQMC.getStream()<<"variance "<<double((DispG0_sum-DispG0_sum00)*i/(weight_sum*weight_sum))<<"\n";
         	CTQMC.getStream()<<"estim. error in G(w0) "<<double(sqrt((DispG0_sum-DispG0_sum00)*tau/(weight_sum*weight_sum)))<<"\n";//((DispG0_sum/weight_sum)*tau/i)
            CTQMC.getStream()<<"Average pert. order "<<Nc_sum/weight_sum << endl;
            CTQMC.getStream()<<"\nthat's all"<<flush; //cout<<"   "<<(double)t_counter/ CLOCKS_PER_SEC;

            ;}
            
        ;}  //end for the "output" block
    ;}  //end for the main loop
cout << "Walker " << CTQMC.getRank() << " has finished." << endl;
return 1;
}

int main (int argc, char **argv)
{
   
   CTQMC.INIT(argc,argv);
   CTQMC.IniProcess();
   int my_rank=CTQMC.getRank();

   INT_RANDOM=time(NULL)+my_rank*500;
   if (my_rank==0) print_logo();
	try{
    	Ini();
	    CTQMC.getStream() << "Walker " << my_rank << " initialized memory." << endl;
        if (int_value("calculate_Gamma4")) 
            { 
                int wc_max=int_value("number_of_Matsubara_frequencies_for_Gamma4");
                dataChi4.assign(n_zone*n_zone*wc_max*wc_max*2*wc_max*n_part*n_part*n_part*n_part,0.0); // [z1][z2][w1][w2_][W][n1][n1_][n2][n2_]
            }
	    cout << "Walker " << my_rank << " initialized memory." << endl;
	   }
	catch (std::bad_alloc&)
	   {
	     CTQMC.getStream() << "Walker " << my_rank << " couldn't allocate memory. Exiting." << endl;
	     cout << "Walker " << my_rank << " couldn't allocate memory. Exiting." << endl;
	     exit(0);
	   };
    worker_routines();
	CTQMC.sync();

    if (CTQMC.getRank()==0) cout << "Gathering data" << endl;

    MeasuredData ComputedData;
    MeasuredData TotalData;

    ComputedData.setFromComplexData(dataGM00,dataGtotal);
    CTQMC.genericReduce(&ComputedData.dataGM00_real[0], &TotalData.dataGM00_real[0], n_zone*wn_max*n_part, MPI_DOUBLE, MPI_SUM, 0); 
    CTQMC.genericReduce(&ComputedData.dataGM00_imag[0], &TotalData.dataGM00_imag[0], n_zone*wn_max*n_part, MPI_DOUBLE, MPI_SUM, 0); 

    CTQMC.genericReduce(&ComputedData.dataGtotal_real[0], &TotalData.dataGtotal_real[0], n_zone*wn_max*n_part*n_part, MPI_DOUBLE, MPI_SUM, 0); 
    CTQMC.genericReduce(&ComputedData.dataGtotal_imag[0], &TotalData.dataGtotal_imag[0], n_zone*wn_max*n_part*n_part, MPI_DOUBLE, MPI_SUM, 0); 

    CTQMC.genericReduce(&weight_sum, &TotalData.weight_sum, 1, MPI_DOUBLE, MPI_SUM, 0);

    if (CTQMC.getRank()==0) {
        std::cout << "Writing data fo files" << endl;
        ofstream str1("Gw.dat");
        str1 << TotalData.serializeGM00();// << std::flush;
        str1.close();

        str1.open("Gw_complete.dat");
        str1 << TotalData.serializeGtotal();// << std::flush;
        str1.close();
        std::cout << "Done. Exiting..." << std::endl;
    };
    return 1;
;}





