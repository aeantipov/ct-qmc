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

ComplexType *** GM00_st = new ComplexType ** [n_zone]; 
ComplexType **** Gtotal_st = new ComplexType ***[n_zone];//[WN_max][n_part];

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
	if (enable_clusters==0 || rnd()>num_clust)	MC_step(); 
	else { 
		cluster_step();
		 }

;}


#include "ini.cpp"
			//Initialisation; memory for dynamical arrays; pre-sets of variables
			//and a number of initial MC steps ...

void gw_gtotal_weight_backup(ComplexType ****Gtotal_from, ComplexType ***Gw_from, n_type weight_from, ComplexType ****Gtotal_to, ComplexType ***Gw_to, n_type& weight_to)
/* Performs a backup from received data to container */
{
   weight_to = weight_from;
   for (int z=0; z<n_zone; z++)
   for (int wn=0; wn<wn_max; wn++)
   for (int j=0; j<n_part; j++)
   {
     Gw_to[z][wn][j]=Gw_from[z][wn][j];
     for (int j2=0; j2<n_part; j2++)
	  {
	    Gtotal_to[z][wn][j][j2]=Gtotal_from[z][wn][j][j2];
	  }
  ;}
}

void Gw_write(ComplexType *****Gtotal_, ComplexType ****Gw_, n_type *weight_)
/* Prints Gw.dat and Gw_complete.dat from Gw and Gtotal backup containers */
{
    ofstream ou_gg("Gw.dat"), ou_total("Gw_complete.dat");
    /* calculate total weight */
    double sumweight = 0.;
    for (int i=1; i<CTQMC.getSize();i++) sumweight+=weight_[i-1];

    /* Sum everything */
    for (int i=1; i<CTQMC.getSize();i++)
    for (int z=0; z<n_zone; z++)
    for (int wn=0; wn<wn_max; wn++)
    for (int j=0; j<n_part; j++)
      {if (i==1) GM00_st[z][wn][j]=Gw_[i-1][z][wn][j]/sumweight;
       else      GM00_st[z][wn][j]+=Gw_[i-1][z][wn][j]/sumweight;
       for (int j2=0; j2<n_part; j2++)
       {
         if (i==1) Gtotal_st[z][wn][j][j2]=Gtotal_[i-1][z][wn][j][j2]/sumweight;
         else      Gtotal_st[z][wn][j][j2]+=Gtotal_[i-1][z][wn][j][j2]/sumweight;
       }
      };

  //    cout << "Total weight: " << sumweight << endl;

     /*Print output*/
      for (int wn=0; wn<wn_max; wn++)
      {
      ou_gg<<wn<<"  ";  ou_total<<wn<<"  ";
      for (int z=0; z<n_zone; z++)
      for (int j=0; j<n_part; j++)
      {
        ou_gg<<double(real(GM00_st[z][wn][j]))<<"  "<<double(imag(GM00_st[z][wn][j]))<<"   ";
        for (int j2=0; j2<n_part; j2++) ou_total<<double(real(Gtotal_st[z][wn][j][j2]))<<"  "<<double(imag(Gtotal_st[z][wn][j][j2]))<<"   ";
      ;}
      ou_gg<<"\n"; ou_total<<"\n";
      ;}         
}

int main_routines()
{   
    cout<<"\n       This is CT-QMC code version 1.6\n       Check www.ct-qmc.org for latest updates.\n";
    cout<<"\n       Using MPI Version.\n";
    cout<<"       N of Markov walkers = " << CTQMC.getSize()-1 << ".\n\n";
    cout << "Preparing alpha_W and G0 in background." << endl; 



   ofstream st("stop.dat"); st<<2<<flush;
   int output_period=int_value("output_period");
   int output_period_time=int_value("output_period_time");

   /* Prepare input containers and allocate memory */

   cout << "Preparing input_containers and allocating memory." << endl; 

   ComplexType **** Gw0_container = new ComplexType *** [CTQMC.getSize()-1];       // Receive container for local Gw
   ComplexType ***** Gtotal_container = new ComplexType **** [CTQMC.getSize()-1];  // Input container for total Gw
   n_type * weights = new n_type [CTQMC.getSize()-1];

   ComplexType **** Gw0_container_backup = new ComplexType *** [CTQMC.getSize()-1];       // Container for local Gw
   ComplexType ***** Gtotal_container_backup = new ComplexType **** [CTQMC.getSize()-1];  // Container for total Gw
   n_type * weights_backup = new n_type [CTQMC.getSize()-1];
   int * finish_recv = new int [CTQMC.getSize()-1];


   MPI_Request **all_Gw = new MPI_Request* [CTQMC.getSize()-1];
   MPI_Request **all_Gtotal = new MPI_Request* [CTQMC.getSize()-1];
   MPI_Request *all_weight = new MPI_Request [CTQMC.getSize()-1];
   MPI_Request *all_finish = new MPI_Request [CTQMC.getSize()-1];


   for (int i=1;i<CTQMC.getSize();i++)
     {
       Gw0_container[i-1] = new ComplexType ** [n_zone];
       Gtotal_container[i-1] = new ComplexType *** [n_zone];
       Gw0_container_backup[i-1] = new ComplexType ** [n_zone];
       Gtotal_container_backup[i-1] = new ComplexType *** [n_zone];
       for (int z=0;z<n_zone;z++)
         {
            Gw0_container[i-1][z] = new ComplexType * [wn_max];
            Gtotal_container[i-1][z] = new ComplexType ** [wn_max];
            Gw0_container_backup[i-1][z] = new ComplexType * [wn_max];
            Gtotal_container_backup[i-1][z] = new ComplexType ** [wn_max];
            for (int w=0;w<wn_max;w++)
	      {
                 Gw0_container[i-1][z][w] = new ComplexType [n_part];
                 Gtotal_container[i-1][z][w] = new ComplexType * [n_part];
                 Gw0_container_backup[i-1][z][w] = new ComplexType [n_part];
                 Gtotal_container_backup[i-1][z][w] = new ComplexType * [n_part];
		 for (int j=0;j<n_part;j++)
		   {
                     Gw0_container_backup[i-1][z][w][j] = 0.;
                     Gtotal_container[i-1][z][w][j] = new ComplexType [n_part];
                     Gtotal_container_backup[i-1][z][w][j] = new ComplexType [n_part];
		     for (int j2=0;j2<n_part;j2++) Gtotal_container_backup[i-1][z][w][j][j2] = 0.;
		   }
	      }
	 }
     };

   cout << "\t Memory allocation for send/recv routines" << endl;
   int mem = 2*sizeof(ComplexType)*(CTQMC.getSize()-1)*n_zone*wn_max*n_part*n_part; int fullmem = mem;
   cout << "Allocated " << mem / 1024 << " Kb for total G container" << endl;
   mem = 2*sizeof(ComplexType)*(CTQMC.getSize()-1)*n_zone*wn_max*n_part; fullmem += mem;
   cout << "Allocated " << mem / 1024 << " Kb for local Gw container" << endl;

   for (int z=0;z<n_zone;z++)
         {
            GM00_st[z] = new ComplexType * [wn_max];
            Gtotal_st[z] = new ComplexType ** [wn_max];
            for (int w=0;w<wn_max;w++)
	      {
                 GM00_st[z][w] = new ComplexType [n_part];
                 Gtotal_st[z][w] = new ComplexType * [n_part];
		 for (int j=0;j<n_part;j++)
		   {
                     Gtotal_st[z][w][j] = new ComplexType [n_part];
		   }
	      }
	 };

   mem= sizeof(ComplexType)*n_zone*wn_max*n_part*(1+n_part); fullmem += mem;

   /* Prepare requests for receiving data */
   for (int i=1;i<CTQMC.getSize();i++)
     {
       all_Gw[i-1]=CTQMC.irecv(i, Gw0_container[i-1],n_zone,wn_max,n_part);
       all_Gtotal[i-1]=CTQMC.irecv(i,Gtotal_container[i-1],n_zone,wn_max,n_part,n_part);
       all_weight[i-1]=CTQMC.irecv_float(i,weights[i-1]);
   	   all_finish[i-1]=CTQMC.irecv_int(i,finish_recv[i-1],FINISH_MESSAGE);
	   finish_recv[i-1]=0;
     };
	int AmountOfFinishedProcesses=0;
    
   mem= sizeof(MPI_Request)*CTQMC.getSize()*(n_zone*wn_max*(n_part+1)*2+1); fullmem += mem;
   cout << "Allocated " << mem / 1024 << " Kb for MPI requests" << endl;
 //  cout << "Totally allocated " << fullmem / 1024 << " Kb " << endl;

   cout << "Syncing..." << flush; 
   CTQMC.sync();
   cout << "Done." << endl; 
   cout << "Walkers started." << endl; 

   int time_step=0; // Timing in MASTER Process

   int check_gw_received[CTQMC.getSize()-1]; // Check flag for receiving Gw. If check_gw_received = N of walkers then some data is received
   for (int i=1;i<CTQMC.getSize();i++) check_gw_received[i-1]=0;

   int check_finished[CTQMC.getSize()-1]; // Check flag for walker finish
   for (int i=1;i<CTQMC.getSize();i++) check_finished[i-1]=0;

//   int tmp; std::cin>>tmp;
   for (;;)
   {
   time_step++;

	for (int i=1;i<CTQMC.getSize();i++){
		/* Test finish of process N i */
		if (CTQMC.test(all_finish[i-1]) ){
	    	check_finished[i-1] = 1;
			while (finish_recv[i-1]<1) {};

            if (finish_recv[i-1]==1){ 
				AmountOfFinishedProcesses++; 
				cout << endl << "Walker " << i << " finished. " << AmountOfFinishedProcesses << " finished in total." << endl; 
				}

	    	finish_recv[i-1]=2;
  	  		};     
	    if (!finish_recv[i-1]) 
			if (CTQMC.testall(all_Gtotal[i-1],n_zone*wn_max*n_part))
        		if (CTQMC.testall(all_Gw[i-1],n_zone*wn_max))
           			if (CTQMC.test(all_weight[i-1])){
	       				gw_gtotal_weight_backup(Gtotal_container[i-1], Gw0_container[i-1], weights[i-1], Gtotal_container_backup[i-1], Gw0_container_backup[i-1],weights_backup[i-1]);
	       				check_gw_received[i-1]=1; // Check flag for data receipt
	       				delete[] all_Gw[i-1]; delete[] all_Gtotal[i-1]; 
	
			 	       	cout << "Received data from process " << i << " " << endl;
	
   	    	       		all_Gw[i-1]=CTQMC.irecv(i, Gw0_container[i-1],n_zone,wn_max,n_part);
   	            		all_Gtotal[i-1]=CTQMC.irecv(i,Gtotal_container[i-1],n_zone,wn_max,n_part,n_part);
   	            		all_weight[i-1]=CTQMC.irecv_float(i,weights[i-1]);
		     			};
		} // end of for (int i=1;i<CTQMC.getSize();i++)
	
	// check finish
    int check_finish=0;
   	for (int i=1;i<CTQMC.getSize();i++) {check_finish+=check_finished[i-1];}
   	if (check_finish==CTQMC.getSize()-1)
    	{
	   		cout << "Exiting..." << flush;
       		return 1;
     	}

	else {
	// Check Gw write
   	int check_gw=0;
   	for (int i=1;i<CTQMC.getSize();i++) {check_gw+=check_gw_received[i-1];}
   	static int curr_t=time(NULL)/output_period_time;
   	if ( (time_step%(20*output_period)==0 && time_step!=0 && check_gw > 0 ) || (check_gw==CTQMC.getSize() - AmountOfFinishedProcesses )|| (time_step%(100000000000*output_period)==0 && time_step!=0)){ 
    	Gw_write(Gtotal_container_backup, Gw0_container_backup, weights_backup);
     	check_gw=0;
    	for (int i=1;i<CTQMC.getSize();i++) check_gw_received[i-1]=0;
	 	cout << "Dumped Gw" << endl;
   		} // end::if for Gw dump
	};
	
}; // end of eternal loop
} // end of all main_routines() procedure


int worker_routines()
{
   CTQMC.getStream()<<n_zone<<"    "<<n_part<<"\n"<<flush;
   n_type s=0; int ii=0, ti=time(NULL);


//   ofstream history("history.dat");

   int scratch_period=int_value("scratch_period");
   int output_period=int_value("output_period");
   int output_period_time=int_value("output_period_time");

   int i_max=alphaW_number;
   int const N_autocorr=int_value("N_autocorr"),N_autocorr_max=500, N_autocorr_factor=20;//int(1+beta*beta*U*U/10);
   n_type autocorr[N_autocorr_max], autocorr_aux[N_autocorr_max];
   n_type autocorr_s[N_autocorr_max], autocorr_s_aux[N_autocorr_max];
   {for (int i=0; i<N_autocorr; i++) {autocorr[i]=0;autocorr_aux[i]=0;autocorr_s[i]=0;autocorr_s_aux[i]=0;};}


CTQMC.sync();
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
//if ( (i%output_period==0 && i!=0 ) || (i%(100*output_period)==0 && i!=0))
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
            {
	    /*================== Gw ================= */
		CTQMC.isend(0,Gtotal_st,n_zone,wn_max,n_part,n_part);
 		CTQMC.isend(0,GM00_st,n_zone,wn_max,n_part);
		CTQMC.isend(0,weight_sum);

		/*
		// test
		stringstream s; 
		s << "Gw_proc" << CTQMC.getRank() << ".dat";
		stringstream s2; 
		s2 << "Gw_complete_proc" << CTQMC.getRank() << ".dat";
		ofstream ou_gg(s.str().c_str());
		ofstream ou_total(s2.str().c_str());
                for (int wn=0; wn<wn_max; wn++)
      		     {
                       ou_gg<<wn<<"  ";  ou_total<<wn<<"  ";
                       for (int z=0; z<n_zone; z++)
                       for (int j=0; j<n_part; j++)
                     {
                       ou_gg<<double(real(GM00_st[z][wn][j]))/weight_sum<<"  "<<double(imag(GM00_st[z][wn][j]))/weight_sum<<"   ";
                       for (int j2=0; j2<n_part; j2++) ou_total<<double(real(Gtotal_st[z][wn][j][j2]))/weight_sum<<"  "<<double(imag(Gtotal_st[z][wn][j][j2]))/weight_sum<<"   ";
    		     ;}
     		       ou_gg<<"\n"; ou_total<<"\n";
     		     ;}         
*/
	    /*================== Gt ================= */
     //          Gt_write();

	    /*================== sigma ================= */
			   /*static int write_sigma=int_value("write_sigma");
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
               ;}*/
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
		 if (CTQMC.getRank() == 1)
		   {
                     cout<<"\n\nESTIMATIONS at "<<i<<" steps, after "<<time(NULL)-ti<<" sec.:\nautocorrelation time "<<int(tau)<<"\n";
                     cout<<"variance "<<double((DispG0_sum-DispG0_sum00)*i/(weight_sum*weight_sum))<<"\n";
              	     cout<<"calc. time "<<int ((DispG0_sum-DispG0_sum00)*tau*(time(NULL)-ti)/sqr(weight_sum*desired_accuracy*beta))<<" sec.\n";
                     cout<<"Average pert. order "<<Nc_sum/weight_sum<<"\n"<<flush;
		   }
		 CTQMC.getStream()<<"\n\nESTIMATIONS at "<<i<<" steps, after "<<time(NULL)-ti<<" sec.:\nautocorrelation time "<<int(tau)<<"\n";
         	 CTQMC.getStream()<<"variance "<<double((DispG0_sum-DispG0_sum00)*i/(weight_sum*weight_sum))<<"\n";
               	 CTQMC.getStream()<<"calc. time "<<int ((DispG0_sum-DispG0_sum00)*tau*(time(NULL)-ti)/sqr(weight_sum*desired_accuracy*beta))<<" sec.\n";
                 CTQMC.getStream()<<"Average pert. order "<<Nc_sum/weight_sum<<endl;
                  jj=1;
            	;}
              ;}

	    /* Check for stop.dat */
            //  ifstream st; st.open("stop.dat"); int k=0; st>>k; st.close();
	     int k=2;

	   /* Output condition */
              if (k!=2 || i>0.95*i_max || (i>1e5 && sqrt((DispG0_sum-DispG0_sum00)*tau/(weight_sum*weight_sum))<desired_accuracy*beta))
              {
	      CTQMC.getStream() << "K = " << k << " i = " << i << endl;
            	      {
		 if (CTQMC.getRank() == 1)
		   {
                     cout<<"\nRESULTS\n"<<i<<"   steps done in "<<time(NULL)-ti<<" sec.";
                     cout<<"\nautocorrelation time "<<int(tau)<<"\n";
          				cout<<"variance "<<double((DispG0_sum-DispG0_sum00)*i/(weight_sum*weight_sum))<<"\n";
         				cout<<"estim. error in G(w0) "<<double(sqrt((DispG0_sum-DispG0_sum00)*tau/(weight_sum*weight_sum)))<<"\n";//((DispG0_sum/weight_sum)*tau/i)
                     cout<<"Average pert. order "<<Nc_sum/weight_sum;
		   };
                     CTQMC.getStream()<<"\nRESULTS\n"<<i<<"   steps done in "<<time(NULL)-ti<<" sec.";
                     CTQMC.getStream()<<"\nautocorrelation time "<<int(tau)<<"\n";
         	     CTQMC.getStream()<<"variance "<<double((DispG0_sum-DispG0_sum00)*i/(weight_sum*weight_sum))<<"\n";
         	     CTQMC.getStream()<<"estim. error in G(w0) "<<double(sqrt((DispG0_sum-DispG0_sum00)*tau/(weight_sum*weight_sum)))<<"\n";//((DispG0_sum/weight_sum)*tau/i)
                     CTQMC.getStream()<<"Average pert. order "<<Nc_sum/weight_sum << endl;
		      			;}
                     CTQMC.getStream()<<"\nthat's all"<<flush; //cout<<"   "<<(double)t_counter/ CLOCKS_PER_SEC;

	   		 CTQMC.isend(0,Gtotal_st,n_zone,wn_max,n_part,n_part);
			 CTQMC.isend(0,GM00_st,n_zone,wn_max,n_part);
			 CTQMC.isend(0,weight_sum);
		     CTQMC.isend_int(0,1,FINISH_MESSAGE);
			 cout << "Sent finish message from process " << CTQMC.getRank() << endl;

                     return 1;
//            			sgn_sum/=(i+1); sgn_sum*=i_max; i=i_max+1
              ;}
            ;}
            
;}  //end for the "output" block
;}  //end for the main loop
/*CTQMC.isend(0,Gtotal_st,n_zone,wn_max,n_part,n_part);
CTQMC.isend(0,GM00_st,n_zone,wn_max,n_part);
CTQMC.isend(0,weight_sum);
*/
cout << "No way! " << CTQMC.getRank() << endl;
CTQMC.isend_int(0,1,FINISH_MESSAGE);
return 1;
}

int main (int argc, char **argv)
{

   CTQMC.INIT(argc,argv);

   CTQMC.IniProcess();




   int my_rank=CTQMC.getRank();

  switch (my_rank) 
  {
  case 0 : 
  { // MASTER Process 0 
      main_routines();
	  CTQMC.sync();
      cout << "That's all" << endl;
    break;
  }
  default : 
  { // WORKER_PROCESSES
    INT_RANDOM=time(NULL)+my_rank*500;
	try{
    	Ini();
	    CTQMC.getStream() << "Process " << my_rank << " initialized memory." << endl;
	    cout << "Process " << my_rank << " initialized memory." << endl;
	   }
	catch (std::bad_alloc&)
	   {
	     CTQMC.getStream() << "Process " << my_rank << " couldn't allocate memory. Exiting." << endl;
	     cout << "Process " << my_rank << " couldn't allocate memory. Exiting." << endl;
	     exit(0);
	   };
    worker_routines();
	CTQMC.sync();
  }
};

return 1;
;}





