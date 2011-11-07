// default.h - default constants are put here
#include <vector>
#include <sstream>

struct default_value // A structure to handle default values
{
 string name;
 string value;
};

class default_dict // A dictionary for all default values
{
  vector<default_value> dict;
public:
  default_dict();			//Default constants initialization
  int int_value(const char *name);      //Returns an integer value, which corresponds to value in @name
  n_type value(const char *name);       //Returns a float value, which corresponds to value in @name
  void print() {for (int i=0;i<dict.size();i++) cout << dict[i].name << "\t" << dict[i].value << endl; }
};

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int default_dict::int_value(const char *name)
//Returns an integer value, which corresponds to value in @name
{
  bool found = false;
  int i=0;
  stringstream input;
  input << name;
  string s = input.str();
  while ((!found) && (i<dict.size()))
    {
     if (dict[i].name==s)
       {
        found = true; 
        stringstream value(dict[i].value);
	int out;
	if (value >> out)
	  {
	    return out;
	  }
	else 
	  {
	    cout << "Conversion error of " << s << endl;
	  }
       };
     i++;
    }
  cout << "What the hell I'm doing here?" << endl;
  cout << "ERROR : no such default value : " + s << endl;
  return 0;
}

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

n_type default_dict::value(const char *name)
//Returns a float value, which corresponds to value in @name
{
  bool found = false;
  int i=0;
  stringstream input;
  input << name;
  string s = input.str();
  while ((!found) && (i<dict.size()))
    {
     if (dict[i].name==s)
       {
        found = true; 
        stringstream value(dict[i].value);
	n_type out; 
	if (value >> out) 
	  {
	    return out;
	  }
	else 
	  {
	    cout << "Conversion error of " << s << endl;
	  }
       };
     i++;
    }
  cout << "What the hell I'm doing here?" << endl;
  cout << "ERROR : no such default value : " + s << endl;
  return 0.0;
}

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

default_dict::default_dict()
//Default constants initialization
{
   default_value current;

   //===================================================================================//
   current.name = "desired_accuracy"; 
   current.value = "0.001";
   /* 
    * desired errorbar in G(\tau)
    */
   dict.push_back(current);
   //===================================================================================//
   current.name = "use_Wang_Landau";
   current.value = "1";
   /*
    * indicates how to deal with quantum Wang Landau reweighting
    * A-priory guerss for WL weights can be read from the file "WL.dat" and then modified by the standard WL procedure
    * The modified weights are stored in file "WL1.dat".
    * Normally, use "WL1.dat" as an a-priory guess for the futher runs. 
    * Possible values
    *    0: do not use WL reveighting at all
    *    1: read from the file and modify
    *    (recommended)
    *    2: do not read from the file; start from the same WL weights for all orders and modify this set
    *    (use it if you do not have any a-priory guess)
    *    3: import from the file, but not modify 
    *    If input file is absent, option 1 and 3 coincide with 2 and 0, respectively.
   */
  dict.push_back(current);
  //===================================================================================//
  current.name = "Initial_steps_without_WL";
  current.value = "3000";
  dict.push_back(current);
  //===================================================================================//
  current.name = "U";
  current.value = "1";
  /*
   * Hubbard U
  */
  dict.push_back(current);
   //===================================================================================//

  current.name="read_delta";
  current.value="0";
  /*
   * whether to read self-energy from the external file
   * Possible values
   * 	0: do not read
   * 	1: read from file "Delta.dat" the same self-energy for all zones
   * 	2: read from file "Delta.dat" different self-energies for different zones
  */
  dict.push_back(current);
   //===================================================================================//

  current.name="number_of_Matsubara_frequencies";
  current.value="21";
  dict.push_back(current);
   //===================================================================================//

  current.name="sparse_Matsubaras";
  current.value="0";
  dict.push_back(current);
   //===================================================================================//

  current.name="sparse_Matsubaras_start_from";
  current.value="10";
  dict.push_back(current);
   //===================================================================================//

  current.name="N_tau";
  current.value="128";
  dict.push_back(current);
   //===================================================================================//

  current.name="rotate_basis";
  current.value="0";
  /*
   * whether to rotate basis for the observables; 
   * file "gg.dat" contains only diagonal elements of the Green function; the same information
   * is used in the sampling procedure/ One may want to deal with, say, momentum representation for 
   * the Green function. In this case, write an appropriate file "rotate.dat". 
   * For simple cases, do not bother youself, an put this flag zero.
  */
  dict.push_back(current);
   //===================================================================================//

  current.name="alpha";
  current.value="-0.03";
  /*
   * \alpha for the interaction operator; put it approximately -0.03
  */
  dict.push_back(current);
   //===================================================================================//

  current.name="alpha_offdiagonal";
  current.value="0.0001";
  dict.push_back(current);
   //===================================================================================//

  current.name="Spin_Flips_Only_In_Global_Move";
  current.value="1";
  dict.push_back(current);
   //===================================================================================//

  current.name="expected_occupancy";
  current.value="0.5";
  dict.push_back(current);
   //===================================================================================//

  current.name="maximum_MC_steps";
  current.value="100000000";
  /*
   * maximal number of MC steps to be performed;
   * also used at certain initial stages of the program, particularly alphaW_accumulation...
   */
  dict.push_back(current);
   //===================================================================================//

  current.name="scratch_period";
  current.value="500";
  /*
   * How often do you recalculate from a scratch
   */
  dict.push_back(current);
   //===================================================================================//

  current.name="output_period";
  current.value="100000";
  /*
   * How often results are written to the files, MC trials
   */
  dict.push_back(current);
   //===================================================================================//

  current.name="output_period_time";
  current.value="5";
  /*
   * Minimal period of the update of external files, seconds
  */
   dict.push_back(current);
   //===================================================================================//

  current.name="WL_factor";
  current.value="0.5";
  /*
   * Maximal order in WL weight calculation is WL_factor*beta*U*n_part
  */
  dict.push_back(current);
   //===================================================================================//

  current.name="WL_initial_factor";
  current.value="0.01";
   dict.push_back(current);
   //===================================================================================//

  current.name="WL_circles_number";
  current.value="4";
   dict.push_back(current);
   //===================================================================================//

  current.name="Additional_Delta_WL";
  current.value="0";
   dict.push_back(current);
   //===================================================================================//

  current.name="WL_Tolerance_factor";
  current.value="0.8";
  /*
   * Fine tune of the procedure of WL weights determination. 
   * See ini.cpp for details...
  */
   dict.push_back(current);
   //===================================================================================//

  current.name="N_autocorr";
  current.value="50";
  /*
   * Number of points in the plot for autocorrelation function ("autocorr.dat")
  */
   dict.push_back(current);
   //===================================================================================//

  current.name="Use_Global_moves";
  current.value="0";
   dict.push_back(current);
   //===================================================================================//

  current.name="W_group_generators";
  current.value="0";
  /*
   * Program works correctly with Use_Global_moves=1 and W_group_generators=0 only if the
   * interaction is spin-indepent. Be careful!
  */
   dict.push_back(current);
   //===================================================================================//

  current.name="calculate_Gamma4";
  current.value="0";
  /*
   * whether to calculate 4-point correlators
  */
   dict.push_back(current);
   //===================================================================================//

  current.name="calculate_Gamma60";
  current.value="0";
  /*
   * whether to calculate 6-point 4-frequency correlators
  */
   dict.push_back(current);
   //===================================================================================//

  current.name="calculate_Gamma6";
  current.value="0";
  /*
   * whether to calculate 6-point 4-frequency correlators
  */
   dict.push_back(current);
   //===================================================================================//

  current.name="chi4_numerical_zero";
  current.value="1e-10";
   dict.push_back(current);
   //===================================================================================//
  current.name="chi6_numerical_zero";
  current.value="1e-10";
   dict.push_back(current);
   //===================================================================================//

  current.name="number_of_Matsubara_frequencies_for_Gamma4";
  current.value="12";
  /*
   * Doubled number of Matsubara's for 4-point correlators
  */
   dict.push_back(current);
   //===================================================================================//
  current.name="Enable_cluster_updates";
  current.value="0";
   dict.push_back(current);
   //===================================================================================//

  current.name="cluster_size";
  current.value="4";
   dict.push_back(current);
   //===================================================================================//

  current.name="number_of_trials_in_cluster";
  current.value="6";
   dict.push_back(current);
   //===================================================================================//

  current.name="Part_of_cluster_steps";
  current.value="0.1";
   dict.push_back(current);
   //===================================================================================//

  current.name="number_of_fields";
  current.value="1";
   dict.push_back(current);
   //===================================================================================//

  current.name="calculate_nn";
  current.value="0";
   dict.push_back(current);
   //===================================================================================//

  current.name="nn_number1";
  current.value="0";
   dict.push_back(current);
   //===================================================================================//

  current.name="nn_number2";
  current.value="0";
   dict.push_back(current);
   //===================================================================================//

  current.name="nn_zone1";
  current.value="0";
   dict.push_back(current);
   //===================================================================================//

  current.name="nn_zone2";
  current.value="0";
   dict.push_back(current);
   //===================================================================================//

  current.name="write_sigma";
  current.value="0";
  /*
   * whether to write "Sigma.dat"
  */
   dict.push_back(current);
   //===================================================================================//



;}





