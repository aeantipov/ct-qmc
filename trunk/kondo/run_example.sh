#!/bin/bash
#===================================================================================
#
# FILE: run_example.sh
#
# USAGE: run_example.sh 
#
# DESCRIPTION: Implementation of the self-consistency loop in a dmft or dual fermion 
# calculation as used in master's thesis by Jonas Sweep. 
#
# kondo22 is the ct-qmc calculation procedure of impurity Green's function (Gw_complete.dat)
# using the s-d model/Kondo Lattice model Hamiltonian. dual21 is the calculation procedure 
# of the hybridisation (Delta.dat) in either dmft or dual fermion formalism depending on the 
# flags `calculate_Gamma4' and `DMFT_sigma'.
#
# OPTIONS: none
# REQUIREMENTS: precompiled ct-qmc-kondo(kondo22), change(../change), submit script(submit_kondo.sh)
# BUGS: ---
# NOTES: ---
# AUTHOR: Jonas Sweep, (jonas*sweep_[at]gmail.com, non-bots, please omit the underscore and asterisk) 
# CREATED: April 2012
# REVISION: 25.04.2012 
#===================================================================================

rm -rfv Gw_complete.dat Delta_old.dat Delta.dat Gamma4.dat
dmft_output_dir=dmft_patch$(date "+%d%m%y%H%M")
N_array=(  0    0   1   10  10   10    6   10   4 )
max_N=10                                # maximum amount of iterations
J_Array=(0.4 0.55 0.6 0.65 0.7 0.75 0.80 0.85 1.0 )
M=9                                     # amount of cycles
acc_array=(40000000 40000000 50000000 55000000 60000000 60000000  60000000 60000000 
60000000 80000000)
#once only input
../change calculate_Gamma4 0 
#Gamma4 is needed for a Dual Fermion calculation so we put it to 1 there

../change DMFT_sigma 1 
#In a Dual Fermion calculation we do not use a DMFT sigma so this one is 0 there
../change U 5.0
../change beta 10.0
mkdir -p $dmft_output_dir
for((m=0; m<M; m++))
do
	#input that changes between cycles
#	rm Delta.dat
	j=${J_Array[$m]}
	../change J $j
	N=${N_array[$m]}
	
	#copy Delta's to incorporate previous calculations
	cp latest_deltas/Delta.dat.$m Delta.dat 
	s=$(echo "$max_N-$N" | bc) #starting value for next loop
	for((i=s; i<max_N; i++))
	do
		echo "dmft iteration $i+1 of $max_N of run $m+1 of $M"
		echo "m=$m signifies J=$j"
		
		../change read_delta 2
		a=${acc_array[$i]}
		../change maximum_MC_steps $a
		cp _input.dat $dmft_output_dir/_input_run$m.$i
		./kondo22>$dmft_output_dir/log_ct-qmc$m.txt.$i
		cp Gw_complete.dat $dmft_output_dir/Gw_complete$m.dat.$i
		#Changing format of the output for compatibility
		cat Gw_complete.dat | awk '{print $1, ($2-($2-$8)/2), $3, 0.0, 0.0, 
			0.0, 0.0, ($8+($2-$8)/2), $9}'>Gw_cond.dat	
		cp Gw_cond.dat Gw_complete.dat
		rm Delta.dat
		./dual21>$dmft_output_dir/log_delta$m.txt.$i
		#Changing format of the output for compatibility
		cat Delta.dat | awk '{print  $1, $2, 0.0, 0.0, 0.0, 0.0, $7, $8, 
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}'>Delta22.dat
		cp Delta.dat Delta_old.dat
		cp Delta22.dat Delta.dat
		cp Delta.dat $dmft_output_dir/Delta$m.dat.$i
	done
#	cp Gamma4.dat $dmft_output_dir/Gamma4.dat.$m	#Gamma is copied in the Dual Fermion case
done
