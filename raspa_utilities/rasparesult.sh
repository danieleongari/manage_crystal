#!/bin/bash
# Extracts data from computed simulations
# Taken from Efrem, modified on 12May2017 to exclude the excess quantities, and check the running statistics if the simulation hasn't finished yet
# Check in order: (1) finished statistics, (2) running average (avg.), (3) snapshot value in the "Init" run

new_data_directory='Values_from_results'

rm -rf $new_data_directory

function getdata() {

    echo " Reading file $1"
  
    unset comp_array
    unset mol_frac_array

    for index in 0 1 2 3 4 5 6 7 8 9
    do
        comp=$(grep -m 1 "Component $index" $1 | awk '{print $3}')
        comp="${comp#[}"
        comp="${comp%]}"
        mol_frac=$(grep -m 1 -A 20 "Component $index" $1 | grep "MolFraction"  | awk '{print $2}')
        mol_frac=${mol_frac:0:5}
        if [ -n "$comp" ];
        then
            comp_array[$index]=$comp
            mol_frac_array[$index]=$mol_frac
        fi
    done
    
    framework=$(grep "Framework name" $1 | awk '{print $3}')
    temp=$(grep "External temperature" $1 | awk '{print $3}')
    heat_des_K=$(grep "Heat of desorption:" -A10 $1 | grep "Average" | awk '{print $2}');       if [ -z "$heat_des_K" ];      then heat_des_K=0.0; fi 
    heat_des_K_err=$(grep "Heat of desorption:" -A10 $1 | grep "Average" | awk '{print $4}');   if [ -z "$heat_des_K_err" ];  then heat_des_K_err=0.0; fi
    heat_des_kJ=$(grep "Heat of desorption:" -A10 $1 | grep "KJ" | awk '{print $1}');           if [ -z "$heat_des_kJ" ];     then heat_des_kJ=0.0; fi
    heat_des_kJ_err=$(grep "Heat of desorption:" -A10 $1 | grep "KJ" | awk '{print $3}');       if [ -z "$heat_des_kJ_err" ]; then heat_des_kJ_err=0.0; fi
    extra1=0.0 #to maintain the same column number as before for the heat of desorption
    extra2=0.0
    ERROR_message=$(grep "ERROR" $1)
    Error_message=$(grep "Error" $1)
    WARNING_message=$(grep "WARNING" $1)
    Warning_message=$(grep "Warning" $1)
    if [ -n "$ERROR_message" -o -n "$Error_message" -o -n "$WARNING_message" -o -n "$Warning_message" ];
    then
        error_messages=$(($error_messages + 1))
    fi

    for comp in "${comp_array[@]}"
    do
       
        pp=$(grep "Component [0-9] \[$comp\] (Adsorbate molecule)" -A50 $1 |\
             grep "Partial pressure" | awk '{print $3}')
        pp=${pp::-13}


        pf=$(grep "Component [0-9] \[$comp\] (Adsorbate molecule)" -A50 $1 |\
                grep "Partial fugacity" | awk '{print $3}')
        pf=${pf::-13}

        fc=$(grep "Component [0-9] \[$comp\] (Adsorbate molecule)" -A50 $1 |\
             grep "Fugacity coefficient" | awk '{print $3}')


#Example
        #[molec/uc] absolute
        avg1=$(grep "Number of molecules:" -A100 $1 | \
               grep "Component [0-9] \[$comp\]" -A27  | \
               grep "Average loading absolute \[molecules" |\
               awk '{print $6}')
 #if not complete, look for the current average 
        if [ -z "$avg1" ]; then
        avg1=$(grep "Component [0-9] ($comp)" -A2  $1 | \
               grep "avg." |\
               tail -2 | head -1 |\
               awk '{print $5}') 
        fi         
 #if not equilibrated read the value of the last snapshot
        if [ -z "$avg1" ]; then
        avg1=$(grep "Component [0-9] ($comp)" -A2  $1 | \
               tail -2 | head -1 |\
               awk '{print $3}') 
        else #remove the last bracket of (avg. 4.3446)
        avg1=${avg1::-1}
        fi    
 #err[molec/uc] absolute
        avg1_err=$(grep "Number of molecules:" -A100 $1 | \
               grep "Component [0-9] \[$comp\]" -A27  | \
               grep "Average loading absolute \[molecules" |\
               awk '{print $8}')
 #print an error equal to zero
        if [ -z "$avg1_err" ]; then avg1_err=0; fi



        #[mmol/g] absolute
        avg2=$(grep "Number of molecules:" -A100 $1 | \
               grep "Component [0-9] \[$comp\]" -A27  | \
               grep "Average loading absolute \[mol/kg" |\
               awk '{print $6}')

        if [ -z "$avg2" ]; then
        avg2=$(grep "Component [0-9] ($comp)" -A2  $1 | \
               grep "avg." |\
               tail -2 | head -1 |\
               awk '{print $9}') 
        fi         
        if [ -z "$avg2" ]; then
        avg2=$(grep "Component [0-9] ($comp)" -A2  $1 | \
               tail -2 | head -1 |\
               awk '{print $5}') 
        else 
        avg2=${avg2::-1}
        fi    
        avg2_err=$(grep "Number of molecules:" -A100 $1 | \
               grep "Component [0-9] \[$comp\]" -A27  | \
               grep "Average loading absolute \[mol/kg" |\
               awk '{print $8}')
        if [ -z "$avg2_err" ]; then avg2_err=0; fi



        #[mg/g] absolute
        avg3=$(grep "Number of molecules:rasparesult_new.sh" -A100 $1 | \
               grep "Component [0-9] \[$comp\]" -A27  | \
               grep "Average loading absolute \[milligram" |\
               awk '{print $6}')
        if [ -z "$avg3" ]; then
        avg3=$(grep "Component [0-9] ($comp)" -A2  $1 | \
               grep "avg." |\
               tail -2 | head -1 |\
               awk '{print $9}') 
        fi         
        if [ -z "$avg3" ]; then
        avg3=$(grep "Component [0-9] ($comp)" -A2  $1 | \
               tail -2 | head -1 |\
               awk '{print $5}') 
        else 
        avg3=${avg3::-1}
        fi    
        avg3_err=$(grep "Number of molecules:" -A100 $1 | \
               grep "Component [0-9] \[$comp\]" -A27  | \
               grep "Average loading absolute \[milligram" |\
               awk '{print $8}')
        if [ -z "$avg3_err" ]; then avg3_err=0; fi

        #[cm3/g] absolute 
        avg4=$(grep "Number of molecules:" -A100 $1 | \
               grep "Component [0-9] \[$comp\]" -A27  | \
               grep "Average loading absolute \[cm^3 (STP)/gr" |\
               awk '{print $7}')
        if [ -z "$avg4" ]; then
        avg4=$(grep "Component [0-9] ($comp)" -A2  $1 | \
               grep "avg." |\
               tail -2 | head -1 |\
               awk '{print $9}') 
        fi         
        if [ -z "$avg4" ]; then
        avg4=$(grep "Component [0-9] ($comp)" -A2  $1 | \
               tail -2 | head -1 |\
               awk '{print $5}') 
        else 
        avg4=${avg4::-1}
        fi    
        avg4_err=$(grep "Number of molecules:" -A100 $1 | \
               grep "Component [0-9] \[$comp\]" -A27  | \
               grep "Average loading absolute \[cm^3 (STP)/gr" |\
               awk '{print $9}')
        if [ -z "$avg4_err" ]; then avg4_err=0; fi

        #[cm3/cm3] absolute
        avg5=$(grep "Number of molecules:" -A100 $1 | \
               grep "Component [0-9] \[$comp\]" -A27  | \
               grep "Average loading absolute \[cm^3 (STP)/cm^3" |\
               awk '{print $7}')
        if [ -z "$avg5" ]; then
        avg5=$(grep "Component [0-9] ($comp)" -A2  $1 | \
               grep "avg." |\
               tail -2 | head -1 |\
               awk '{print $9}') 
        fi         
        if [ -z "$avg5" ]; then
        avg5=$(grep "Component [0-9] ($comp)" -A2  $1 | \
               tail -2 | head -1 |\
               awk '{print $5}') 
        else 
        avg5=${avg5::-1}
        fi    
        avg5_err=$(grep "Number of molecules:" -A100 $1 | \
               grep "Component [0-9] \[$comp\]" -A27  | \
               grep "Average loading absolute \[cm^3 (STP)/cm^3" |\
               awk '{print $9}')
        if [ -z "$avg5_err" ]; then avg5_err=0; fi

        echo $framework $temp $comp $pp $pf $fc $avg1 $avg1_err $avg2 $avg2_err $avg3 $avg3_err $avg4 $avg4_err $avg5 $avg5_err $extra1 $extra2 $heat_des_K $heat_des_K_err $heat_des_kJ $heat_des_kJ_err >> tmp
     
        data_to_file $comp
        rm tmp
   done
}

function data_to_file() {
    temp_short=${temp:0:5}
    filename='adsiso_'$framework'_T'$temp_short'_numcomps'${#comp_array[@]}
    index=0
    for component in ${comp_array[@]}
    do
        mol_frac=${mol_frac_array[$index]:0:5}
        filename=$filename'_'$component'-'$mol_frac
        index=$(($index + 1))
    done
    filename=$filename'_'$1'_tmp.txt'
    cat tmp >> $new_data_directory/$filename
}

function sortdata() {
    filename=${1%_tmp.txt}'.txt'
    cat $new_data_directory/$1 | { 
     echo "#0framework  1temp  2comp  3partial_p  4partial_f  5fugacity_coeff  6avg_load[molec/UC]  7avg_err[molec/UC]  8avg_load[mol/kg]  9avg_err[mol/kg] 10avg_load[mg/g] 11avg_err[mg/g] 12avg_load[cm^3(STP)/g] 13avg_err[cm^3(STP)/g]  14avg_load[cm^3(STP)/cm^3] 15avg_err[cm^3(STP)/cm^3] 16extra1 17extra2 18heat_of_des[K] 19heat_of_des_err[K] 20heat_of_des[kJ/mol] 21heat_of_des_err[kJ/mol]";
        sort -k 3,3 -k 4,4n
    } | column -t  > $new_data_directory/$filename
    rm $new_data_directory/$1    #unsorted temporary
}


mkdir $new_data_directory

error_messages=0
for file in $(find ./Results/ -name "output_*.data")
do
    getdata $file
done

for file in $(ls $new_data_directory)
do
    sortdata $file
done

echo "Done"
echo "ERROR OR WARNING MESSAGES: "$error_messages

exit 0
