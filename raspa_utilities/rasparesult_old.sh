#!/bin/bash
#Extracts data from computed simulations

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
    heat_desorption_K=$(grep "Heat of desorption:" -A10 $1 | grep "Average" | awk '{print $2}')
    heat_desorption_K_error=$(grep "Heat of desorption:" -A10 $1 | grep "Average" | awk '{print $4}')
    heat_desorption_kJ=$(grep "Heat of desorption:" -A10 $1 | grep "KJ" | awk '{print $1}')
    heat_desorption_kJ_error=$(grep "Heat of desorption:" -A10 $1 | grep "KJ" | awk '{print $3}')
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

        fc=$(grep "Component [0-9] \[$comp\] (Adsorbate molecule)" -A50 $1 |\
             grep "Fugacity coefficient" | awk '{print $3}')

        pf=$(grep "Component [0-9] \[$comp\] (Adsorbate molecule)" -A50 $1 |\
                grep "Partial fugacity" | awk '{print $3}')

        avg1=$(grep "Number of molecules:" -A100 $1 | \
               grep "Component [0-9] \[$comp\]" -A27  | \
               grep "Average loading absolute \[molecules" |\
               awk '{print $6}')

        avg2=$(grep "Number of molecules:" -A100 $1 | \
               grep "Component [0-9] \[$comp\]" -A27  | \
               grep "Average loading absolute \[molecules" |\
               awk '{print $8}')

        avg3=$(grep "Number of molecules:" -A100 $1 | \
               grep "Component [0-9] \[$comp\]" -A27  | \
               grep "Average loading absolute \[mol/kg" |\
               awk '{print $6}')

        avg4=$(grep "Number of molecules:" -A100 $1 | \
               grep "Component [0-9] \[$comp\]" -A27  | \
               grep "Average loading absolute \[mol/kg" |\
               awk '{print $8}')

        avg5=$(grep "Number of molecules:" -A100 $1 | \
               grep "Component [0-9] \[$comp\]" -A27  | \
               grep "Average loading absolute \[milligram" |\
               awk '{print $6}')

        avg6=$(grep "Number of molecules:" -A100 $1 | \
               grep "Component [0-9] \[$comp\]" -A27  | \
               grep "Average loading absolute \[milligram" |\
               awk '{print $8}')



        ave1=$(grep "Number of molecules:" -A100 $1 | \
               grep "Component [0-9] \[$comp\]" -A27  | \
               grep "Average loading excess \[cm^3 (STP)/gr" |\
               awk '{print $7}')
        
        ave2=$(grep "Number of molecules:" -A100 $1 | \
               grep "Component [0-9] \[$comp\]" -A27  | \
               grep "Average loading excess \[cm^3 (STP)/gr" |\
               awk '{print $9}')

        ave3=$(grep "Number of molecules:" -A100 $1 | \
               grep "Component [0-9] \[$comp\]" -A27  | \
               grep "Average loading excess \[cm^3 (STP)/cm^3" |\
               awk '{print $7}')

        ave4=$(grep "Number of molecules:" -A100 $1 | \
               grep "Component [0-9] \[$comp\]" -A27  | \
               grep "Average loading excess \[cm^3 (STP)/cm^3" |\
               awk '{print $9}')




#        ave1=$(grep "Number of molecules:" -A100 $1 | \
#               grep "Component [0-9] \[$comp\]" -A27  | \
#               grep "Average loading excess \[molecules" |\
#               awk '{print $6}')
#        
#        ave2=$(grep "Number of molecules:" -A100 $1 | \
#               grep "Component [0-9] \[$comp\]" -A27  | \
#               grep "Average loading excess \[molecules" |\
#               awk '{print $8}')
#
#        ave3=$(grep "Number of molecules:" -A100 $1 | \
#               grep "Component [0-9] \[$comp\]" -A27  | \
#               grep "Average loading excess \[mol/kg" |\
#               awk '{print $6}')
#
#       ave4=$(grep "Number of molecules:" -A100 $1 | \
#                grep "Component [0-9] \[$comp\]" -A27  | \
#               grep "Average loading excess \[mol/kg" |\
#               awk '{print $8}')

        ave5=$(grep "Number of molecules:" -A100 $1 | \
               grep "Component [0-9] \[$comp\]" -A27  | \
               grep "Average loading excess \[milligram" |\
               awk '{print $6}')

        ave6=$(grep  "Number of molecules:" -A100 $1 | \
               grep "Component [0-9] \[$comp\]" -A27  | \
               grep "Average loading excess \[milligram" |\
               awk '{print $8}')

        echo ${framework:-"not_ready"}   ${temp:-"not_ready"}   ${comp:-"not_ready"}   ${pp:-"not_ready"}   ${pf:-"not_ready"}  ${fc:-"not_ready"}   ${avg1:-"not_ready"}   ${avg2:-"not_ready"}   ${avg3:-"not_ready"}   ${avg4:-"not_ready"}   ${avg5:-"not_ready"}   ${avg6:-"not_ready"}   ${ave1:-"not_ready"}   ${ave2:-"not_ready"}   ${ave3:-"not_ready"}   ${ave4:-"not_ready"}   ${ave5:-"not_ready"}   ${ave6:-"not_ready"}   ${heat_desorption_K:-"not_ready"}   ${heat_desorption_K_error:-"not_ready"}   ${heat_desorption_kJ:-"not_ready"}   ${heat_desorption_kJ_error:-"not_ready"}>> tmp
     
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
#    echo "#framework  temp  comp  partial_p  partial_f  fugacity_coeff  avg_load[molec/unit_cell]  avg_err[molec/unit_cell]  avg_load[mol/kg]  avg_err[mol/kg] avg_load[mg/g] avg_err[mg/g] exs_load[molec/unit_cell] exs_err[molec/unit_cell]  exs_load[mol/kg] exs_err[mol/kg] exs_load[mg/g] exc_err[mg/g] heat_of_des[K] heat_of_des_err[K] heat_of_des[kJ/mol] heat_of_des_err[kJ/mol]";
    cat $new_data_directory/$1 | { 
     echo "#framework  temp  comp  partial_p  partial_f  fugacity_coeff  avg_load[molec/unit_cell]  avg_err[molec/unit_cell]  avg_load[mol/kg]  avg_err[mol/kg] avg_load[mg/g] avg_err[mg/g] avg_load[cm^3(STP)/g] avg_err[cm^3(STP)/g]  avg_load[cm^3(STP)/cm^3] avg_err[cm^3(STP)/cm^3] exs_load[mg/g] exc_err[mg/g] heat_of_des[K] heat_of_des_err[K] heat_of_des[kJ/mol] heat_of_des_err[kJ/mol]";
        sort -k 3,3 -k 4,4n
    } | column -t  > $new_data_directory/$filename
    rm $new_data_directory/$1 
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
