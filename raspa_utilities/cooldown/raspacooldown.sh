#!/bin/bash
# This is a BASH code to sistematically run a cota cool-down single molecule calculation
#Daniele Ongari 2Dic2015, 8Jun2016, 9Sep2017

run_raspa='/home/daniele/Programs/raspa-2.0/bin/simulate'

if [ $# -lt 3 ]
  then
    echo "Not enough arguments: frame, molec, chargemethod required!"
    exit 1
fi

################# INPUT info ################################ CHECK!!!!
temp_start=300
temp_step=20
temp_finish=1 #this can not be 0 or raspa would crash!
ncycles=1000

ucx=1
ucy=1
ucz=1
frame=$1
molec=$2
chargemethod=$3 #Ewald or None

Nmolecules=1

#input needed in this folder#####################################: molecule.def, framework.cif, force_field.def + force_field_mixing_rules.def + pseudo_atoms.def. 

#############################################################
#for the other inputs, check the cat>simulation.input file below
#remember that Movies=yes @ last run

#Start removing old files
rm -rf 1_first* 2_inter* 3_last*
rm -f potential*

##################################### Run first simulation

echo " -------------------------------------------------- First step @ ${temp_start}K" 

cat > simulation.input <<EOF
SimulationType                MonteCarlo
NumberOfCycles                ${ncycles}
NumberOfInitializationCycles  0
PrintEvery                    $((ncycles/10))

RestartFile                   no
RestartStyle                  Raspa

CutOff                        13.00 

Forcefield                    Local
UseChargesFromCIFFile         yes

Charge                        ${chargemethod}
EwaldPrecision                1e-6

Framework 0
FrameworkName                 ${frame}
InputFileType                 cif
UnitCells                     1 1 1
HeliumVoidFraction            0.00
Movies                        no 
WriteMoviesEvery              ${ncycles}

ExternalTemperature           $temp_start

Component 0 MoleculeName             ${molec}
            MoleculeDefinition       Local
            TranslationProbability   1.0
            RotationProbability      1.0
            ReinsertionProbability   2.0
            SwapProbability          0.0
            CreateNumberOfMolecules  ${Nmolecules}
EOF

          newfolder="1_first_${temp_start}K"
mkdir ./${newfolder}
cd    ./${newfolder}
cp ../simulation.input ../*def ../*cif .
$run_raspa
rm *def *cif


##################################### Run intermediate cooling

temp=$((temp_start-temp_step))
i=0

while [  $temp -gt 0  ]              
do

i=$((i+1))
ii=$(printf %02d $i)

echo " -------------------------------------------------- Step number $ii @ ${temp}K" 
           newfolder="2_inter${ii}_${temp}K"
mkdir ../${newfolder}
cp simulation.input ../${newfolder}/
mkdir ../${newfolder}/RestartInitial
mkdir ../${newfolder}/RestartInitial/System_0
                                                                 restart_name="restart_${frame}_${ucx}.${ucy}.${ucz}_${temp}.000000_0"
cp Restart/System_0/*  ../${newfolder}/RestartInitial/System_0/${restart_name}

cd ../${newfolder}/

grep -rl " RestartFile " simulation.input                   | xargs sed -i "s/^.*RestartFile.*$/RestartFile                      yes/g" simulation.input
grep -rl " ExternalTemperature " simulation.input           | xargs sed -i "s/^.*ExternalTemperature.*$/ExternalTemperature    ${temp}/g" simulation.input
grep -rl " CreateNumberOfMolecules " simulation.input       | xargs sed -i "s/^.*CreateNumberOfMolecules.*$/CreateNumberOfMolecules          0/g" simulation.input
cp ../*def ../*cif .
$run_raspa
rm *def *cif

temp=$((temp-temp_step))
done

##################################### Run final with Movie

temp=$temp_finish && echo " -------------------------------------------------- Last step @ ${temp}K" 
           newfolder="3_last_${temp_finish}K"
mkdir ../${newfolder}
cp simulation.input ../${newfolder}/
mkdir ../${newfolder}/RestartInitial
mkdir ../${newfolder}/RestartInitial/System_0
                                                                 restart_name="restart_${frame}_${ucx}.${ucy}.${ucz}_${temp}.000000_0"
cp Restart/System_0/*  ../${newfolder}/RestartInitial/System_0/${restart_name}

cd ../${newfolder}/

grep -rl " RestartFile " simulation.input                   | xargs sed -i "s/^.*RestartFile.*$/RestartFile                      yes/g" simulation.input
grep -rl " ExternalTemperature " simulation.input           | xargs sed -i "s/^.*ExternalTemperature.*$/ExternalTemperature    ${temp}/g" simulation.input
grep -rl " Movies " simulation.input                        | xargs sed -i "s/^.*Movies .*$/Movies                           yes/g" simulation.input
cp ../*def ../*cif .
$run_raspa
rm *def *cif

cd ..
################## ANALYSIS

grep  'Current total potential energy:  ' */Output/*/*  | cut -d' ' --complement -f1 | sed -r 's/.{30}//' | sed -r 's/.{34}$//' >> pot_tot.txt
egrep 'Current Host-Host energy:         ' */Output/*/* | cut -d' ' --complement -f1 | sed -r 's/.{25}//' | sed -r 's/.{34}$//'  >> pot_hh.txt
egrep 'Current Host-Adsorbate energy:    ' */Output/*/* | cut -d' ' --complement -f1 | sed -r 's/.{25}//' | sed -r 's/.{34}$//'  >> pot_ha.txt
echo "#Total_pot(K) Host-host_pot(K) Host-Ads_pot(K)"   
paste pot_tot.txt pot_hh.txt pot_ha.txt | column -s $'\t' -t >> potential.txt
rm pot_*

cat > potential.gnuplot <<EOF
set term png 
set output 'potential.png'
set xlabel 'Cooling Steps'
set ylabel 'Potential Energy (kJ/mol)' 
plot 'potential.txt' u (\$1/120) w l lw 3 title "Total", 'potential.txt' u (\$2/120) w l title "Host-Host", 'potential.txt' u (\$3/120) w l title "Host-Adsorbate"
EOF

gnuplot potential.gnuplot
gnome-open potential.png 

echo "----------FINAL ENERGY (K)------------"
tail potential.txt

echo 
echo "---- print movie.vmd to visualize it with vmd ---"







