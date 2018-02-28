#!/bin/bash

#example: axsf_manage_crystal.sh xxxx.axsf "-x 2 -y 3"
#makes the folder: axsf_xxxx

baseinput=$(basename "$1")
title="${baseinput%.*}"

animsteps=$(sed '1q;d' $1 | awk '{print $2}')
atoms=$(sed '8q;d' $1 | awk '{print $1}')
options=$2

echo title=$title
echo animsteps=$animsteps
echo atoms=$atoms
echo options=$options

mkdir axsf_$title
cp $1 axsf_$title
cd axsf_$title


for a in $(seq 1 $animsteps); do 
	echo ANIMSTEPS 1 				> ${title}_animstep_${a}.axsf
	head -8 $1 | tail -7 				>> ${title}_animstep_${a}.axsf
	head -$((6+(2+atoms)*a)) $1 | tail -${atoms}	>> ${title}_animstep_${a}.axsf

	manage_crystal.py ${title}_animstep_${a}.axsf -o ${title}_animstep_${a}_managed.axsf $options
	atoms_managed=$(sed '8q;d' ${title}_animstep_${a}_managed.axsf | awk '{print $1}')

	if [ $a -eq 1 ]; then
		echo ANIMSTEPS $animsteps 					> ${title}_managed.axsf
		head -6 ${title}_animstep_${a}_managed.axsf | tail -5		>> ${title}_managed.axsf 
		fi	

	echo PRIMCOORD   1 							>> ${title}_managed.axsf
	echo ${atoms_managed} 1							>> ${title}_managed.axsf
	tail -${atoms_managed} ${title}_animstep_${a}_managed.axsf		>> ${title}_managed.axsf

	done


