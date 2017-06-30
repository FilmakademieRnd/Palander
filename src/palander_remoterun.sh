#!/bin/bash

# This source file contains the Palander remote running script, to be used in
# conjunction with the Blender add-on palander_addon.py, the simulation engine
# palander_engine.cpp and the script base palander_scriptbase. All files are
# available at
# <https://github.com/FilmakademieRnD/Palander>
#
# Copyright (c) 2017 Animationsinstitut, Filmakademie Baden-Wuerttemberg
# <http://research.animationsinstitut.de/>
#
# This project has been realized in the scope of the Media Solution Center
# project, funded by the state of Baden-Wuerttemberg, Germany:
# <http://msc-bw.de/en>
#
# While it is meant to work with the open-source programs Palabos and Blender,
# this file is not dependent on either program. Therefore, this file is released
# under the Beerware License, Revision 42:
#
# "THE BEER-WARE LICENSE" (Revision 42):
# <ari.harjunmaa@filmakademie.de> wrote this file. As long as you retain this notice
# you can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return.  - Ari Harjunmaa
#

echo --------------------------------------
echo HLRS: Building submission directory...
if ! cd ~/palander 2>/dev/null; then
	echo HLRS: Remote directory missing! Aborting...
	exit 1
fi
id=$(date +%d%m%y_%H%M%S)
mkdir run_$id
cd run_$id
mkdir tmp
cp ../palander_engine .
mv ../object*.stl ../toPalabos.xml .

mask=$(awk '/maskFile/ {print($2)}' toPalabos.xml)
nod=$(awk '/nodes/ {print($2)}' toPalabos.xml)
hrs=$(awk '/runhrs/ {print($2)}' toPalabos.xml)
mins=$(awk '/runmins/ {print($2)}' toPalabos.xml)
runtime=$(printf "%02d:%02d:00" $hrs $mins)

echo HLRS: Using $nod nodes for $hrs hours and $mins minutes.
awk '{if(/BASE_NODES/) {str = $0; sub("BASE_NODES",'$nod',str); print(str); cpus = 24 * '$nod'}
      else if(/BASE_CPUS/) {str = $0; sub("BASE_CPUS",cpus,str); print(str)}
      else if(/BASE_TIME/) {str = $0; sub("BASE_TIME","'$runtime'",str); print(str)}
      else if(/BASE_QUEUE/) {if('$hrs' == 0 && '$mins' <= 25 && '$nod' <= 384) {str = $0; sub("BASE_QUEUE","test",str); print(str)}}
      else {print}}' ../palander_scriptbase > ./script

if [ $mask ]; then
	if [ -f ../$mask ]; then
		cp ../$mask .
	else
		echo "cp $mask \$PBS_O_WORKDIR/.." >> ./script
	fi
fi

echo HLRS: Submitting to batch queue...
qsub script > /dev/null 2>&1
echo HLRS: Done!
echo --------------------------------------
