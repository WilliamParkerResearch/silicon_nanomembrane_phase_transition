#!/bin/bash

PREFIX=Si.I4_1amd


nks_betasn=$(grep "&plot" $PREFIX.bands.freq| awk '{print $5}')
nbnd_betasn=$(grep "&plot" $PREFIX.bands.freq| awk '{print $3}'| tr -d ',')


awk '/^            /{n+=1}{print > "kpoint_betasn_"n""}' $PREFIX.freq

for n in `seq 1 1 $nks_betasn`
do
	sed '1d' kpoint_betasn_$n | paste -s > ordered_kpoint_betasn_$n
done

rm kpoint_betasn_*

paste_files_betasn=$(seq -f "ordered_kpoint_betasn_%g" $nks_betasn| tr '\n' ' ')

awk 'NF' $paste_files_betasn > frequencies_unform_betasn

rm ordered_kpoint_betasn_*

kpoints_betasn=$(seq 1 1 $nks_betasn | paste -s -d ',') 



grep "      " $PREFIX.freq | nl >> kpoint_list_betasn

path_kpoints_betasn=$(sed -n '/^\/$/ { :a; n; p; ba; }' $PREFIX.matdyn.in | sed 1d | awk '{print $1"-"$2"-"$3}' | sed "s/\./\\\./g" | sed 's/\-/.*/g')

for pattern in $path_kpoints_betasn
do
	grep "$pattern" kpoint_list_betasn >> path_kpoints_betasn
done

rm kpoint_list_betasn

path_kpoint_tik_betasn=$(sort path_kpoints_betasn | uniq | awk '{print $1}' | paste -s -d ',')
path_kpoint_labels_betasn=$(sed -n '/^\/$/ { :a; n; p; ba; }' $PREFIX.matdyn.in | sed 1d | awk '{print "\""$NF"\""}'| paste -s -d ',')



dos_frequency_betasn=$(awk '{print $1}' $PREFIX.dos | sed 1d | paste -s -d ',')
dos_betasn=$(awk '{print $2}' $PREFIX.dos | sed 1d | paste -s -d ',')


#------------------------Output to QEData.py------------------------------------
echo " " > QEData.py

echo "nks_betasn = $nks_betasn" >> QEData.py
echo "nbnd_betasn = $nbnd_betasn" >> QEData.py
echo "path_kpoint_tik_betasn = [$path_kpoint_tik_betasn]" >> QEData.py
echo "path_kpoint_lables_betasn = [$path_kpoint_labels_betasn]" >> QEData.py

echo "dos_frequencies_betasn = np.array([$dos_frequency_betasn])" >> QEData.py
echo "dos_betasn = np.array([$dos_betasn])" >> QEData.py

echo "kpoints_betasn = np.array([$kpoints_betasn])" >> QEData.py

for index in `seq 1 1 $nbnd_betasn`
do
        freq_betasn=$(awk "{print \$$index}" frequencies_unform_betasn | paste -s -d ',')
        echo "frequency_betasn_$index = np.array([$freq_betasn])" >> QEData.py
done
rm frequencies_unform_betasn
