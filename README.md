# Project_psyche
Scripts for project psyche

# Democon scripts:

calc_4d_pi.py

Run time: ~15 minutes


# Histex

while read species; do
mkdir -p $species/autosomes
Histex -G ${species}_Table | ~/apps/GENESCOPE.FK/GeneScopeFK.R -o $species/autosomes -k 31 > $species/autosomes/${species}_autosomes.log
done < ../Species_list.txt
