#!/bin/csh
#
grep "cluster" $1 |grep water | awk '{print $5}'  |sed -e 's/\./\,/g'  >cluster-water-rank.txt
grep "cluster" $1 |grep water | awk '{print $8}'  |sed -e 's/\./\,/g'  >cluster-water-rmsd.txt
grep "cluster" $1 |grep water | awk '{print $10}' |sed -e 's/\./\,/g'  >cluster-water-rmsd-sd.txt
grep "cluster" $1 |grep water | awk '{print $16}' |sed -e 's/\./\,/g'  >cluster-water-score.txt
grep "cluster" $1 |grep water | awk '{print $18}' |sed -e 's/\./\,/g'  >cluster-water-score-sd.txt
grep "cluster" $1 |grep it1 | awk '{print $5}'  |sed -e 's/\./\,/g'  >cluster-it1-rank.txt
grep "cluster" $1 |grep it1 | awk '{print $8}'  |sed -e 's/\./\,/g'  >cluster-it1-rmsd.txt
grep "cluster" $1 |grep it1 | awk '{print $10}' |sed -e 's/\./\,/g'  >cluster-it1-rmsd-sd.txt
grep "cluster" $1 |grep it1 | awk '{print $16}' |sed -e 's/\./\,/g'  >cluster-it1-score.txt
grep "cluster" $1 |grep it1 | awk '{print $18}' |sed -e 's/\./\,/g'  >cluster-it1-score-sd.txt
grep -A4 "water single structure i-RMSD" $1 | grep "water first acceptable" | awk '{print $7}'  |sed -e 's/\./\,/g' >water-first-acc-rank.txt
grep -A4 "water single structure i-RMSD" $1 | grep "water first acceptable" | awk '{print $9}'  |sed -e 's/\./\,/g' > water-first-acc-rmsd.txt
grep -A4 "water single structure i-RMSD" $1 | grep "water best acceptable"  | awk '{print $8}'  |sed -e 's/\./\,/g' > water-best-acc-rank.txt
grep -A4 "water single structure i-RMSD" $1 | grep "water best acceptable"  | awk '{print $10}' |sed -e 's/\./\,/g' > water-best-acc-rmsd.txt
grep -A4 "it1 single structure i-RMSD" $1 | grep "it1 first acceptable" | awk '{print $7}'  |sed -e 's/\./\,/g' >it1-first-acc-rank.txt
grep -A4 "it1 single structure i-RMSD" $1 | grep "it1 first acceptable" | awk '{print $9}'  |sed -e 's/\./\,/g' > it1-first-acc-rmsd.txt
grep -A4 "it1 single structure i-RMSD" $1 | grep "it1 best acceptable"  | awk '{print $8}'  |sed -e 's/\./\,/g' > it1-best-acc-rank.txt
grep -A4 "it1 single structure i-RMSD" $1 | grep "it1 best acceptable"  | awk '{print $10}' |sed -e 's/\./\,/g' > it1-best-acc-rmsd.txt
grep -A4 "it0 single structure i-RMSD" $1 | grep "it0 first acceptable" | awk '{print $7}'  |sed -e 's/\./\,/g' >it0-first-acc-rank.txt
grep -A4 "it0 single structure i-RMSD" $1 | grep "it0 first acceptable" | awk '{print $9}'  |sed -e 's/\./\,/g' > it0-first-acc-rmsd.txt
grep -A4 "it0 single structure i-RMSD" $1 | grep "it0 best acceptable or higher" | awk '{print $8}' |sed -e 's/\./\,/g' > it0-best-acc-rank.txt
grep -A4 "it0 single structure i-RMSD" $1 | grep "it0 best acceptable or higher" | awk '{print $10}' |sed -e 's/\./\,/g' > it0-best-acc-rmsd.txt
grep "it0:   structures with i-RMSD<2A" $1 | awk '{print $NF}' |sed -e 's/\./\,/g' > it0.txt
grep "it0:   structures within best200 with i-RMSD<1A:" $1 | awk '{print $NF}' |sed -e 's/\./\,/g' > it0-top200-1.txt
grep "it0:   structures within best200 with i-RMSD<2A:" $1 | awk '{print $NF}' |sed -e 's/\./\,/g' > it0-top200-2.txt
grep "it1:   structures with i-RMSD<1A:" $1 | awk '{print $NF}' |sed -e 's/\./\,/g' > it1-1.txt
grep "it1:   structures with i-RMSD<2A:" $1 | awk '{print $NF}' |sed -e 's/\./\,/g' > it1-2.txt
grep "water: structures with i-RMSD<1A:" $1 | awk '{print $NF}' |sed -e 's/\./\,/g' > water-1.txt
grep "water: structures with i-RMSD<2A:" $1 | awk '{print $NF}' |sed -e 's/\./\,/g' > water-2.txt
paste cluster-water-rank.txt cluster-water-rmsd.txt cluster-water-rmsd-sd.txt cluster-water-score.txt cluster-it1-rank.txt cluster-it1-rmsd.txt cluster-it1-rmsd-sd.txt cluster-it1-score.txt >clusters.txt
\rm cluster-*.txt
paste water-best-acc-rank.txt water-best-acc-rmsd.txt water-first-acc-rank.txt water-first-acc-rmsd.txt >water-struc.txt
paste it1-best-acc-rank.txt it1-best-acc-rmsd.txt it1-first-acc-rank.txt it1-first-acc-rmsd.txt >it1-struc.txt
paste it0-best-acc-rank.txt it0-best-acc-rmsd.txt it0-first-acc-rank.txt it0-first-acc-rmsd.txt >it0-struc.txt
paste water-struc.txt it1-struc.txt it0-struc.txt >structures.txt
paste it0.txt it0-top200-2.txt it1-2.txt water-2.txt it0-top200-1.txt it1-1.txt water-1.txt >acceptable.txt
\rm water-* it1* it0*

