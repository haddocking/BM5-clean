#!/bin/csh
#
set target=pdb
echo '######################################################'
echo '  '
echo '  HADDOCK scoring water'
echo '______________________________________________________'
foreach i ($argv)
  set clufi=$i/structures/it1/water/clusters_haddock-sorted.stat_best4
  echo '===========' $i '============='
  echo ' '
  echo '===> Cluster stats water <===='
  if (-e $clufi) then
    awk '{if ($6-$7<=4) {print "  ",NR-1,$1,$6,$7,$2,$3}}' $clufi | sed -e 's/file\.nam_//'| grep -v "#" >ttmp1
    cat /dev/null >ttmp2
    foreach j (`grep clust ttmp1 |awk '{print $2}'`)
      if (-e $i/structures/it1/water/clusters_best4.fnat) then
        grep $j"\ " $i/structures/it1/water/clusters_best4.fnat >>ttmp2
      endif
    end
    if (`wc -l ttmp1 |awk '{print $1}'` > 0) then
      if (-e $i/structures/it1/water/clusters_best4.fnat) then
        paste ttmp1 ttmp2 |awk '{printf "water acceptable clusters: Rank %2d %-8s i-RMSD %5.2f +/- %5.2f  Fnat= %5.2f +/- %5.2f HADDOCKscore= %5.1f +/- %5.1f\n",$1,$2,$3,$4,$8,$9,$5,$6}'
      else
        cat ttmp1 |awk '{printf "water acceptable clusters: Rank %2d %-8s i-RMSD %5.2f +/- %5.2f - - - - - - - -\n",$1,$2,$3,$4}'
      endif
    else
      awk '{print "  ",NR-1,$1,$6,$7,$2,$3}' $clufi | sed -e 's/file\.nam_//'| grep -v "#" >ttmp1
      cat /dev/null >ttmp2
      foreach j (`grep clust ttmp1 |awk '{print $2}'`)
        if (-e $i/structures/it1/water/clusters_best4.fnat) then
          grep $j"\ " $i/structures/it1/water/clusters_best4.fnat >>ttmp2
        endif
      end
      if (`wc -l ttmp1 |awk '{print $1}'` > 0) then
        if (-e $i/structures/it1/water/clusters_best4.fnat) then
          paste ttmp1 ttmp2 |awk '{printf "water best cluster: Rank %2d %-8s i-RMSD %5.2f +/- %5.2f  Fnat= %5.2f +/- %5.2f HADDOCKscore= %5.1f +/- %5.1f\n",$1,$2,$3,$4,$8,$9,$5,$6}' | sort -nk8 |head -1
        else
          cat ttmp1 |awk '{printf "water best cluster: Rank %2d %-8s i-RMSD %5.2f +/- %5.2f - - - - - - - -\n",$1,$2,$3,$4}' | sort -nk8 |head -1
        endif 
      else
        echo "no water clusters generated - - - - - - - - - - - - - -"
      endif
    endif
  endif
  echo ' '
  set clufi=$i/structures/it1/clusters_haddock-sorted.stat_best4
  echo '===> Cluster stats it1 <===='
  if (-e $clufi) then
    awk '{if ($6-$7<=4) {print "  ",NR-1,$1,$6,$7,$2,$3}}' $clufi | sed -e 's/file\.nam_//'| grep -v "#" >ttmp1
    cat /dev/null >ttmp2
    foreach j (`grep clust ttmp1 |awk '{print $2}'`)
      if (-e $i/structures/it1/clusters_best4.fnat) then
        grep $j"\ " $i/structures/it1/clusters_best4.fnat >>ttmp2
      endif
    end
    if (`wc -l ttmp1 |awk '{print $1}'` > 0) then
      if (-e $i/structures/it1/clusters_best4.fnat) then
        paste ttmp1 ttmp2 |awk '{printf "it1   acceptable clusters: Rank %2d %-8s i-RMSD %5.2f +/- %5.2f  Fnat= %5.2f +/- %5.2f HADDOCKscore= %5.1f +/- %5.1f\n",$1,$2,$3,$4,$8,$9,$5,$6}'
      else
        cat ttmp1 |awk '{printf "it1   acceptable clusters: Rank %2d %-8s i-RMSD %5.2f +/- %5.2f - - - - - - - -\n",$1,$2,$3,$4}'
      endif
    else
      awk '{print "  ",NR-1,$1,$6,$7,$2,$3}' $clufi | sed -e 's/file\.nam_//'| grep -v "#" >ttmp1
      cat /dev/null >ttmp2
      foreach j (`grep clust ttmp1 |awk '{print $2}'`)
        if (-e $i/structures/it1/clusters_best4.fnat) then
          grep $j"\ " $i/structures/it1/clusters_best4.fnat >>ttmp2
        endif
      end
      if (`wc -l ttmp1 |awk '{print $1}'` > 0) then
        if (-e $i/structures/it1/clusters_best4.fnat) then
          paste ttmp1 ttmp2 |awk '{printf "it1   best cluster: Rank %2d %-8s i-RMSD %5.2f +/- %5.2f  Fnat= %5.2f +/- %5.2f HADDOCKscore= %5.1f +/- %5.1f\n",$1,$2,$3,$4,$8,$9,$5,$6}' | sort -nk8 |head -1
        else
          cat ttmp1 |awk '{printf "it1   best cluster: Rank %2d %-8s i-RMSD %5.2f +/- %5.2f - - - - - - - -\n",$1,$2,$3,$4}' | sort -nk8 |head -1
        endif 
      else
        echo "no it1   clusters generated - - - - - - - - - - - - - -"
      endif
    endif
  endif
  #
  # check i-RMSDs at water
  #
  echo ' '
  echo '===> water single structure i-RMSD stats <===='
  set accept=0
  set rmsfile=$i/structures/it1/water/i-RMSD.dat
  set rmsbestfile=$i/structures/it1/water/i-RMSD-sorted.dat
  set fnatfile=$i/structures/it1/water/file.nam_fnat
  grep $target $rmsfile | awk '{if ($2<=2) {print NR,$1,$2}}' | head -1 >ttmp1
  if (`grep $target ttmp1 | wc -l |awk '{print $1}'` > 0) then
    grep `awk '{print $2}' ttmp1` $fnatfile >ttmp2
    paste ttmp1 ttmp2 |awk '{printf  "water first medium <2A          : Rank %4d i-RMSD: %5.2f  Fnat: %5.2f\n",$1,$3,$5}'
    set accept=1
  else
    echo "no water first medium <2A - - - - - - -"
  endif
  grep $target $rmsfile | awk '{if ($2<=4) {print NR,$1,$2}}' | head -1 >ttmp1
  if (`grep $target ttmp1 | wc -l |awk '{print $1}'` > 0) then
    grep `awk '{print $2}' ttmp1` $fnatfile >ttmp2
    paste ttmp1 ttmp2 |awk '{printf  "water first acceptable <4A      : Rank %4d i-RMSD: %5.2f  Fnat: %5.2f\n",$1,$3,$5}'
    set accept=1
  else
    echo "no water first acceptable <4A - - - - - - -"
  endif
  grep $target $rmsfile | awk '{if ($2<=4) {print NR,$1,$2}}' | head -1 >ttmp1
  if (`grep $target ttmp1 | wc -l |awk '{print $1}'` > 0) then
    grep $target $rmsbestfile | head -1 | awk '{if ($2<=4) {print $1}}' | head -1 >ttmp0
    grep -n `cat ttmp0` $rmsfile |sed -e 's/:/\ /g' >ttmp1
    grep `awk '{print $2}' ttmp1` $fnatfile >ttmp2
    paste ttmp1 ttmp2 |awk '{printf  "water best acceptable or higher : Rank %4d i-RMSD: %5.2f  Fnat: %5.2f\n",$1-1,$3,$5}'
    set accept=1
  else
    echo "no water best acceptable or higher - - - - - - -"
  endif
  #
  # check i-RMSDs at it1
  #
  echo ' '
  echo '===> it1 single structure i-RMSD stats <===='
  set rmsfile=$i/structures/it1/i-RMSD.dat
  set rmsbestfile=$i/structures/it1/i-RMSD-sorted.dat
  set fnatfile=$i/structures/it1/file.nam_fnat
  grep $target $rmsfile | awk '{if ($2<=2) {print NR,$1,$2}}' | head -1 >ttmp1
  if (`grep $target ttmp1 | wc -l |awk '{print $1}'` > 0) then
    grep `awk '{print $2}' ttmp1` $fnatfile >ttmp2
    paste ttmp1 ttmp2 |awk '{printf  "it1 first medium <2A            : Rank %4d i-RMSD: %5.2f  Fnat: %5.2f\n",$1,$3,$5}'
    set accept=1
  else
    echo "no it1 acceptable <2A - - - - - - -"
  endif
  grep $target $rmsfile | awk '{if ($2<=4) {print NR,$1,$2}}' | head -1 >ttmp1
  if (`grep $target ttmp1 | wc -l |awk '{print $1}'` > 0) then
    grep `awk '{print $2}' ttmp1` $fnatfile >ttmp2
    paste ttmp1 ttmp2 |awk '{printf  "it1 first acceptable <4A        : Rank %4d i-RMSD: %5.2f  Fnat: %5.2f\n",$1,$3,$5}'
    set accept=1
  else
    echo "no it1 first acceptable <4A - - - - - - -"
  endif
  grep $target $rmsfile | awk '{if ($2<=4) {print NR,$1,$2}}' | head -1 >ttmp1
  if (`grep $target ttmp1 | wc -l |awk '{print $1}'` > 0) then
    grep $target $rmsbestfile | head -1 | awk '{if ($2<=4) {print $1}}' | head -1 >ttmp0
    grep -n `cat ttmp0` $rmsfile |sed -e 's/:/\ /g' >ttmp1
    grep `awk '{print $2}' ttmp1` $fnatfile >ttmp2
    paste ttmp1 ttmp2 |awk '{printf  "it1 best acceptable or higher   : Rank %4d i-RMSD: %5.2f  Fnat: %5.2f\n",$1-1,$3,$5}'
    set accept=1
  else
    echo "no it1 best acceptable or higher - - - - - - -"
  endif
  #
  # check i-RMSDs at it0
  #
  echo ' '
  echo '===> it0 single structure i-RMSD stats <===='
  set rmsfile=$i/structures/it0/i-RMSD.dat
  set rmsbestfile=$i/structures/it0/i-RMSD-sorted.dat
  set fnatfile=$i/structures/it0/file.nam_fnat
  grep $target $rmsfile | awk '{if ($2<=2) {print NR,$1,$2}}' | head -1 >ttmp1
  if (`grep $target ttmp1 | wc -l |awk '{print $1}'` > 0) then
    grep `awk '{print $2}' ttmp1` $fnatfile >ttmp2
    paste ttmp1 ttmp2 |awk '{printf  "it0 first medium <2A            : Rank %4d i-RMSD: %5.2f  Fnat: %5.2f\n",$1,$3,$5}'
    set accept=1
  else
    echo "no it0 acceptable <2A - - - - - - -"
  endif
  grep $target $rmsfile | awk '{if ($2<=4) {print NR,$1,$2}}' | head -1 >ttmp1
  if (`grep $target ttmp1 | wc -l |awk '{print $1}'` > 0) then
    grep `awk '{print $2}' ttmp1` $fnatfile >ttmp2
    paste ttmp1 ttmp2 |awk '{printf  "it0 first acceptable <4A        : Rank %4d i-RMSD: %5.2f  Fnat: %5.2f\n",$1,$3,$5}'
    set accept=1
  else
    echo "no it0 first acceptable <4A - - - - - - -"
  endif
  grep $target $rmsfile | awk '{if ($2<=4) {print NR,$1,$2}}' | head -1 >ttmp1
  if (`grep $target ttmp1 | wc -l |awk '{print $1}'` > 0) then
    grep $target $rmsbestfile | head -1 | awk '{if ($2<=4) {print $1}}' | head -1 >ttmp0
    grep -n `cat ttmp0` $rmsfile |sed -e 's/:/\ /g' >ttmp1
    grep `awk '{print $2}' ttmp1` $fnatfile >ttmp2
    paste ttmp1 ttmp2 |awk '{printf  "it0 best acceptable or higher   : Rank %4d i-RMSD: %5.2f  Fnat: %5.2f\n",$1-1,$3,$5}'
    set accept=1
  else
    echo "no it0 best acceptable or higher - - - - - - -"
  endif

  echo ' '
  echo '===> Overall number of acceptable or better models <==='
  echo 'it0: structures with i-RMSD<4A:' |awk '{printf "%s   %s %s %s %s ",$1,$2,$3,$4,$5}'
  grep pdb $i/structures/it0/i-RMSD.dat | awk '$2<=4' | wc -l
  echo 'it0: structures within best200 with i-RMSD<4A:' |awk '{printf "%s   %s %s %s %s %s ",$1,$2,$3,$4,$5,$6}'
  grep pdb $i/structures/it0/i-RMSD.dat | head -200 | awk '$2<=4' | wc -l
  echo 'it0: structures within best200 with i-RMSD<2A:' |awk '{printf "%s   %s %s %s %s %s ",$1,$2,$3,$4,$5,$6}'
  grep pdb $i/structures/it0/i-RMSD.dat | head -200 | awk '$2<=2' | wc -l
  echo 'it0: structures within best200 with i-RMSD<1A:' |awk '{printf "%s   %s %s %s %s %s ",$1,$2,$3,$4,$5,$6}'
  grep pdb $i/structures/it0/i-RMSD.dat | head -200 | awk '$2<=1' | wc -l
  echo 'it1: structures with i-RMSD<4A:' |awk '{printf "%s   %s %s %s %s ",$1,$2,$3,$4,$5}'
  grep pdb $i/structures/it1/i-RMSD.dat | awk '$2<=4' | wc -l
  echo 'it1: structures with i-RMSD<2A:' |awk '{printf "%s   %s %s %s %s ",$1,$2,$3,$4,$5}'
  grep pdb $i/structures/it1/i-RMSD.dat | awk '$2<=2' | wc -l
  echo 'it1: structures with i-RMSD<21:' |awk '{printf "%s   %s %s %s %s ",$1,$2,$3,$4,$5}'
  grep pdb $i/structures/it1/i-RMSD.dat | awk '$2<=1' | wc -l
  echo 'water: structures with i-RMSD<4A:' |awk '{printf "%s %s %s %s %s ",$1,$2,$3,$4,$5}'
  grep pdb $i/structures/it1/water/i-RMSD.dat | awk '$2<=4' | wc -l
  echo 'water: structures with i-RMSD<2A:' |awk '{printf "%s %s %s %s %s ",$1,$2,$3,$4,$5}'
  grep pdb $i/structures/it1/water/i-RMSD.dat | awk '$2<=2' | wc -l
  echo 'water: structures with i-RMSD<1A:' |awk '{printf "%s %s %s %s %s ",$1,$2,$3,$4,$5}'
  grep pdb $i/structures/it1/water/i-RMSD.dat | awk '$2<=1' | wc -l
end

\rm ttmp[012]
