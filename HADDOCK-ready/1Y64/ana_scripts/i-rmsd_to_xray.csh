#!/bin/csh
#
source ~abonvin/haddock2.4/haddock_configure.csh
set WDIR=/Users/abonvin/Desktop/benchmark5/docking-ready/1Y64/ana_scripts
set refe=$WDIR/target.pdb
set izone=$WDIR/target.izone
set atoms='CA,C,N,O'
#
# Define the location of profit
#
if ( `printenv |grep PROFIT | wc -l` == 0) then
  set found=`which profit |wc -l`
  if ($found == 0) then
     echo 'PROFIT environment variable not defined'
     echo '==> no rmsd calculations '
  else
     setenv PROFIT `which profit`
  endif
endif

cat /dev/null >rmsd-interface_xray.disp

foreach i (`cat file.nam`)
  if ($i:e == "gz") then
    gzip -dc $i > $i:t:r.tmp2
  else
    \cp $i $i:t:r.tmp2
  endif
  $HADDOCKTOOLS/pdb_segid-to-chain $i:t:r.tmp2 |sed -e 's/BB/CA/' >$i:t:r.tmp1
  echo $i >>rmsd-interface_xray.disp
  $PROFIT <<_Eod_ |grep RMS |tail -1 >>rmsd-interface_xray.disp
    refe $refe
    mobi $i:t:r.tmp1
    `cat $izone`
    atom $atoms
    fit
    quit
_Eod_
\rm $i:t:r.tmp1
\rm $i:t:r.tmp2
end
echo "#struc i-RMSD" >i-RMSD.dat
awk '{if ($1 == "RMS:") {printf "%8.3f ",$2} else {printf "\n %s ",$1}}' rmsd-interface_xray.disp |grep pdb |awk '{print $1,$2}' >> i-RMSD.dat
head -1 i-RMSD.dat >i-RMSD-sorted.dat
grep pdb i-RMSD.dat |sort -nk2  >> i-RMSD-sorted.dat
\rm rmsd-interface_xray.disp
