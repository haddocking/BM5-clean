#!/bin/csh
#
#
# Define the location of profit
#
if ( `printenv |grep PROFIT | wc -l` == 0) then
  set found=`which profit |wc -l`
  if ($found == 0) then
     echo 'PROFIT environment variable not defined'
     echo '==> no rmsd calculations '
     goto exit
  else
     setenv PROFIT `which profit`
  endif
endif

set atoms='CA,C,N,O'
cat /dev/null >rmsd-interface_xray.disp
cat /dev/null >rmsd-interface_r_xray.disp
cat /dev/null >rmsd-interface_l_xray.disp

foreach i ($argv)
  set WDIR=$i
  set refe=$WDIR/ana_scripts/target.pdb
  set unbound=$WDIR/ana_scripts/target-unbound.pdb
  set izone=$WDIR/ana_scripts/target.izone
  set izoneA=$WDIR/ana_scripts/target.izoneA
  set izoneB=$WDIR/ana_scripts/target.izoneB
  set strucr=${i}/${i}_r_u.pdb
  set strucl=${i}/${i}_l_u.pdb
  echo ${i}/ana_script/target-unbound.pdb >>rmsd-interface_xray.disp
  echo $strucr >>rmsd-interface_r_xray.disp
  echo $strucl >>rmsd-interface_l_xray.disp

  $PROFIT <<_Eod_ |grep RMS |tail -n1 >>rmsd-interface_xray.disp
    refe $refe
    mobi $unbound
    `cat $izone`
    atom $atoms
    fit
    quit
_Eod_

  $PROFIT <<_Eod_ |grep RMS |tail -n1 >>rmsd-interface_r_xray.disp
    refe $refe
    mobi $strucr
    `cat $izoneA`
    atom $atoms
    fit
    quit
_Eod_

  $PROFIT <<_Eod_ |grep RMS |tail -n1 >>rmsd-interface_l_xray.disp
    refe $refe
    mobi $strucl
    `cat $izoneB`
    atom $atoms
    fit
    quit
_Eod_

end

echo "#struc i-RMSD" >i-RMSD.dat
awk '{if ($1 == "RMS:") {printf "%8.3f ",$2} else {printf "\n %s ",$1}}' rmsd-interface_xray.disp |grep pdb |awk '{print $1,$2}' >> i-RMSD.dat
head -n1 i-RMSD.dat >i-RMSD-sorted.dat
grep pdb i-RMSD.dat |sort -nk2  >> i-RMSD-sorted.dat
\rm rmsd-interface_xray.disp

echo "#struc i-RMSD" >i-RMSD_r.dat
awk '{if ($1 == "RMS:") {printf "%8.3f ",$2} else {printf "\n %s ",$1}}' rmsd-interface_r_xray.disp |grep pdb |awk '{print $1,$2}' >> i-RMSD_r.dat
head -n1 i-RMSD_r.dat >i-RMSD_r-sorted.dat
grep pdb i-RMSD_r.dat |sort -nk2  >> i-RMSD_r-sorted.dat
\rm rmsd-interface_r_xray.disp

echo "#struc i-RMSD" >i-RMSD_l.dat
awk '{if ($1 == "RMS:") {printf "%8.3f ",$2} else {printf "\n %s ",$1}}' rmsd-interface_l_xray.disp |grep pdb |awk '{print $1,$2}' >> i-RMSD_l.dat
head -n1 i-RMSD_l.dat >i-RMSD_l-sorted.dat
grep pdb i-RMSD_l.dat |sort -nk2  >> i-RMSD_l-sorted.dat
\rm rmsd-interface_l_xray.disp

exit:
