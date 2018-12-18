#!/bin/csh
#
foreach i ($argv)
  set rootname=`echo $i |sed -e 's/_r_u\.pdb//g'`
  set dirname=../docking-ready/$rootname
  mkdir $dirname
  mkdir $dirname/ana_scripts
  \cp ${rootname}_r_b-matched.pdb $dirname
  \cp ${rootname}_l_b-matched.pdb $dirname
  \cp ${rootname}_r_u.pdb $dirname
  \cp ${rootname}_l_u.pdb $dirname
  cat ${rootname}_r_b-matched.pdb ${rootname}_l_b-matched.pdb |pdb_tidy |pdb_chainxseg> $dirname/ana_scripts/target.pdb
  cat ${rootname}_r_u.pdb ${rootname}_l_u.pdb |pdb_tidy |pdb_chainxseg> $dirname/ana_scripts/target-unbound.pdb
  cd $dirname/ana_scripts
  ~/haddock2.4/tools/contact-chainID target.pdb 3.9 >target.contacts3.9
  ~/haddock2.4/tools/contact-chainID target.pdb 5 >target.contacts5
  ~/haddock2.4/tools/contact-chainID target.pdb 10 >target.contacts10
  ~/haddock2.4/tools/contact-chainID target.pdb 200 >target.contacts200
  ~/haddock_git/haddock-CSB-tools/CAPRI_analysis/make_izone.csh target.contacts10 >tmp
  ~/haddock_git/haddock-CSB-tools/CAPRI_analysis/join_zones.py tmp >target.izone
  ~/haddock_git/haddock-CSB-tools/CAPRI_analysis/make_izone.csh target.contacts200 >tmp
  ~/haddock_git/haddock-CSB-tools/CAPRI_analysis/join_zones.py tmp >tt
  grep "zone\ A" tt >target.lzone
  echo "fit" >>target.lzone
  grep "zone\ B" tt | sed -e 's/zone/rzone/g' >>target.lzone
  \rm tmp tt target.contacts200
  cd ../../../matched
end
