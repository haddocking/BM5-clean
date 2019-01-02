#!/bin/csh
#
source ~abonvin/haddock2.4/haddock_configure.csh
set WDIR=/home/abonvin/docking/BM5-clean/HADDOCK-ready/1XQS/ana_scripts
cd $1/structures/it1/water/analysis
if ( -e cluster.out.gz ) then
  gunzip cluster.out.gz
endif
\rm *out.gz
cd ..
$WDIR/i-rmsd_to_xray.csh 
$WDIR/l-rmsd_to_xray.csh 
$HADDOCKTOOLS/ana_structures.csh
$WDIR/fraction-native.csh 
$HADDOCKTOOLS/ana_clusters.csh -best 4 analysis/cluster.out
$WDIR/cluster-fnat.csh 4
cd ../analysis
if ( -e cluster.out.gz ) then
  gunzip cluster.out.gz
endif
\rm *out.gz
cd ..
$WDIR/i-rmsd_to_xray.csh 
$WDIR/l-rmsd_to_xray.csh 
$HADDOCKTOOLS/ana_structures.csh
$WDIR/fraction-native.csh 
$HADDOCKTOOLS/ana_clusters.csh -best 4 analysis/cluster.out
$WDIR/cluster-fnat.csh 4
cd ../it0
$WDIR/i-rmsd_to_xray.csh 
$WDIR/l-rmsd_to_xray.csh 
$HADDOCKTOOLS/ana_structures.csh
$WDIR/fraction-native.csh 
cd ../../../





