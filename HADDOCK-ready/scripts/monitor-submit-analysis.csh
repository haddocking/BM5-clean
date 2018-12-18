#!/bin/csh
#
cd /home/abonvin/docking/BM5
set numjobs=`qstat -u abonvin |grep anal |wc -l |awk '{print $1}'`
if ($numjobs < 10) then
  date
  ./setup-analysis.csh [1-9BC]* |grep "Launching ana" 
endif
