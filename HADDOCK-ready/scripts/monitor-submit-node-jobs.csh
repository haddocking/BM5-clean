#!/bin/csh
#
cd /home/abonvin/docking/BM5
set numjobs=`qstat -u abonvin |grep run1-cm |wc -l |awk '{print $1}'`
if ($numjobs < 12) then
  date
  ./setup-run-node-cm.csh [1-9BC]* |grep "Launching run" 
#else
#  date
#  echo "Maximum number of runs reached"
endif
