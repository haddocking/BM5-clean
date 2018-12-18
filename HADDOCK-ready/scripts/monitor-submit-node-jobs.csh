#!/bin/csh
#
cd /home/abonvin/docking/BM5
set numjobs=`qstat -u abonvin |grep run1-cm |wc -l |awk '{print $1}'`
if ($numjobs < 10) then
  date
  ./setup-run-node.csh [1-7BC]* |grep "Launching run" 
endif
