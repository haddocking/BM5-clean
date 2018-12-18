#!/bin/csh
#
set echo
cd /home/abonvin/docking/BM5
set numjobs=`ps -ef |grep python |grep abonvin |grep grid |wc -l |awk '{print $1}'`
if ($numjobs < 20) then
  date
  ./setup-run-grid.csh [1-9BC]* |grep "Launching run" 
endif
