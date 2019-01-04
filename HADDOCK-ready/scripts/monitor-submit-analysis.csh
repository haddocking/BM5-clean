#!/bin/csh
#
cd /home/abonvin/docking/BM5
set ncount=0
while ($ncount < 10)
  @ ncount+=1
  set numjobs=`qstat -u abonvin |grep anal |wc -l |awk '{print $1}'`
  if ($numjobs < 10) then
    date
    ./setup-analysis-ti.csh [1-9BC]* |grep "Launching ana" 
  endif
  set numjobs=`qstat -u abonvin |grep anal |wc -l |awk '{print $1}'`
  if ($numjobs < 10) then
    date
    ./setup-analysis-re.csh [1-9BC]* |grep "Launching ana" 
  endif
#  set numjobs=`qstat -u abonvin |grep anal |wc -l |awk '{print $1}'`
#  if ($numjobs < 10) then
#    date
#    ./setup-analysis-cm.csh [1-9BC]* |grep "Launching ana" 
#  endif
#  set numjobs=`qstat -u abonvin |grep anal |wc -l |awk '{print $1}'`
#  if ($numjobs < 10) then
#    date
#    ./setup-analysis-ra.csh [1-9BC]* |grep "Launching ana" 
#  endif
end

