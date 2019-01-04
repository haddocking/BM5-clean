#!/bin/csh
#
cd /home/abonvin/docking/BM5
set ncount=0
while ($ncount < 30)
  @ ncount+=1
  set numjobs=`ps -ef |grep python |grep abonvin |grep -v grid |grep haddock2 | wc -l |awk '{print $1}'`
  if ($numjobs < 30) then
    date
    ./setup-run-local-re.csh [1-9BC]* |grep "Launching run"
  endif
end
