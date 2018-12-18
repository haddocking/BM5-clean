#!/bin/csh
#
foreach i ($argv)
  cd $i/
  cat /dev/null >hbonds.tbl
  cat ${i}*.tbl >hbonds.tbl
  cat restraint-bodies.tbl >>hbonds.tbl
  if ( `wc -l hbonds.tbl | awk '{print $1}'` == 0 ) then
    \rm hbonds.tbl
  endif
  cd ../
end
