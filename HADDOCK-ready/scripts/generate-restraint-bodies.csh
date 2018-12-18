#!/bin/csh
#
foreach i ($argv)
  cd $i/
  cat /dev/null >restraint-bodies.tbl
  foreach j (*_u.pdb)
    ~abonvin/haddock_git/haddock-tools/restrain_bodies.py $j >>restraint-bodies.tbl
  end
  if ( `wc -l restraint-bodies.tbl | awk '{print $1}'` == 0 ) then
    \rm restraint-bodies.tbl
  endif
  cd ../
end
