#!/bin/csh
#
foreach i ( 1IRAtrunc_r_*.info )
  set numchain=`grep chains $i |awk '{print $3}'`
  if ( $numchain == 1 ) then
    pdb_chain -A $i:r.pdb |pdb_chainxseg > $i:r.pdbA
  endif
end
foreach i ( 1IRAtrunc_l_*.info )
  set numchain=`grep chains $i |awk '{print $3}'`
  if ( $numchain == 1 ) then
    pdb_chain -B $i:r.pdb |pdb_chainxseg > $i:r.pdbB
  endif
end

