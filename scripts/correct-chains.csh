#!/bin/csh
#
set newchain=A
foreach i ( $argv )
  set numchain=`grep chains $i |awk '{print $3}'`
  if ( $numchain == 2 ) then
    set chain1=`grep ATOM $i:r.pdb |head -n1 |awk '{print substr($0,22,1)}'`
    set chain2=`grep ATOM $i:r.pdb |tail -n1 |awk '{print substr($0,22,1)}'`
    set start1=`pdb_selchain -$chain1 $i:r.pdb |grep ATOM |head -n1 |awk '{print substr($0,23,4)}'`
    set start2=`pdb_selchain -$chain2 $i:r.pdb |grep ATOM |head -n1 |awk '{print substr($0,23,4)}'`
    set end1=`pdb_selchain -$chain1 $i:r.pdb |grep ATOM |tail -n1 |awk '{print substr($0,23,4)}'`
    set end2=`pdb_selchain -$chain2 $i:r.pdb |grep ATOM |tail -n1 |awk '{print substr($0,23,4)}'`
    echo $i $chain1 $chain2 $start1 $end1 $start2 $end2
    if ($end1 < 4000) then
      set shift=4000
    endif
    if ($end1 < 3000) then
      set shift=3000
    endif
    if ($end1 < 2000) then
      set shift=2000
    endif
    if ($end1 < 1000) then
      set shift=1000
    endif
    if ($end1 < 500) then
      set shift=500
    endif
    switch ( "$i" )
      case "*_r_*":
        set newchain=A
        breaksw
      case "*_l_*":
        set newchain=B
        breaksw
    endsw
    echo "setting chain to "$newchain
    pdb_selchain -$chain1 $i:r.pdb |pdb_chain -$newchain > tmp
    pdb_selchain -$chain2 $i:r.pdb |pdb_shiftres -$shift | pdb_chain -$newchain >> tmp
    pdb_tidy tmp >$i:r.pdb-chain-fixed
  endif

  if ( $numchain == 3 ) then
    set chain1=`grep ATOM $i:r.pdb |awk '{print substr($0,22,1)}' |sort |uniq |head -n1`
    set chain2=`grep ATOM $i:r.pdb |awk '{print substr($0,22,1)}' |sort |uniq |head -n2 |tail -n1`
    set chain3=`grep ATOM $i:r.pdb |awk '{print substr($0,22,1)}' |sort |uniq |tail -n1`
    set start1=`pdb_selchain -$chain1 $i:r.pdb |grep ATOM |head -n1 |awk '{print substr($0,23,4)}'`
    set start2=`pdb_selchain -$chain2 $i:r.pdb |grep ATOM |head -n1 |awk '{print substr($0,23,4)}'`
    set start3=`pdb_selchain -$chain3 $i:r.pdb |grep ATOM |head -n1 |awk '{print substr($0,23,4)}'`
    set end1=`pdb_selchain -$chain1 $i:r.pdb |grep ATOM |tail -n1 |awk '{print substr($0,23,4)}'`
    set end2=`pdb_selchain -$chain2 $i:r.pdb |grep ATOM |tail -n1 |awk '{print substr($0,23,4)}'`
    set end3=`pdb_selchain -$chain3 $i:r.pdb |grep ATOM |tail -n1 |awk '{print substr($0,23,4)}'`
    echo $i $chain1 $chain2 $chain3 $start1 $end1 $start2 $end2 $start3 $end3
    if ($end1 < 4000) then
      set shift1=4000
    endif
    if ($end1 < 3000) then
      set shift1=3000
    endif
    if ($end1 < 2000) then
      set shift1=2000
    endif
    if ($end1 < 1000) then
      set shift1=1000
    endif
    if ($end1 < 500) then
      set shift1=500
    endif
    if ($end2 < 4000) then
      set shift2=4000
    endif
    if ($end2 < 3000) then
      set shift2=3000
    endif
    if ($end2 < 2000) then
      set shift2=2000
    endif
    if ($end2 < 1000) then
      set shift2=1000
    endif
    if ($end2 < 500) then
      set shift2=500
    endif
    @ shift2+=$shift1
    switch ( "$i" )
      case "*_r_*":
        set newchain=A
        breaksw
      case "*_l_*":
        set newchain=B
        breaksw
    endsw
    echo "setting chain to "$newchain
    pdb_selchain -$chain1 $i:r.pdb |pdb_chain -$newchain > tmp
    pdb_selchain -$chain2 $i:r.pdb |pdb_shiftres -$shift1 | pdb_chain -$newchain >> tmp
    pdb_selchain -$chain3 $i:r.pdb |pdb_shiftres -$shift2 | pdb_chain -$newchain >> tmp
    pdb_tidy tmp >$i:r.pdb-chain-fixed
  endif

  if ( $numchain == 4 ) then
    set chain1=`grep ATOM $i:r.pdb |awk '{print substr($0,22,1)}' |sort |uniq |head -n1`
    set chain2=`grep ATOM $i:r.pdb |awk '{print substr($0,22,1)}' |sort |uniq |head -n2 |tail -n1`
    set chain3=`grep ATOM $i:r.pdb |awk '{print substr($0,22,1)}' |sort |uniq |head -n3 |tail -n1`
    set chain4=`grep ATOM $i:r.pdb |awk '{print substr($0,22,1)}' |sort |uniq |tail -n1`
    set start1=`pdb_selchain -$chain1 $i:r.pdb |grep ATOM |head -n1 |awk '{print substr($0,23,4)}'`
    set start2=`pdb_selchain -$chain2 $i:r.pdb |grep ATOM |head -n1 |awk '{print substr($0,23,4)}'`
    set start3=`pdb_selchain -$chain3 $i:r.pdb |grep ATOM |head -n1 |awk '{print substr($0,23,4)}'`
    set start4=`pdb_selchain -$chain4 $i:r.pdb |grep ATOM |head -n1 |awk '{print substr($0,23,4)}'`
    set end1=`pdb_selchain -$chain1 $i:r.pdb |grep ATOM |tail -n1 |awk '{print substr($0,23,4)}'`
    set end2=`pdb_selchain -$chain2 $i:r.pdb |grep ATOM |tail -n1 |awk '{print substr($0,23,4)}'`
    set end3=`pdb_selchain -$chain3 $i:r.pdb |grep ATOM |tail -n1 |awk '{print substr($0,23,4)}'`
    set end4=`pdb_selchain -$chain4 $i:r.pdb |grep ATOM |tail -n1 |awk '{print substr($0,23,4)}'`
    echo $i $chain1 $chain2 $chain3 $chain4 $start1 $end1 $start2 $end2 $start3 $end3 $start4 $end4
    if ($end1 < 4000) then
      set shift1=4000
    endif
    if ($end1 < 3000) then
      set shift1=3000
    endif
    if ($end1 < 2000) then
      set shift1=2000
    endif
    if ($end1 < 1000) then
      set shift1=1000
    endif
    if ($end1 < 500) then
      set shift1=500
    endif
    if ($end2 < 4000) then
      set shift2=4000
    endif
    if ($end2 < 3000) then
      set shift2=3000
    endif
    if ($end2 < 2000) then
      set shift2=2000
    endif
    if ($end2 < 1000) then
      set shift2=1000
    endif
    if ($end2 < 500) then
      set shift2=500
    endif
    if ($end3 < 4000) then
      set shift3=4000
    endif
    if ($end3 < 3000) then
      set shift3=3000
    endif
    if ($end3 < 2000) then
      set shift3=2000
    endif
    if ($end3 < 1000) then
      set shift3=1000
    endif
    if ($end3 < 500) then
      set shift3=500
    endif
    @ shift2+=$shift1
    @ shift3+=$shift2

    switch ( "$i" )
      case "*_r_*":
        set newchain=A
        breaksw
      case "*_l_*":
        set newchain=B
        breaksw
    endsw
    echo "setting chain to "$newchain
    pdb_selchain -$chain1 $i:r.pdb |pdb_chain -$newchain > tmp
    pdb_selchain -$chain2 $i:r.pdb |pdb_shiftres -$shift1 | pdb_chain -$newchain >> tmp
    pdb_selchain -$chain3 $i:r.pdb |pdb_shiftres -$shift2 | pdb_chain -$newchain >> tmp
    pdb_selchain -$chain4 $i:r.pdb |pdb_shiftres -$shift3 | pdb_chain -$newchain >> tmp
    pdb_tidy tmp >$i:r.pdb-chain-fixed
  endif

  if ( $numchain == 6 ) then
    set chain1=`grep ATOM $i:r.pdb |awk '{print substr($0,22,1)}' |sort |uniq |head -n1`
    set chain2=`grep ATOM $i:r.pdb |awk '{print substr($0,22,1)}' |sort |uniq |head -n2 |tail -n1`
    set chain3=`grep ATOM $i:r.pdb |awk '{print substr($0,22,1)}' |sort |uniq |head -n3 |tail -n1`
    set chain4=`grep ATOM $i:r.pdb |awk '{print substr($0,22,1)}' |sort |uniq |head -n4 |tail -n1`
    set chain5=`grep ATOM $i:r.pdb |awk '{print substr($0,22,1)}' |sort |uniq |head -n5 |tail -n1`
    set chain6=`grep ATOM $i:r.pdb |awk '{print substr($0,22,1)}' |sort |uniq |tail -n1`
    set start1=`pdb_selchain -$chain1 $i:r.pdb |grep ATOM |head -n1 |awk '{print substr($0,23,4)}'`
    set start2=`pdb_selchain -$chain2 $i:r.pdb |grep ATOM |head -n1 |awk '{print substr($0,23,4)}'`
    set start3=`pdb_selchain -$chain3 $i:r.pdb |grep ATOM |head -n1 |awk '{print substr($0,23,4)}'`
    set start4=`pdb_selchain -$chain4 $i:r.pdb |grep ATOM |head -n1 |awk '{print substr($0,23,4)}'`
    set start5=`pdb_selchain -$chain5 $i:r.pdb |grep ATOM |head -n1 |awk '{print substr($0,23,4)}'`
    set start6=`pdb_selchain -$chain6 $i:r.pdb |grep ATOM |head -n1 |awk '{print substr($0,23,4)}'`
    set end1=`pdb_selchain -$chain1 $i:r.pdb |grep ATOM |tail -n1 |awk '{print substr($0,23,4)}'`
    set end2=`pdb_selchain -$chain2 $i:r.pdb |grep ATOM |tail -n1 |awk '{print substr($0,23,4)}'`
    set end3=`pdb_selchain -$chain3 $i:r.pdb |grep ATOM |tail -n1 |awk '{print substr($0,23,4)}'`
    set end4=`pdb_selchain -$chain4 $i:r.pdb |grep ATOM |tail -n1 |awk '{print substr($0,23,4)}'`
    set end5=`pdb_selchain -$chain5 $i:r.pdb |grep ATOM |tail -n1 |awk '{print substr($0,23,4)}'`
    set end6=`pdb_selchain -$chain6 $i:r.pdb |grep ATOM |tail -n1 |awk '{print substr($0,23,4)}'`
    echo $i $chain1 $chain2 $chain3 $chain4 $chain5 $chain6 $start1 $end1 $start2 $end2 $start3 $end3 $start4 $end4 $start5 $end5 $start6 $end6
    if ($end1 < 4000) then
      set shift1=4000
    endif
    if ($end1 < 3000) then
      set shift1=3000
    endif
    if ($end1 < 2000) then
      set shift1=2000
    endif
    if ($end1 < 1000) then
      set shift1=1000
    endif
    if ($end1 < 500) then
      set shift1=500
    endif
    if ($end2 < 4000) then
      set shift2=4000
    endif
    if ($end2 < 3000) then
      set shift2=3000
    endif
    if ($end2 < 2000) then
      set shift2=2000
    endif
    if ($end2 < 1000) then
      set shift2=1000
    endif
    if ($end2 < 500) then
      set shift2=500
    endif
    if ($end3 < 4000) then
      set shift3=4000
    endif
    if ($end3 < 3000) then
      set shift3=3000
    endif
    if ($end3 < 2000) then
      set shift3=2000
    endif
    if ($end3 < 1000) then
      set shift3=1000
    endif
    if ($end3 < 500) then
      set shift3=500
    endif
    if ($end4 < 4000) then
      set shift4=4000
    endif
    if ($end4 < 3000) then
      set shift4=3000
    endif
    if ($end4 < 2000) then
      set shift4=2000
    endif
    if ($end4 < 1000) then
      set shift4=1000
    endif
    if ($end4 < 500) then
      set shift4=500
    endif
    if ($end5 < 4000) then
      set shift5=4000
    endif
    if ($end5 < 3000) then
      set shift5=3000
    endif
    if ($end5 < 2000) then
      set shift5=2000
    endif
    if ($end5 < 1000) then
      set shift5=1000
    endif
    if ($end5 < 500) then
      set shift5=500
    endif
    @ shift2+=$shift1
    @ shift3+=$shift2
    @ shift4+=$shift3
    @ shift5+=$shift4

    switch ( "$i" )
      case "*_r_*":
        set newchain=A
        breaksw
      case "*_l_*":
        set newchain=B
        breaksw
    endsw
    echo "setting chain to "$newchain
    pdb_selchain -$chain1 $i:r.pdb |pdb_chain -$newchain > tmp
    pdb_selchain -$chain2 $i:r.pdb |pdb_shiftres -$shift1 | pdb_chain -$newchain >> tmp
    pdb_selchain -$chain3 $i:r.pdb |pdb_shiftres -$shift2 | pdb_chain -$newchain >> tmp
    pdb_selchain -$chain4 $i:r.pdb |pdb_shiftres -$shift3 | pdb_chain -$newchain >> tmp
    pdb_selchain -$chain5 $i:r.pdb |pdb_shiftres -$shift4 | pdb_chain -$newchain >> tmp
    pdb_selchain -$chain6 $i:r.pdb |pdb_shiftres -$shift5 | pdb_chain -$newchain >> tmp
    pdb_tidy tmp >$i:r.pdb-chain-fixed
  endif
end

