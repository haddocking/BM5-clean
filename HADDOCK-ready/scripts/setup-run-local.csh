#!/bin/csh
#
source /home/abonvin/haddock2.4/haddock_configure.csh
set runname=ti5
set maxrun=1
set counter=0
foreach i ($argv)
  if ($counter < $maxrun) then
    echo "AMBIG_TBL=./ambig5.tbl" >${i}/run.param
    echo "HADDOCK_DIR=/home/abonvin/haddock_git/haddock2.4" >>${i}/run.param
    echo "N_COMP=2" >>${i}/run.param
    echo "PDB_FILE1=./"${i}_r_u.pdb >>${i}/run.param
    echo "PDB_FILE2=./"${i}_l_u.pdb >>${i}/run.param
    echo "PROJECT_DIR=./" >>${i}/run.param
    echo "PROT_SEGID_1=A" >>${i}/run.param
    echo "PROT_SEGID_2=B" >>${i}/run.param
    echo "RUN_NUMBER=1-"$runname >>${i}/run.param
    if ( -e ${i}/hbonds.tbl ) then
      echo "HBOND_FILE=hbonds.tbl" >>${i}/run.param
    endif
    if (! -e ${i}/run1-$runname ) then
      echo "Launching run for "$i" "$runname
      cd $i
      haddock2.4
      \cp ligand* run1-$runname/toppar >&/dev/null
      cd run1-$runname
      patch -p0 -i ../../data/run.cns.patch-$runname
      haddock2.4 >&haddock.out &
      @ counter+=1
      cd ../../
    else
      if (! -e $i/run1-$runname/structures/it1/water/analysis/DONE ) then
        if (`grep -i finishing $i/run1-$runname/haddock.out |wc -l |awk '{print $1}'` > 0 ) then
          echo "Detecting crashed run for "$i" "$runname" - restarting"
          cd $i/run1-$runname/
	  ./tools/haddock-clean
          haddock2.4 >&haddock.out &
          @ counter+=1
	  cd ../..
	else
          if (`grep -i thread $i/run1-$runname/haddock.out |grep -iv exception |grep -v File |wc -l |awk '{print $1}'` > 0 ) then
            echo "Detecting crashed run for "$i" "$runname" - restarting"
            cd $i/run1-$runname/
	    ./tools/haddock-clean
            haddock2.4 >&haddock.out &
            @ counter+=1
	    cd ../..
          else  
            echo "Run for "$i" "$runname" already launched"
          endif
	endif
      else
        echo "Run for "$i" "$runname" finished"
      endif
    endif
  else
    echo "Maximum number of runs ("$maxrun") reached - stopping"
    exit
  endif
end
