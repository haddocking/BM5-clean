#!/bin/csh
#
source /home/abonvin/haddock2.4-node/haddock_configure.csh
set runname=cm
set maxrun=1
set counter=0
set wdir=`pwd`
foreach i ($argv)
  if ($counter < $maxrun) then
    echo "HADDOCK_DIR=/home/abonvin/haddock_git/haddock2.4-node" >${i}/run.param
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
      sed -e 's/ssub\ short/csh/g' -e 's/cpunumber_1=200/cpunumber_1=50/g' run.cns > run.cns.node
      \mv run.cns run.cns.store
      \mv run.cns.node run.cns
      cd ..
      echo "#\!/bin/csh" >${i}-run1-$runname.csh
      echo "#" >>${i}-run1-$runname.csh
      echo "#PBS -S /bin/csh" >>${i}-run1-$runname.csh
      echo "#PBS -q medium" >>${i}-run1-$runname.csh
      echo "#PBS -l nodes=1:ppn=48" >>${i}-run1-$runname.csh
      echo "source /home/abonvin/haddock2.4-node/haddock_configure.csh" >>${i}-run1-$runname.csh
      echo "cd "$wdir/${i}/run1-$runname >>${i}-run1-$runname.csh
      echo "haddock2.4 >&haddock.out" >>${i}-run1-$runname.csh
      qsub ${i}-run1-$runname.csh
      @ counter+=1
      cd ..
    else
      if (! -e $i/run1-$runname/structures/it1/water/analysis/DONE ) then
        if (`grep -i finishing $i/run1-$runname/haddock.out |wc -l |awk '{print $1}'` > 0 ) then
          echo "Detecting crashed run for "$i" "$runname" - restarting"
          cd $i/run1-$runname/
	  ./tools/haddock-clean
	  cd ..
          qsub ${i}-run1-$runname.csh
          @ counter+=1
	  cd ..
	else
          if (`grep -i thread $i/run1-$runname/haddock.out |grep -iv exception |grep -v File |wc -l |awk '{print $1}'` > 0 ) then
            echo "Detecting crashed run for "$i" "$runname" - restarting"
            cd $i/run1-$runname/
	    ./tools/haddock-clean
	    cd ..
            qsub ${i}-run1-$runname.csh
            @ counter+=1
	    cd ..
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
