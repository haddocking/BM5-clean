#!/bin/csh
#
source /home/abonvin/haddock2.4/haddock_configure.csh
set runname=refb
set maxrun=10
set counter=0
set clean=0
foreach i ($argv)
  if ($counter < $maxrun) then
    if ( -e ${i}/run1-$runname/structures/it1/water/analysis/CLUST_DONE && -e ${i}/run1-$runname/structures/it1/analysis/CLUST_DONE ) then
      if ( -e ${i}/run1-$runname/structures/it0/file.nam_fnat ) then
        echo "Analysis for "$i" "$runname" already performed"
	if (! -e ${i}/run1-$runname.stats ) then
	  cd $i
	  ../results-stats.csh run1-$runname > run1-$runname.stats
	  cd ..
	endif
	if ( $clean == 1 ) then
	  \rm ${i}/run1-$runname/structures/it0/*disp >&/dev/null
	  \rm ${i}/run1-$runname/structures/it0/*seed >&/dev/null
	  \rm ${i}/run1-$runname/structures/it1/*seed >&/dev/null
	  \rm ${i}/run1-$runname/structures/it1/analysis/*contacts >&/dev/null
	  \rm ${i}/run1-$runname/structures/it1/analysis/*pdb >&/dev/null
	  \rm ${i}/run1-$runname/structures/it1/analysis/*crd >&/dev/null
	  \rm ${i}/run1-$runname/structures/it1/water/*seed >&/dev/null
	  \rm ${i}/run1-$runname/structures/it1/water/analysis/*contacts >&/dev/null
	  \rm ${i}/run1-$runname/structures/it1/water/analysis/*pdb >&/dev/null
	  \rm ${i}/run1-$runname/structures/it1/water/analysis/*crd >&/dev/null
	endif
      else
        if (! -e ${i}/${i}-analyse-${runname}.csh ) then
          echo "Launching analysis for "$i" "$runname
  	  cd $i
	  echo "#\!/bin/csh" >${i}-analyse-${runname}.csh
	  echo "#" >>${i}-analyse-${runname}.csh
	  echo "source /home/abonvin/haddock2.4/haddock_configure.csh" >>${i}-analyse-${runname}.csh
	  echo "cd "`pwd` >>${i}-analyse-${runname}.csh
          echo "./ana_scripts/run_all.csh run1-"$runname >>${i}-analyse-${runname}.csh
          echo "../results-stats.csh run1-"$runname" > run1-"$runname".stats" >>${i}-analyse-${runname}.csh
	  echo "../clean-run-multi.csh run1-"$runname" >&/dev/null " >>${i}-analyse-${runname}.csh
	  ssub medium ${i}-analyse-${runname}.csh
	  @ counter+=1
	  cd ..
	else
	  echo "Analysis for "$i" "$runname" already launched"
	endif
      endif
    else
      echo "Run for "$i" "$runname" not yet completed"
    endif
  else
    echo "Maximum number of analysis runs ("$maxrun") reached - stopping"
    exit
  endif
end
