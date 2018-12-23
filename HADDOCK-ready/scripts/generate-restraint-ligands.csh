#!/bin/csh
#
foreach i ($argv)
  cd $i/
  cat /dev/null >restraint-ligand.tbl
  foreach j (*_u.pdb)
    if (`grep HETATM $j:r.info | awk '{print $NF}'` > 0 ) then
      foreach k (`grep HETATM $j |sed -e 's/ATM/\ /g' |awk '{print $4}' |sort |uniq`)
        set lname=$k
        echo "Detected ligand "$lname" in "$j
        ~abonvin/haddock_git/haddock-tools/restrain_ligand.py -l $lname $j >>restraint-ligand.tbl
        if (! -e ligand.top ) then
          echo "Warning missing ligand top and param files for "$j
        endif
      end
    endif
  end
  if ( `wc -l restraint-ligand.tbl | awk '{print $1}'` == 0 ) then
    \rm restraint-ligand.tbl
  endif
  cd ../
end
