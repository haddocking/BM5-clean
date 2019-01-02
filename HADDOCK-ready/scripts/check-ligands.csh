#!/bin/csh
#
foreach i ($argv)
  cd $i/
  foreach j (*_u.pdb)
    if (`grep HETATM $j:r.info | awk '{print $NF}'` > 0 ) then
      foreach k (`grep HETATM $j |sed -e 's/ATM/\ /g' |awk '{print $4}' |sort |uniq`)
        set lname=$k
        echo "Detected ligand "$lname" in "$j
      end 
    endif
  end
  cd ../
end
