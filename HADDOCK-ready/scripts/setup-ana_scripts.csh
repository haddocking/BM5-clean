#!/bin/csh
#
foreach i ($argv)
  cd $i/ana_scripts
  grep A target.izone >target.izoneA
  grep B target.izone >target.izoneB
  foreach j (../../ana_scripts/*csh)
     echo "#\!/bin/csh" > `basename $j`
     echo "#" >> `basename $j`
     echo "source ~abonvin/haddock2.4/haddock_configure.csh" >> `basename $j`
     echo "set WDIR="`pwd` >> `basename $j`
     cat $j >> `basename $j`
  end
  chmod +x *csh
  cd ../../
end
