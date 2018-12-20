#!/bin/csh
#
foreach i ($argv)
  cd $i/ana_scripts
  awk '{print $1}' target.contacts3.9 | sort -n | uniq | awk '{printf "%s+",$1}' > interface_r.lis
  awk '{print $4}' target.contacts3.9 | sort -n | uniq | awk '{printf "%s+",$1}' > interface_l.lis
  ../../scripts/generate-act-act.sh A `cat interface_r.lis` B `cat interface_l.lis` >../ambig.tbl
  cd ../../
end
