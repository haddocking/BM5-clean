#!/bin/csh

foreach i ($argv)
if ( -d $i ) then
set curdir=`pwd`
cd $i
echo "Now cleaning "$i
echo "cleaning up directories..."
./tools/haddock-clean
gunzip structures/it1/analysis/cluster.out.gz structures/it1/water/analysis/cluster.out.gz >&/dev/null
\rm  structures/it1/analysis/*out.gz structures/it1/water/analysis/*out.gz >&/dev/null
\rm  structures/it1/analysis/*crd structures/it1/water/analysis/*crd >&/dev/null
\rm  structures/it1/analysis/*pdb structures/it1/water/analysis/*pdb >&/dev/null
\rm  begin/*out begin/*_1 >&/dev/null
\rm  structures/it0/*_[1-9] structures/it1/*_[1-9] structures/it1/water/*_[1-9] >&/dev/null
\rm  structures/it0/*[0-1].seed* >&/dev/null
\rm  structures/it0/*[2-3].seed* >&/dev/null
\rm  structures/it0/*[4-5].seed* >&/dev/null
\rm  structures/it0/*[6-7].seed* >&/dev/null
\rm  structures/it0/*[8-9].seed* >&/dev/null
\rm  structures/it1/*[0-1].seed* >&/dev/null
\rm  structures/it1/*[2-3].seed* >&/dev/null
\rm  structures/it1/*[4-5].seed* >&/dev/null
\rm  structures/it1/*[6-7].seed* >&/dev/null
\rm  structures/it1/*[8-9].seed* >&/dev/null
\rm  structures/it0/*[0-1].disp* >&/dev/null
\rm  structures/it0/*[2-3].disp* >&/dev/null
\rm  structures/it0/*[4-5].disp* >&/dev/null
\rm  structures/it0/*[6-7].disp* >&/dev/null
\rm  structures/it0/*[8-9].disp* >&/dev/null
\rm  structures/it1/*[0-1].disp* >&/dev/null
\rm  structures/it1/*[2-3].disp* >&/dev/null
\rm  structures/it1/*[4-5].disp* >&/dev/null
\rm  structures/it1/*[6-7].disp* >&/dev/null
\rm  structures/it1/*[8-9].disp* >&/dev/null
\rm  structures/it1/water/*.pdb_* >&/dev/null
\rm  structures/it1/*.pdb_* >&/dev/null
\rm  structures/it0/*.pdb_* >&/dev/null
\rm  *.job.* >&/dev/null
gzip -f *.out >&/dev/null
echo "now compressing all tructures..."
gzip -f  begin/*.pdb structures/it1/*.pdb structures/it1/water/*.pdb structures/it1/analysis/*.pdb structures/it1/water/analysis*/*.pdb >&/dev/null
gzip -f  structures/it0/*[0-1].pdb >&/dev/null
gzip -f  structures/it0/*[2-3].pdb >&/dev/null
gzip -f  structures/it0/*[4-5].pdb >&/dev/null
gzip -f  structures/it0/*[6-7].pdb >&/dev/null
gzip -f  structures/it0/*[8-9].pdb >&/dev/null
cd $curdir
endif
end
  
