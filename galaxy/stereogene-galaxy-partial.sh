#!/bin/bash
#$1 is chrom file
#$2 is query file
#$3 is query type
#$4 is reference file 
#$5 is reference type 
#$6 is window size
#$7 is query long Galaxy name
#$8 is reference long Galaxy name
#$9 is pdf name as Galaxy asks for
#$10 is NONE/BASE wigg?
#$11 is wigg name for correlation track
#$12 is confounder file
#$13 is confounder type 
#$14 is confounder long Galaxy name
qname=${2}.${3}
rname=${4}.${5}
conname=${12}.${13}
#stereogene looks at the extenson of the file
ln -f ${2} $qname
ln -f ${4} $rname
ln -f ${12} $conname
#debug
#echo StereoGene -s -chrom $1 $qname $rname >> ~/tmp/log-galaxy/simplecmd.log
StereoGene -chrom $1 -wSize $6 -outWig ${10} -R -pcorProfile $conname $rname $qname > out.txt
outname=`awk -e '/^out=\"(.+)\"/{match($_,"^out=\"(.+)\"",a);print a[1]}' out.txt`
Rscript ${outname}_report.r "\"${7}\"" "\"${8}\"" "\"${14}\"" 2>&1 > /dev/null
unlink $qname
unlink $rname
unlink $conname
mv ${outname}.html $9
if [ ${10} != 'NONE' ] 
then
  mv ${outname}.wig ${11}
fi
