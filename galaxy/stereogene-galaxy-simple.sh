#!/bin/bash
#$1 is chrom file
#$2 is query
#$3 is query type
#$4 is reference
#$5 is reference type 
#$6 is window size 
#$7 is pdf name as Galaxy asks for
qname=${2}.${3}
rname=${4}.${5}
#stereogene looks at the extenson of the file
ln $2 $qname
ln $4 $rname
#debug
#echo StereoGene -s -chrom $1 $qname $rname >> ~/tmp/log-galaxy/simplecmd.log
StereoGene -chrom $1 -wSize $6 -R $rname $qname > out.txt
outname=`awk -e '/^out=\"(.+)\"/{match($_,"^out=\"(.+)\"",a);print a[1]}' out.txt`
Rscript ${outname}.r
unlink $qname
unlink $rname
mv ${outname}.pdf $7
