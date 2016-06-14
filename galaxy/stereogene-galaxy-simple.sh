#!/bin/bash
#$1 is chrom file
#$2 is query
#$3 is query type
#$4 is reference
#$5 is reference type 
#$6 is statistics as Galaxy asks for
#$7 is statistics as Galaxy asks for
rm -f statictics
#not to increse old statistics, we start new one!
qname=${2}.${3}
rname=${4}.${5}
#stereogene looks at the extenson of the file
ln $2 $qname
ln $4 $rname
#debug
#echo StereoGene -s -chrom $1 $qname $rname >> ~/tmp/log-galaxy/simplecmd.log
StereoGene -chrom $1 $rname $qname > out.txt
unlink $qname
unlink $rname
mv out.txt $6
mv ~/example/fw4.pdf $7
