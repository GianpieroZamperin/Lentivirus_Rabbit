# Lentivirus_Rabbit
python script used to filter low quality and duplicated NGS reads

USAGE:
python filter_reads.py read_R1.fastq.gz read_R2.fastq.gz 2>>log | sort -k 2,2 -k 3,3 -T TEMP_SORT | awk -F "\t" 'BEGIN {fil=0; tot=0;} {if ($2"\t"$3 != id) {print "@"$1"\n"$2"\n+\n"$4 |& "gzip > F.nodup.gz"; print "@"$1"\n"$3"\n+\n"$5 |& "gzip > R.nodup.gz"; fil=fil+1} id=$2"\t"$3; tot=tot+1;} END {print "\t"fil" out of "tot" fragments passed the filtering" >> "log";}'
