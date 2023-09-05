#!/bin/bash


cbusco3


# create config files:
for ASSEMBLY in *fasta ; \
do \
pd_config.py \
-n ${ASSEMBLY%%.fasta}.config.json \
-s /home/ubuntu/avas/haplomerging/Aphrocallistes_short-read.list \
$ASSEMBLY \
/home/ubuntu/avas/haplomerging/Aphrocallistes_long-read.list ; \
done

for file in *json ; do sed -i 's/\"core\"\: 12/\"core\"\: 30/' $file ; done
for file in *json ; do sed -i 's/\"mem\"\: 10000/\"mem\"\: 100000/' $file ; done
for file in *json ; do sed -i 's/\"mem\"\: 20000/\"mem\"\: 100000/' $file ; done
for file in *json ; do sed -i 's/\"mem\"\: 30000/\"mem\"\: 100000/' $file ; done



# haplo-purging
for ASSEMBLY in *fasta ; \
do \
run_purge_dups.py -p bash \
${ASSEMBLY%.fasta}.config.json \
/home/ubuntu/tools/purge_dups/bin \
Aphrocallistes ; \
done


cdeactivate