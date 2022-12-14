#!/bin/bash




BAMLONG1=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.ONT-RNA_minimap2.redone.sorted.bam
BAMLONG2=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.ONT-RNA_gmap.sorted.redone.k16.bam



# ONT-RNA_minimap2

bamutils removeclipping -f $BAMLONG1 ${BAMLONG1%.bam}.cl.bam
bamutils filter ${BAMLONG1%.bam}.cl.bam ${BAMLONG1%.bam}.cl.f1000.bam -minlen 1000 -mapped -nosecondary -noqcfail
bamutils filter ${BAMLONG1%.bam}.cl.bam ${BAMLONG1%.bam}.cl.f500.bam -minlen 500 -mapped -nosecondary -noqcfail
bamutils filter ${BAMLONG1%.bam}.cl.bam ${BAMLONG1%.bam}.cl.f300.bam -minlen 300 -mapped -nosecondary -noqcfail
bamutils filter ${BAMLONG1%.bam}.cl.bam ${BAMLONG1%.bam}.cl.f100.bam -minlen 100 -mapped -nosecondary -noqcfail



# ONT-RNA_gmap2

bamutils removeclipping -f $BAMLONG2 ${BAMLONG2%.bam}.cl.bam
bamutils filter ${BAMLONG2%.bam}.cl.bam ${BAMLONG2%.bam}.cl.f1000.bam -minlen 1000 -mapped -nosecondary -noqcfail
bamutils filter ${BAMLONG2%.bam}.cl.bam ${BAMLONG2%.bam}.cl.f500.bam -minlen 500 -mapped -nosecondary -noqcfail
bamutils filter ${BAMLONG2%.bam}.cl.bam ${BAMLONG2%.bam}.cl.f300.bam -minlen 300 -mapped -nosecondary -noqcfail
bamutils filter ${BAMLONG2%.bam}.cl.bam ${BAMLONG2%.bam}.cl.f100.bam -minlen 100 -mapped -nosecondary -noqcfail
