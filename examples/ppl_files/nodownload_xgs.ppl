#import bam_analysis
#
create_genref_for_bam       cpus=1  mem=8G      time=4:00:00  processdeps=none
bedtools_genomecov_norm_bam cpus=1  mem=1024    time=4:00:00  processdeps=none
bedtools_genomecov_tum_bam  cpus=1  mem=1024    time=4:00:00  processdeps=none
strelka_germline            cpus=8  mem=6G      time=6:00:00,12:00:00  processdeps=afterok:create_genref_for_bam
manta_somatic               cpus=8  mem=8G      time=8:00:00,16:00:00  processdeps=afterok:create_genref_for_bam
strelka_somatic             cpus=8  mem=6G      time=8:00:00,16:00:00  processdeps=afterok:create_genref_for_bam,afterok:manta_somatic
msisensor                   cpus=2  mem=8G      time=8:00:00,24:00:00  processdeps=afterok:create_genref_for_bam
snp_pileup_plus_facets      cpus=1  mem=8G      time=8:00:00,24:00:00  processdeps=afterok:create_genref_for_bam
platypus_germline           cpus=1  mem=4096    time=8:00:00,16:00:00  processdeps=afterok:create_genref_for_bam
