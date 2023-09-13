#import bam_analysis
#
download_aws_norm_bam cpus=1  mem=2048 time=48:00:00 processdeps=none
download_aws_tum_bam  cpus=1  mem=2048 time=48:00:00 processdeps=none
index_norm_bam        cpus=1  mem=1024 time=48:00:00 processdeps=afterok:download_aws_norm_bam
index_tum_bam         cpus=1  mem=1024 time=48:00:00 processdeps=afterok:download_aws_tum_bam
manta_somatic         cpus=8  mem=8G   time=48:00:00 processdeps=afterok:index_norm_bam,afterok:index_tum_bam
strelka_somatic       cpus=8  mem=4096 time=48:00:00 processdeps=afterok:index_norm_bam,afterok:index_tum_bam,afterok:manta_somatic
msisensor             cpus=2  mem=10G  time=48:00:00 processdeps=afterok:index_norm_bam,afterok:index_tum_bam
clear_datadir         cpus=1  mem=1024 time=48:00:00 processdeps=afterok:manta_somatic,afterok:strelka_somatic
