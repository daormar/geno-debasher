#import bam_analysis
#
copy_norm_bam               cpus=1  mem=2048    time=24:00:00 processdeps=none
copy_tum_bam                cpus=1  mem=2048    time=24:00:00 processdeps=none
index_norm_bam              cpus=1  mem=1024    time=4:00:00  processdeps=afterok:copy_norm_bam
index_tum_bam               cpus=1  mem=1024    time=4:00:00  processdeps=afterok:copy_tum_bam
create_genref_for_bam       cpus=1  mem=8G      time=4:00:00  processdeps=afterok:index_norm_bam
