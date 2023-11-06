# *- bash -*

load_panpipe_module "genop_bam_analysis"

ega_asp_basic_pipeline()
{
    add_panpipe_process "download_ega_asp_norm_bam" "cpus=1  mem=2048 time=48:00:00" "processdeps=none"
    add_panpipe_process "download_ega_asp_tum_bam"  "cpus=1  mem=2048 time=48:00:00" "processdeps=none"
    add_panpipe_process "index_norm_bam"            "cpus=1  mem=1024 time=48:00:00" "processdeps=afterok:download_ega_asp_norm_bam"
    add_panpipe_process "index_tum_bam"             "cpus=1  mem=1024 time=48:00:00" "processdeps=afterok:download_ega_asp_tum_bam"
    add_panpipe_process "manta_somatic"             "cpus=2  mem=8G   time=48:00:00" "processdeps=afterok:index_norm_bam,afterok:index_tum_bam"
    add_panpipe_process "strelka_somatic"           "cpus=2  mem=4096 time=48:00:00" "processdeps=afterok:index_norm_bam,afterok:index_tum_bam,afterok:manta_somatic"
    add_panpipe_process "msisensor"                 "cpus=2  mem=10G  time=48:00:00" "processdeps=afterok:index_norm_bam,afterok:index_tum_bam"
    add_panpipe_process "clear_datadir"             "cpus=1  mem=1024 time=48:00:00" "processdeps=afterok:manta_somatic,afterok:strelka_somatic,afterok:msisensor"
}
