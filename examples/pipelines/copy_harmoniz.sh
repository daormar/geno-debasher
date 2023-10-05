# *- bash -*

load_panpipe_module "bam_analysis"

copy_harmoniz_pipeline()
{
    add_panpipe_process "copy_norm_bam"         "cpus=1  mem=2048 time=48:00:00" "processdeps=none"
    add_panpipe_process "copy_tum_bam"          "cpus=1  mem=2048 time=48:00:00" "processdeps=none"
    add_panpipe_process "norm_bam_to_ubam"      "cpus=1  mem=4096 time=48:00:00" "processdeps=afterok:copy_norm_bam"
    add_panpipe_process "align_norm_ubam"       "cpus=6  mem=8192 time=48:00:00" "processdeps=afterok:norm_bam_to_ubam"
    add_panpipe_process "tum_bam_to_ubam"       "cpus=1  mem=4096 time=48:00:00" "processdeps=afterok:copy_tum_bam"
    add_panpipe_process "align_tum_ubam"        "cpus=6  mem=8192 time=48:00:00" "processdeps=afterok:tum_bam_to_ubam"
    add_panpipe_process "index_norm_bam"        "cpus=1  mem=1024 time=48:00:00" "processdeps=afterok:align_norm_ubam"
    add_panpipe_process "index_tum_bam"         "cpus=1  mem=1024 time=48:00:00" "processdeps=afterok:align_tum_ubam"
    add_panpipe_process "strelka_germline"      "cpus=4  mem=6G   time=48:00:00" "processdeps=afterok:index_norm_bam"
    add_panpipe_process "platypus_germline"     "cpus=1  mem=4096 time=48:00:00" "processdeps=afterok:index_norm_bam"
    add_panpipe_process "gatk_haplotypecaller"  "cpus=4  mem=4096 time=48:00:00" "processdeps=afterok:index_norm_bam"
    add_panpipe_process "manta_somatic"         "cpus=2  mem=4096 time=48:00:00" "processdeps=afterok:index_norm_bam,afterok:index_tum_bam"
    add_panpipe_process "strelka_somatic"       "cpus=4  mem=6G   time=48:00:00" "processdeps=afterok:index_norm_bam,afterok:index_tum_bam"
    add_panpipe_process "mutect2_somatic"       "cpus=4  mem=8G   time=48:00:00" "processdeps=afterok:index_norm_bam,afterok:index_tum_bam"
    add_panpipe_process "lofreq_somatic"        "cpus=4  mem=4096 time=48:00:00" "processdeps=afterok:index_norm_bam,afterok:index_tum_bam"
    add_panpipe_process "msisensor"             "cpus=2  mem=4096 time=48:00:00" "processdeps=afterok:index_norm_bam,afterok:index_tum_bam"
    add_panpipe_process "concat_germline_snvs"  "cpus=1  mem=2048 time=48:00:00" "processdeps=afterok:strelka_germline,afterok:platypus_germline"
}
