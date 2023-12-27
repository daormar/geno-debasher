# *- bash -*

load_debasher_module "genodb_bam_analysis"

genodb_aws_basic_program()
{
    add_debasher_process "download_aws_norm_bam" "cpus=1  mem=2048 time=48:00:00" "processdeps=none"
    add_debasher_process "download_aws_tum_bam"  "cpus=1  mem=2048 time=48:00:00" "processdeps=none"
    add_debasher_process "index_norm_bam"        "cpus=1  mem=1024 time=48:00:00" "processdeps=afterok:download_aws_norm_bam"
    add_debasher_process "index_tum_bam"         "cpus=1  mem=1024 time=48:00:00" "processdeps=afterok:download_aws_tum_bam"
    add_debasher_process "manta_somatic"         "cpus=8  mem=8G   time=48:00:00" "processdeps=afterok:index_norm_bam,afterok:index_tum_bam"
    add_debasher_process "strelka_somatic"       "cpus=8  mem=4096 time=48:00:00" "processdeps=afterok:index_norm_bam,afterok:index_tum_bam,afterok:manta_somatic"
    add_debasher_process "msisensor"             "cpus=2  mem=10G  time=48:00:00" "processdeps=afterok:index_norm_bam,afterok:index_tum_bam"
    add_debasher_process "clear_datadir"         "cpus=1  mem=1024 time=48:00:00" "processdeps=afterok:manta_somatic,afterok:strelka_somatic"
}
