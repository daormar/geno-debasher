# *- bash -*

load_panpipe_module "bam_analysis"

copy_debug_pipeline()
{
    add_panpipe_process "copy_norm_bam"         "cpus=1  mem=2048    time=24:00:00"
    add_panpipe_process "copy_tum_bam"          "cpus=1  mem=2048    time=24:00:00"
    add_panpipe_process "index_norm_bam"        "cpus=1  mem=1024    time=4:00:00"
    add_panpipe_process "index_tum_bam"         "cpus=1  mem=1024    time=4:00:00"
    add_panpipe_process "create_genref_for_bam" "cpus=1  mem=8G      time=4:00:00"
}
