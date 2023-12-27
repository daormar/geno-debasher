# *- bash -*

load_debasher_module "genop_bam_analysis"

genop_copy_debug_program()
{
    add_debasher_process "copy_norm_bam"         "cpus=1  mem=2048    time=24:00:00"
    add_debasher_process "copy_tum_bam"          "cpus=1  mem=2048    time=24:00:00"
    add_debasher_process "index_norm_bam"        "cpus=1  mem=1024    time=4:00:00"
    add_debasher_process "index_tum_bam"         "cpus=1  mem=1024    time=4:00:00"
    add_debasher_process "create_genref_for_bam" "cpus=1  mem=8G      time=4:00:00"
}
