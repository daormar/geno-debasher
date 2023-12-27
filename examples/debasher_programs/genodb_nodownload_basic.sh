# *- bash -*

load_debasher_module "genop_bam_analysis"

genop_nodownload_basic_program()
{
    add_debasher_process "manta_somatic"   "cpus=2  mem=4096 time=48:00:00" "processdeps=none"
    add_debasher_process "strelka_somatic" "cpus=2  mem=4096 time=48:00:00" "processdeps=none"
    add_debasher_process "msisensor"       "cpus=2  mem=4096 time=48:00:00" "processdeps=none"
}
