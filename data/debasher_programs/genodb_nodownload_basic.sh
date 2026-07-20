# *- bash -*

load_debasher_module "genodb_bam_analysis"

genodb_nodownload_basic_program()
{
    add_debasher_process "manta_somatic"   "cpus=2  mem=4096 time=48:00:00" "processdeps=none"
    add_debasher_process "strelka_somatic" "cpus=2  mem=4096 time=48:00:00" "processdeps=none"
    add_debasher_process "msisensor_pro"   "cpus=2  mem=4096 time=48:00:00" "processdeps=none"
}
