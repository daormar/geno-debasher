# *- bash -*

load_debasher_module "genodb_bam_analysis"

genodb_nodownload_basic_program()
{
    add_debasher_process "genodb.manta_somatic"   "cpus=2  mem=4096 time=48:00:00" "processdeps=none"
    add_debasher_process "genodb.strelka_somatic" "cpus=2  mem=4096 time=48:00:00" "processdeps=none"
    add_debasher_process "genodb.msisensor_pro"   "cpus=2  mem=4096 time=48:00:00" "processdeps=none"
}
