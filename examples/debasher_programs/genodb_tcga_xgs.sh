# *- bash -*

load_debasher_module "genodb_bam_analysis"

genodb_tcga_xgs_program()
{
    add_debasher_process "download_gdc_norm_bam"       "cpus=1  mem=2048    time=24:00:00"          "processdeps=none"
    add_debasher_process "download_gdc_tum_bam"        "cpus=1  mem=2048    time=24:00:00"          "processdeps=none"
    add_debasher_process "index_norm_bam"              "cpus=1  mem=1024    time=4:00:00"           "processdeps=afterok:download_gdc_norm_bam"
    add_debasher_process "index_tum_bam"               "cpus=1  mem=1024    time=4:00:00"           "processdeps=afterok:download_gdc_tum_bam"
    add_debasher_process "create_genref_for_bam"       "cpus=1  mem=8G      time=4:00:00"           "processdeps=afterok:download_gdc_norm_bam"
    add_debasher_process "bedtools_genomecov_norm_bam" "cpus=1  mem=1024    time=4:00:00"           "processdeps=afterok:index_norm_bam"
    add_debasher_process "bedtools_genomecov_tum_bam"  "cpus=1  mem=1024    time=4:00:00"           "processdeps=afterok:index_tum_bam"
    add_debasher_process "strelka_germline"            "cpus=8  mem=6G      time=6:00:00,12:00:00"  "processdeps=afterok:create_genref_for_bam,afterok:index_norm_bam"
    add_debasher_process "manta_somatic"               "cpus=8  mem=8G      time=8:00:00,16:00:00"  "processdeps=afterok:create_genref_for_bam,afterok:index_norm_bam,afterok:index_tum_bam"
    add_debasher_process "strelka_somatic"             "cpus=8  mem=6G      time=8:00:00,16:00:00"  "processdeps=afterok:create_genref_for_bam,afterok:manta_somatic"
    add_debasher_process "msisensor_pro"               "cpus=2  mem=8G      time=8:00:00,24:00:00"  "processdeps=afterok:create_genref_for_bam,afterok:index_norm_bam,afterok:index_tum_bam"
    add_debasher_process "snp_pileup_plus_facets"      "cpus=1  mem=8G      time=8:00:00,24:00:00"  "processdeps=afterok:create_genref_for_bam,afterok:index_norm_bam,afterok:index_tum_bam"
    add_debasher_process "platypus_germline"           "cpus=1  mem=4096    time=8:00:00,16:00:00"  "processdeps=afterok:create_genref_for_bam,afterok:index_norm_bam"
}
