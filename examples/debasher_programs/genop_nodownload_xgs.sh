# *- bash -*

load_panpipe_module "genop_bam_analysis"

genop_nodownload_xgs_program()
{
    add_panpipe_process "create_genref_for_bam"       "cpus=1  mem=8G      time=4:00:00"  "processdeps=none"
    add_panpipe_process "bedtools_genomecov_norm_bam" "cpus=1  mem=1024    time=4:00:00"  "processdeps=none"
    add_panpipe_process "bedtools_genomecov_tum_bam"  "cpus=1  mem=1024    time=4:00:00"  "processdeps=none"
    add_panpipe_process "strelka_germline"            "cpus=8  mem=6G      time=6:00:00,12:00:00"  "processdeps=afterok:create_genref_for_bam"
    add_panpipe_process "manta_somatic"               "cpus=8  mem=8G      time=8:00:00,16:00:00"  "processdeps=afterok:create_genref_for_bam"
    add_panpipe_process "strelka_somatic"             "cpus=8  mem=6G      time=8:00:00,16:00:00"  "processdeps=afterok:create_genref_for_bam,afterok:manta_somatic"
    add_panpipe_process "msisensor"                   "cpus=2  mem=8G      time=8:00:00,24:00:00"  "processdeps=afterok:create_genref_for_bam"
    add_panpipe_process "snp_pileup_plus_facets"      "cpus=1  mem=8G      time=8:00:00,24:00:00"  "processdeps=afterok:create_genref_for_bam"
    add_panpipe_process "platypus_germline"           "cpus=1  mem=4096    time=8:00:00,16:00:00"  "processdeps=afterok:create_genref_for_bam"
}
