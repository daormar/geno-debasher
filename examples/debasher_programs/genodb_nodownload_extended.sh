# *- bash -*

load_debasher_module "genodb_bam_analysis"

genodb_nodownload_extended_program()
{
    add_debasher_process "create_genref_for_bam"     "cpus=1  mem=8G      time=4:00:00"  "processdeps=none"
    add_debasher_process "parallel_split_norm_bam"   "cpus=1  mem=2G      time=5:00:00,10:00:00  throttle=16" "processdeps=none"
    add_debasher_process "parallel_split_tum_bam"    "cpus=1  mem=2G      time=5:00:00,10:00:00  throttle=16" "processdeps=none"
    add_debasher_process "parallel_samtools_mpileup_norm_bam" "cpus=1  mem=4G  time=8:00:00,16:00:00 throttle=24" "processdeps=afterok:create_genref_for_bam,aftercorr:parallel_split_norm_bam"
    add_debasher_process "parallel_samtools_mpileup_tum_bam"  "cpus=1  mem=4G  time=8:00:00,16:00:00 throttle=24" "processdeps=afterok:create_genref_for_bam,aftercorr:parallel_split_tum_bam"
    add_debasher_process "gen_sequenza_gcc"          "cpus=1  mem=1G      time=1:00:00"  "processdeps=afterok:create_genref_for_bam"
    add_debasher_process "parallel_bam2seqz"         "cpus=1  mem=2G      time=5:00:00  throttle=16" "processdeps=afterok:gen_sequenza_gcc,aftercorr:parallel_samtools_mpileup_norm_bam,aftercorr:parallel_samtools_mpileup_tum_bam"
    add_debasher_process "seqzmerge_plus_sequenza"   "cpus=1  mem=10G     time=5:00:00,10:00:00"  "processdeps=afterok:parallel_bam2seqz"
    add_debasher_process "parallel_delly"            "cpus=1  mem=10G,30G time=5:00:00  throttle=16" "processdeps=afterok:create_genref_for_bam,aftercorr:parallel_split_norm_bam,aftercorr:parallel_split_tum_bam"
    add_debasher_process "parallel_lumpy"            "cpus=2  mem=10G,30G time=5:00:00,10:00:00  throttle=16" "processdeps=aftercorr:parallel_split_norm_bam,aftercorr:parallel_split_tum_bam"
    add_debasher_process "parallel_svtyper"          "cpus=1  mem=8G,16G  time=6:00:00,24:00:00  throttle=16" "processdeps=aftercorr:parallel_lumpy"
    add_debasher_process "strelka_germline"          "cpus=8  mem=6G      time=6:00:00,12:00:00"  "processdeps=afterok:create_genref_for_bam"
    add_debasher_process "manta_somatic"             "cpus=8  mem=8G      time=8:00:00,16:00:00"  "processdeps=afterok:create_genref_for_bam"
    add_debasher_process "strelka_somatic"           "cpus=8  mem=6G      time=8:00:00,16:00:00"  "processdeps=afterok:create_genref_for_bam,afterok:manta_somatic"
    add_debasher_process "msisensor"                 "cpus=2  mem=10G,20G time=8:00:00,24:00:00"  "processdeps=afterok:create_genref_for_bam"
    add_debasher_process "snp_pileup_plus_facets"    "cpus=1  mem=25G     time=8:00:00,24:00:00"  "processdeps=afterok:create_genref_for_bam"
    add_debasher_process "cnvkit"                    "cpus=8  mem=30G     time=6:00:00,24:00:00"  "processdeps=afterok:create_genref_for_bam"
    add_debasher_process "platypus_germline"         "cpus=1  mem=4096    time=8:00:00,16:00:00"  "processdeps=afterok:create_genref_for_bam"
}
