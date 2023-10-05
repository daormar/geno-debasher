# *- bash -*

load_panpipe_module "bam_analysis"

copy_extended_pipeline()
{
    add_panpipe_process "copy_norm_bam"               "cpus=1  mem=2048    time=24:00:00" "processdeps=none"
    add_panpipe_process "copy_tum_bam"                "cpus=1  mem=2048    time=24:00:00" "processdeps=none"
    add_panpipe_process "index_norm_bam"              "cpus=1  mem=1024    time=4:00:00"  "processdeps=afterok:copy_norm_bam"
    add_panpipe_process "index_tum_bam"               "cpus=1  mem=1024    time=4:00:00"  "processdeps=afterok:copy_tum_bam"
    add_panpipe_process "bedtools_genomecov_norm_bam" "cpus=1  mem=2G      time=5:00:00" "processdeps=afterok:index_norm_bam"
    add_panpipe_process "bedtools_genomecov_tum_bam"  "cpus=1  mem=2G      time=5:00:00"  "processdeps=afterok:index_tum_bam"
    add_panpipe_process "create_genref_for_bam"       "cpus=1  mem=8G      time=4:00:00"  "processdeps=afterok:index_norm_bam"
    add_panpipe_process "parallel_split_norm_bam"     "cpus=1  mem=2G      time=5:00:00,10:00:00  throttle=16" "processdeps=afterok:index_norm_bam"
    add_panpipe_process "parallel_split_tum_bam"      "cpus=1  mem=2G      time=5:00:00,10:00:00  throttle=16" "processdeps=afterok:index_tum_bam"
    add_panpipe_process "parallel_samtools_mpileup_norm_bam" "cpus=1  mem=4G  time=8:00:00,16:00:00 throttle=24" "processdeps=afterok:create_genref_for_bam,aftercorr:parallel_split_norm_bam"
    add_panpipe_process "parallel_samtools_mpileup_tum_bam"  "cpus=1  mem=4G  time=8:00:00,16:00:00 throttle=24" "processdeps=afterok:create_genref_for_bam,aftercorr:parallel_split_tum_bam"
    add_panpipe_process "gen_sequenza_gcc"            "cpus=1  mem=1G      time=01:00:00" "processdeps=afterok:create_genref_for_bam"
    add_panpipe_process "parallel_bam2seqz"           "cpus=1  mem=2G      time=5:00:00  throttle=16" "processdeps=afterok:gen_sequenza_gcc,aftercorr:parallel_samtools_mpileup_norm_bam,aftercorr:parallel_samtools_mpileup_tum_bam"
    add_panpipe_process "seqzmerge_plus_sequenza"     "cpus=1  mem=10G     time=5:00:00,10:00:00"  "processdeps=afterok:parallel_bam2seqz"
    add_panpipe_process "parallel_delly"              "cpus=1  mem=10G,30G time=5:00:00  throttle=16" "processdeps=afterok:create_genref_for_bam,aftercorr:parallel_split_norm_bam,aftercorr:parallel_split_tum_bam"
    add_panpipe_process "parallel_lumpy"              "cpus=2  mem=10G,30G time=5:00:00,10:00:00  throttle=16" "processdeps=aftercorr:parallel_split_norm_bam,aftercorr:parallel_split_tum_bam"
    add_panpipe_process "parallel_svtyper"            "cpus=1  mem=8G,16G  time=6:00:00,24:00:00  throttle=16" "processdeps=aftercorr:parallel_lumpy"
    add_panpipe_process "strelka_germline"            "cpus=8  mem=6G      time=6:00:00,12:00:00"  "processdeps=afterok:create_genref_for_bam,afterok:index_norm_bam"
    add_panpipe_process "manta_somatic"               "cpus=8  mem=8G      time=8:00:00,16:00:00"  "processdeps=afterok:create_genref_for_bam,afterok:index_norm_bam,afterok:index_tum_bam"
    add_panpipe_process "strelka_somatic"             "cpus=8  mem=6G      time=8:00:00,16:00:00"  "processdeps=afterok:create_genref_for_bam,afterok:index_norm_bam,afterok:index_tum_bam,afterok:manta_somatic"
    add_panpipe_process "msisensor"                   "cpus=2  mem=10G,20G time=8:00:00,24:00:00"  "processdeps=afterok:create_genref_for_bam,afterok:index_norm_bam,afterok:index_tum_bam"
    add_panpipe_process "snp_pileup_plus_facets"      "cpus=1  mem=25G     time=8:00:00,24:00:00"  "processdeps=afterok:create_genref_for_bam,afterok:index_norm_bam,afterok:index_tum_bam"
    add_panpipe_process "cnvkit"                      "cpus=8  mem=30G     time=6:00:00,24:00:00"  "processdeps=afterok:create_genref_for_bam,afterok:index_norm_bam,afterok:index_tum_bam"
    add_panpipe_process "platypus_germline"           "cpus=1  mem=4096    time=8:00:00,16:00:00"  "processdeps=afterok:create_genref_for_bam,afterok:index_norm_bam"
    add_panpipe_process "clear_datadir"               "cpus=1  mem=1024    time=0:10:00"  "processdeps=afterok:parallel_delly,afterok:parallel_lumpy,afterok:parallel_svtyper,afterok:manta_somatic,afterok:strelka_germline,afterok:strelka_somatic,afterok:msisensor,afterok:cnvkit,afterok:snp_pileup_plus_facets,afterok:platypus_germline,afterok:seqzmerge_plus_sequenza"
}
