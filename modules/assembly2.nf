process assembly2 {
    input:
        val mypath
    output:
        //stdout
        val mypath
        //path "pyoutputs.txt", emit: pyoutputs
        
    
    """ 
    samplename=\$(echo ${mypath} | rev | cut -d "/" -f 1 | rev)
    #Call variants
    mkdir ${mypath}/variants
    singularity exec docker://staphb/samtools:1.12 samtools mpileup -A -d 8000 --reference ${params.reference}/NC_001474.2_DENV2.fasta -B -Q 0 ${mypath}/alignment/\${samplename}.primertrim.sorted.bam | ivar variants -r ${params.reference}/NC_001474.2_DENV2.fasta -m 10 -p ${mypath}/variants/\${samplename}.variants -q 20 -t 0.25 -g ${params.reference}/Dengue2_GCF_000871845.1_ViralProj20183_genomic.gff.gz

    #Generate consensus assembly
    mkdir ${mypath}/assembly
    singularity exec docker://staphb/samtools:1.12 samtools mpileup -A -B -d 8000 --reference ${params.reference}/NC_001474.2_DENV2.fasta -Q 0 ${mypath}/alignment/\${samplename}.primertrim.sorted.bam | ivar consensus -t 0 -m 10 -n N -p ${mypath}/assembly/\${samplename}.consensus
    """
}