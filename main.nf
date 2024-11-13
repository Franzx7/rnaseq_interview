/* Author: Franz AKE */
/* Nextflow script for RNA-Seq analysis for Interview */


// required params
params.accession = null
params.genome_gtf = null
params.genome_fasta = null


// checks
if (params.accession==null)
        error( "Provide acession IDs... [--accession]")

if (params.genome_fasta==null)
        error( "Provide org genome FASTA reference... [--genome_fasta]")

if (params.genome_gtf==null)
        error( "Provide org GTF reference... [--genome_gtf]")



process fetch_samples {
        tag "${sample}"
        conda 'bioconda::sra-tools=3.1.1'
        input:
                val sample
        output:
                file "*.fastq"
        script:
        """
                 fastq-dump --split-files ${sample}
        """
}


process qc_fastqc {
        tag "${fastqs}"
        conda 'bioconda::fastqc'
        publishDir './qc_reports/', overwrite: true, mode: 'copy'
        input:
            file fastqs
        output:
                file "*_fastqc*"
        script:
        """
                fastqc ${fastqs}

        """
}


process generate_genome_index_STAR {
        tag "${genome_fasta}, ${genome_gtf}"
        conda 'bioconda::star'
        input:
            path genome_fasta
            path genome_gtf
        output:
            file "genome_index"
        script:
        """
        mkdir genome_index/
        STAR --runThreadN 6 --runMode genomeGenerate --genomeDir genome_index --genomeFastaFiles ${genome_fasta} --sjdbGTFfile ${genome_gtf} --sjdbOverhang 99
        """
}


process alignment_to_genome {
        tag "${genome_index}, ${fastqs}"
        publishDir "./STAR_aligned_bam_files", overwrite: true, mode: 'copy'
        conda 'bioconda::star'
        input:
                path genome_index
                path fastqs
        output:
            file "*.out*"
        script:
        """
        STAR --genomeDir ${genome_index} --runThreadN 50 --readFilesIn ${fastqs} --outFileNamePrefix ${fastqs[1]}_mapped --outSAMtype BAM SortedByCoordinate --outSAMunmapped With>
        """
}



process generate_countMat {
        tag "${sample_bams}, ${genome_gtf}"
        publishDir "./final_counts", overwrite: true, mode: 'copy'
        conda 'bioconda::subread'
        input:
            file sample_bams
            path genome_gtf
        output:
            file 'counts.txt'
        script:
        """
        featureCounts -O -t gene -g gene_name -a ${genome_gtf} -o counts.txt ./*.bam -T 20 -p
        """
}



// Main
workflow {

        fetch_samples( Channel.of(params.accession).splitCsv().flatten() )
        qc_fastqc( fetch_samples.out )
        generate_genome_index_STAR( params.genome_fasta, params.genome_gtf )
        alignment_to_genome( generate_genome_index_STAR.out, fetch_samples.out )
        generate_countMat( alignment_to_genome.out.collect(), params.genome_gtf )
}

