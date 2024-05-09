version 1.0


# Call repeats using
task trgt {
    input {
        File bam
        File bam_index
        String pname
        String sex
        File ref_fasta
        File ref_fasta_index
        File tandem_repeat_bed

        Int threads
    }

    Boolean sex_defined = defined(sex)
    String karyotype = if select_first([sex, "FEMALE"]) == "MALE" then "XY" else "XX"
    String bam_basename = basename(bam, ".bam")
    Int disk_size = ceil((size(bam, "GB") + size(ref_fasta, "GB")) * 2 + 20)

    command <<<
        set -euxo pipefail

        echo "Running TRGT on ~{pname} "

        echo ~{if sex_defined then "" else "Sex is not defined for ~{pname}.  Defaulting to karyotype XX for TRGT."}

        trgt --version

        trgt \
            --threads ~{threads} \
            --karyotype ~{karyotype} \
            --genome ~{ref_fasta} \
            --repeats ~{tandem_repeat_bed} \
            --reads ~{bam} \
            --output-prefix ~{pname}_trgt/~{bam_basename}.trgt

        bcftools --version

        bcftools sort \
            --output-type z \
            --output ~{pname}_trgt/~{bam_basename}.trgt.sorted.vcf.gz \
            ~{pname}_trgt/~{bam_basename}.trgt.vcf.gz

        bcftools index \
            --threads ~{threads - 1} \
            --tbi \
            ~{pname}_trgt/~{bam_basename}.trgt.sorted.vcf.gz
        
        samtools --version

        samtools sort \
            -@ ~{threads - 1} \
            -o ~{bam_basename}.trgt.spanning.sorted.bam \
            ~{pname}_trgt/~{bam_basename}.trgt.spanning.bam

        samtools index \
            -@ ~{threads - 1} \
            ~{pname}_trgt/~{bam_basename}.trgt.spanning.sorted.bam
    >>>

    output {
        File spanning_reads = "~{pname}_trgt/~{bam_basename}.trgt.spanning.sorted.bam"
        File spanning_reads_index = "~{pname}_trgt/~{bam_basename}.trgt.spanning.sorted.bam.bai"
        File repeat_vcf = "~{pname}_trgt/~{bam_basename}.trgt.sorted.vcf.gz"
        File repeat_vcf_index = "~{pname}_trgrt/~{bam_basename}.trgt.sorted.vcf.gz.tbi"
    }

    runtime {
        docker: "~quay.io/pacbio/trgt@sha256:8c9f236eb3422e79d7843ffd59e1cbd9b76774525f20d88cd68ca64eb63054eb"
        cpu: threads
        memory: "4 GB"
        disk: disk_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
        maxRetries: 2
    }
}

