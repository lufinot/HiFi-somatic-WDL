version 1.0


task trgt {
	input {
		File bam
        File bam_index

        String pname
        File ref_fasta
        File ref_fasta_index
		File tandem_repeat_bed

        Boolean is_tumor

        Int threads

		RuntimeAttributes runtime_attributes
	}

	Boolean sex_defined = defined(sex)
	String karyotype = if select_first([sex, "FEMALE"]) == "MALE" then "XY" else "XX"
	String bam_basename = basename(bam, ".bam")
	Int disk_size = ceil((size(bam, "GB") + size(ref_fasta, "GB")) * 2 + 20)
    String label = if is_tumor then "tumor" else "control"

	command <<<
		set -euo pipefail

		echo ~{if sex_defined then "" else "Sex is not defined for ~{sample_id}.  Defaulting to karyotype XX for TRGT."}

		trgt --version

		trgt \
			--threads ~{threads} \
			--karyotype ~{karyotype} \
			--genome ~{ref_fasta} \
			--repeats ~{tandem_repeat_bed} \
			--reads ~{bam} \
			--output-prefix ~{bam_basename}.trgt

		bcftools --version

		bcftools sort \
			--output-type z \
			--output ~{bam_basename}.trgt.sorted.vcf.gz \
			~{bam_basename}.trgt.vcf.gz

		bcftools index \
			--threads ~{threads - 1} \
			--tbi \
			~{bam_basename}.trgt.sorted.vcf.gz
		
		samtools --version

		samtools sort \
			-@ ~{threads - 1} \
			-o ~{bam_basename}.trgt.spanning.sorted.bam \
			~{bam_basename}.trgt.spanning.bam

		samtools index \
			-@ ~{threads - 1} \
			~{bam_basename}.trgt.spanning.sorted.bam
	>>>

	output {
		File spanning_reads = pname + "_trgrt/" + label + "~{bam_basename}.trgt.spanning.sorted.bam"
		File spanning_reads_index = pname + "_trgrt/" + label + "~{bam_basename}.trgt.spanning.sorted.bam.bai"
		File repeat_vcf = pname + "_trgrt/" + label + "~{bam_basename}.trgt.sorted.vcf.gz"
		File repeat_vcf_index = pname + "_trgrt/" + label + "~{bam_basename}.trgt.sorted.vcf.gz.tbi"
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

