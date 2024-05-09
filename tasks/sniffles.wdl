version 1.0

# Use Sniffles to call somatic SVs
task sniffles {
	input {
		File bam
        File bam_index

        String pname
        File ref_fasta
        File ref_fasta_index

        Int threads
	}

	String bam_basename = basename(bam, ".bam")
	Int disk_size = ceil((size(bam, "GB") + size(ref_fasta, "GB")) * 2 + 20)

	command <<<
		set -euxo pipefail

		sniffles --version

		sniffles \ 
			--input ~{bam}\
			--vcf ~{pname}_sniffles/~{bam_basename}.sniffles.vcf \
			--reference ~{ref_fasta} \
			--mosaic \
			--threads ~{threads - 1}
	>>>

	output {
		File sv_vcf = "~{pname}_sniffles/~{bam_basename}.sniffles.vcf"
	}

	runtime {
		docker: "~quay.io/biocontainers/sniffles@sha256:feb1c41eae608ebc2c2cb594144bb3c221b87d9bb691d0af6ad06b49fd54573a"
		cpu: threads
		memory: "4 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: 1
		maxRetries: 2
	}
}

