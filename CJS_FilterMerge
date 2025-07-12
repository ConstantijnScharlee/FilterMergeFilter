version 1.0

workflow FilterAndMergeVCFs {
  input {
    File sample_list_tsv
    String tissue_name
    String gcs_input_prefix  # e.g. gs://bucket/samples
    String gcs_output_prefix # e.g. gs://bucket/outputtedsamples

    Int extract_cpu
    String extract_memory

    Int fix_cpu
    String fix_memory

    Int merge_cpu
    String merge_memory
  }

  Array[String] sample_ids = read_lines(sample_list_tsv)

  scatter (sample_id in sample_ids) {
    call DownloadVCF {
      input:
        sample_id = sample_id,
        gcs_input_prefix = gcs_input_prefix
    }

    call FilterVCF {
      input:
        input_vcf = DownloadVCF.vcf,
        sample_id = sample_id,
        cpu = extract_cpu,
        memory = extract_memory
    }

    call FixPloidyAndBgzip {
      input:
        filtered_vcf = FilterVCF.filtered_vcf,
        cpu = fix_cpu,
        memory = fix_memory
    }

    call UploadToGCS {
      input:
        file_to_upload = FixPloidyAndBgzip.bgzipped_vcf,
        sample_id = sample_id,
        gcs_output_prefix = gcs_output_prefix
    }
  }

  call MergeAndIndexVCFs {
    input:
      vcfs_to_merge = FixPloidyAndBgzip.bgzipped_vcf,
      output_name = tissue_name + "_merged.vcf.gz",
      cpu = merge_cpu,
      memory = merge_memory
  }

  output {
    File merged_vcf = MergeAndIndexVCFs.merged_vcf
    File merged_vcf_index = MergeAndIndexVCFs.merged_vcf_index
  }
}

task DownloadVCF {
  input {
    String sample_id
    String gcs_input_prefix
  }

  command <<<
    gsutil cp ${gcs_input_prefix}/${sample_id}/${sample_id}.hard-filtered.vcf.gz .
    gsutil cp ${gcs_input_prefix}/${sample_id}/${sample_id}.hard-filtered.vcf.gz.tbi .
  >>>

  output {
    File vcf = "${sample_id}.hard-filtered.vcf.gz"
    File index = "${sample_id}.hard-filtered.vcf.gz.tbi"
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
    memory: "1G"
    cpu: 1
  }
}

task FilterVCF {
  input {
    File input_vcf
    String sample_id
    Int cpu
    String memory
  }

  command <<<
    bcftools view \
      -i 'FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) > 0.35 && FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) < 0.65' \
      ${input_vcf} -Oz -o ${sample_id}_HET.vcf.gz
    tabix -p vcf ${sample_id}_HET.vcf.gz
  >>>

  output {
    File filtered_vcf = "${sample_id}_HET.vcf.gz"
  }

  runtime {
    docker: "quay.io/biocontainers/bcftools:1.20--h7e0c702_0"
    memory: memory
    cpu: cpu
  }
}

task FixPloidyAndBgzip {
  input {
    File filtered_vcf
    Int cpu
    String memory
  }

  command <<<
    bcftools +fixploidy ${filtered_vcf} -o fixed.vcf
    bgzip -c fixed.vcf > fixed_HET_FP.vcf.gz
    tabix -p vcf fixed_HET_FP.vcf.gz
  >>>

  output {
    File bgzipped_vcf = "fixed_HET_FP.vcf.gz"
  }

  runtime {
    docker: "quay.io/biocontainers/bcftools:1.20--h7e0c702_0"
    memory: memory
    cpu: cpu
  }
}

task UploadToGCS {
  input {
    File file_to_upload
    String sample_id
    String gcs_output_prefix
  }

  command <<<
    gsutil cp ${file_to_upload} ${gcs_output_prefix}/${sample_id}_HET_FP.vcf.gz
    gsutil cp ${file_to_upload}.tbi ${gcs_output_prefix}/${sample_id}_HET_FP.vcf.gz.tbi
  >>>

  output {
    File uploaded_vcf = "${file_to_upload}"
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
    memory: "1G"
    cpu: 1
  }
}

task MergeAndIndexVCFs {
  input {
    Array[File] vcfs_to_merge
    String output_name
    Int cpu
    String memory
  }

  command <<<
    bcftools merge -Oz -o ${output_name} ${sep=' ' vcfs_to_merge}
    tabix -p vcf ${output_name}
  >>>

  output {
    File merged_vcf = output_name
    File merged_vcf_index = output_name + ".tbi"
  }

  runtime {
    docker: "quay.io/biocontainers/bcftools:1.20--h7e0c702_0"
    memory: memory
    cpu: cpu
  }
}
