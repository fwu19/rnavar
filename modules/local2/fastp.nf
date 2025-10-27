process FASTP {
    module = ['fastp/0.23.4-GCC-13.2.0']

    cpus = 6
    memory = 36G
    time = 4.h

    tag "FASTP on ${meta.id}"

    publishDir "${params.outdir}/trimmed_fastq/", pattern: '*.gz', mode: 'copy'

    input:
    tuple val(meta), path(reads, stageAs: "input*/*")
    

    output:
    tuple val(meta), path("trimmed_fastq/*"), emit: reads
    tuple val(meta), path( "*.fastp.json" ), emit: js
    tuple val(meta), path( "*.fastp.html" ), emit: html

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readList = reads instanceof List ? reads.collect { it.toString() } : [reads.toString()]
    def read1 = readList[0]
    def read2 = readList[1]
    //readList.eachWithIndex { v, ix -> (ix & 1 ? read2 : read1) << v }

    if (params.adapters){
        adapter_list = params.adapters.split(',').collect()
        if ( adapter_list.size == 1){
            adapter1 = adapter_list[0]
            adapter2 = adapter_list[0]
        } else {
            adapter1 = adapter_list[0]
            adapter2 = adapter_list[1]
        }

        """
        mkdir trimmed_fastq
        fastp -w ${task.cpus} \
        $args \
        --adapter_sequence $adapter1 --adapter_sequence_r2 $adapter2 \
        --in1 $read1 --in2 $read2 \
        --out1 trimmed_fastq/${prefix}_1.fastq.gz --out2 trimmed_fastq/${prefix}_2.fastq.gz \
        -j ${prefix}.fastp.json -h ${prefix}.fastp.html \
        --detect_adapter_for_pe -l 20 -g
        """

    }else{

        """
        mkdir trimmed_fastq
        fastp -w ${task.cpus} \
        $args \
        --in1 $read1 --in2 $read2 \
        --out1 trimmed_fastq/${prefix}_1.fastq.gz --out2 trimmed_fastq/${prefix}_2.fastq.gz \
        -j ${prefix}.fastp.json -h ${prefix}.fastp.html \
        --detect_adapter_for_pe -l 20 -g
        """
    }
}
