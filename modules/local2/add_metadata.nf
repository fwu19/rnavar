process ADD_METADATA {
    container "docker://fwu19/r-libs:4.1.2" 

    label 'process_single'

    tag "Generate ${out_csv}"

    input:
    path ( in_csv, stageAs: "input.csv" )
    path ( metadata )
    val ( out_csv )

    output:
    path ( out_csv ), emit: csv
    path  ( "versions.yml" ), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    add_metadata.r input.csv $metadata $out_csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1)
    END_VERSIONS
    """
}
