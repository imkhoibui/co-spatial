process CELLRANGER_MKREF {
    tag "${params.meta}"
    label "process_high"

    container "quay.io/cumulus/cellranger-arc:2.0.0"

    output:
    path("${params.meta}/*")                 , emit: build

    script:
    def config_in           = task.ext.config_in ?: ""
    def version             = task.ext.version ?: ""
    def cratac_txt          = task.ext.cratac_txt ?: ""
    """
    #!/bin/bash
    cellranger-arc mkref --ref-version="$version" \
        --config="$config_in"
    """
}