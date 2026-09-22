#!/usr/bin/env nextflow
process PREPARE_EGGNOG {
    label 'eggnog_prepare'
    tag 'Prepare eggNOG database'

    output:
    path "eggnog_db", emit: database_dir

    script:
    """
    set -euo pipefail

    mkdir -p eggnog_db

    export MAMBA_SKIP_ACTIVATE=""
    source /usr/local/bin/_activate_current_env.sh

    DOWNLOADER=\$(command -v download_eggnog_data.py || true)

    if [ -z "\$DOWNLOADER" ]; then
        echo "ERROR: download_eggnog_data.py was not found inside the eggNOG container." >&2
        exit 1
    fi

    cp "\$DOWNLOADER" ./download_eggnog_data.py

    # eggNOG-mapper 2.1.13 still points to the obsolete eggnogdb.embl.de host.
    # Keep the original mapper version, but update only the download host.
    sed -i \
        's#eggnogdb.embl.de#eggnog5.embl.de#g' \
        download_eggnog_data.py

    python download_eggnog_data.py \
        -y \
        -M \
        --data_dir eggnog_db

    # The original downloader may exit successfully even when wget fails,
    # therefore validate the required database files explicitly.
    echo "Validating eggNOG database..."

    test -s eggnog_db/eggnog.db || {
	echo "ERROR: eggnog.db is missing." >&2
	exit 1
    }

    test -s eggnog_db/eggnog.taxa.db || {
	echo "ERROR: eggnog.taxa.db is missing." >&2
	exit 1
    }

    test -s eggnog_db/eggnog_proteins.dmnd || {
	echo "ERROR: eggnog_proteins.dmnd is missing." >&2
	exit 1
    }

    test -s eggnog_db/mmseqs/mmseqs.db || {
	echo "ERROR: mmseqs/mmseqs.db is missing." >&2
	exit 1
    }

     if ! find eggnog_db -maxdepth 1 -iname '*mmseqs*' -print -quit | grep -q .; then
        echo "ERROR: eggNOG MMseqs2 database was not downloaded." >&2
        exit 1
     fi
     """
   }
}


process EGGNOG_PROKARYOTE {
    label "eggnog_mapper_v2"
    tag "${fasta_file.baseName}"
    publishDir "data/reproduced/prokaryote_output_tools/eggnog", mode: 'copy'

    input:
    path fasta_file
    path eggnog_db

    output:
    path("${fasta_file.baseName}_eggnog/*"), emit: eggnog_results

    script:
    """
    set -euo pipefail

    export MAMBA_SKIP_ACTIVATE=""
    export EGGNOG_DATA_DIR=/database
    export NXT_TASK_MONITOR=0

    source /usr/local/bin/_activate_current_env.sh

    # Use a local temp directory to avoid polluting the work directory
    mkdir -p tmp
    mkdir -p ${fasta_file.baseName}_eggnog
    emapper.py \
    --itype genome \
    --genepred prodigal \
    --decorate_gff yes \
    -i ${fasta_file} \
    -o ${fasta_file.baseName} \
    --data_dir /database \
    -m mmseqs \
    --cpu ${task.cpus} \
    --output_dir ${fasta_file.baseName}_eggnog \
    --dbmem
    """

}

process EGGNOG_EUKARYOTE {
    label "eggnog_mapper_v2"
    tag "$sample_id"
    publishDir "data/reproduced/eukaryote_output_tools/eggnog", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)
    path eggnog_db

    output:
    tuple val(sample_id), path ("${sample_id}_eggnog/*"), emit: eggnog_results

    script:
    """
    set -euo pipefail
    # Activate the container's environment for eggNOG-mapper.
    export MAMBA_SKIP_ACTIVATE=""
    export EGGNOG_DATA_DIR=/database
    export NXT_TASK_MONITOR=0

    # Use a local temp directory to avoid polluting the work directory
    mkdir -p tmp
    mkdir -p ${sample_id}_eggnog

    source /usr/local/bin/_activate_current_env.sh

    emapper.py \
        -i ${proteins} \
        --itype proteins \
        -m mmseqs \
        --cpu ${task.cpus} \
        --data_dir /database \
	--output ${sample_id} \
	--output_dir ${sample_id}_eggnog \
        --temp_dir ./tmp
        """
}
