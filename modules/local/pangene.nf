/*
 * PANGENE: one pangenome for the whole run.
 *
 * pangene builds its graph from one set of proteins aligned to every genome.
 * The set is the reference proteome and the genes of each sample that no
 * reference gene matches (bin/pangenome.py); miniprot aligns it to every
 * sample's genome, the sample's own included, one alignment per genome. Each
 * gene is then one row of the matrix and each genome one column, named after
 * its sample, so the genes the samples share are the rows found in every
 * column.
 *
 * The genomes are <name>.fasta, one per sample, as PATHOTYPR writes them or,
 * with --skip_typing, PANGENE_GENOME below. Miniprot jobs run in parallel
 * batches, bounded by task.cpus to avoid oversubscription.
 */
process PANGENE {
    tag "${genomes.size()} genomes"
    // The AMAP code this task runs, in its cache key (see docs/about/contributing.md)
    ext code: { files("${projectDir}/bin/pangenome.py") }

    publishDir "${params.outdir}/PANGENE", mode: 'link', overwrite: true,
               pattern: "pangenome.{gfa,bubbles.txt,Rtab,stat}"

    input:
    // Every sample's genome, <name>.fasta: a list even for one sample
    path genomes, stageAs: 'genomes/*', arity: '1..*'
    // Every sample's <name>.AMAP.tsv and <name>.AMAP.faa, for the genes no
    // reference gene matches
    path corrected, stageAs: 'corrected/*', arity: '1..*'
    // In a folder of its own: a proteome on the host may have any name,
    // proteins.faa included, which this task writes
    path ref_proteome, stageAs: 'reference/*'

    output:
    path "pangenome.Rtab", emit: rtab
    path "pangenome.{gfa,bubbles.txt,stat}"

    script:
    """
    # AMAP code in this task's cache key: ${task.ext.code.collect { f -> f.name }.join(' ')}
    set -euo pipefail

    # 1) The proteins every genome is searched with
    python3 ${projectDir}/bin/pangenome.py proteins --reference ${ref_proteome} \\
        --out proteins.faa corrected/*.AMAP.tsv

    # 2) Align them to every genome, the sample's own included, in parallel
    #    with bounded concurrency
    genomes=(genomes/*.fasta)
    echo "==> aligning \$(grep -c '>' proteins.faa) proteins to \${#genomes[@]} genomes (${task.cpus} CPUs)..."
    THREADS_PER_JOB=2
    MAX_JOBS=\$(( ${task.cpus} / THREADS_PER_JOB ))
    [ \$MAX_JOBS -lt 1 ] && MAX_JOBS=1

    mkdir -p paf
    pids=()
    fail=0
    for genome in "\${genomes[@]}"; do
        miniprot --outs=0.97 --no-cs -S -Iut\${THREADS_PER_JOB} \\
                 "\$genome" proteins.faa \\
                 > "paf/\$(basename "\$genome" .fasta).paf" &
        pids+=(\$!)
        if [ \${#pids[@]} -ge \$MAX_JOBS ]; then
            for p in "\${pids[@]}"; do
                wait \$p || { echo "ERROR: miniprot PID \$p failed" >&2; fail=1; }
            done
            pids=()
        fi
    done
    # Guarded: an empty array is unbound under set -u in the bash before 4.4
    # that macOS has as /bin/bash
    [ \${#pids[@]} -eq 0 ] || for p in "\${pids[@]}"; do
        wait \$p || { echo "ERROR: miniprot PID \$p failed" >&2; fail=1; }
    done
    [ \$fail -eq 0 ] || exit 1

    # 3) One graph for the whole run; pangene names each genome after its PAF
    pangene -p0 paf/*.paf > pangenome.gfa
    k8 /usr/local/bin/pangene.js call       pangenome.gfa > pangenome.bubbles.txt
    k8 /usr/local/bin/pangene.js gfa2matrix pangenome.gfa > genes_by_genome.Rtab
    gfatools stat pangenome.gfa > pangenome.stat

    # 4) One column per genome, named after its sample, a genome none of the
    #    proteins was found in included
    python3 ${projectDir}/bin/pangenome.py matrix --out pangenome.Rtab \\
        genes_by_genome.Rtab "\${genomes[@]}"

    rm -rf paf
    """
}

/*
 * The genome PANGENE aligns against when lineage typing is skipped: the
 * sample sheet's own, copied to <sampleId>.fasta. As given, two samples' files
 * can share a name (assembly.fasta from two assemblies), which PANGENE cannot
 * stage side by side, and a name that does not end in .fasta is never aligned
 * against.
 */
process PANGENE_GENOME {
    tag "$sampleId"

    input:
    // Staged in a directory of its own: a genome already called
    // <sampleId>.fasta would otherwise be copied onto itself, through the link
    // Nextflow stages it as.
    tuple val(sampleId), path(fasta, stageAs: 'input/*')

    output:
    tuple val(sampleId), path("${sampleId}.fasta")

    script:
    """
    cp "${fasta}" "${sampleId}.fasta"
    """
}
