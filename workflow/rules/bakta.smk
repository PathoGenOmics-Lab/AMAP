rule run_bakta:
    input:
        fasta = "/home/ruizro/bakta/sequences/{sample}.fasta",
        proteins = "sequence.gb"
    output:
        directory("/home/ruizro/bakta/{sample}/")
    conda:
        "../envs/bakta.yaml"
    params:
        db = "/home/ruizro/bakta/db/",
        prefix = "{sample}",
        locus_tag = "{sample}"
    shell:
        """
        bakta --db {params.db} --verbose --output {output} --prefix {params.prefix} --locus-tag {params.locus_tag} {input.fasta} --proteins {input.proteins}
        """