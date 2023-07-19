rule download_bakta_db:
    output:
        directory(config['workdir'] + "/db_dir")
    conda:
        "../envs/bakta.yaml"
    shell:
        """
        if not config.get('db_dir'):
            bakta_db download --output {output} --type {config['db_type']}
        """