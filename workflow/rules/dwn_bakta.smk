rule download_bakta_db:
    output:
        directory(config['workdir'] + "/db_dir")
    conda:
        "../envs/bakta.yaml"
    run:
        if not config.get('db_dir'):
            shell("bakta_db download --output {output} --type {config['db_type']}")