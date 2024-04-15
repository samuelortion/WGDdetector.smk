def log_file_name(rule_name, ext="stdout"):
    return (
        lambda wildcards: logdir
        / rule_name
        / f"{os.path.basename(wildcards.name)}.{ext}"
    )
