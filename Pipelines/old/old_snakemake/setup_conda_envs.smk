#setup one conda environment given an "env" file listing the packages
#touch an output file to flag that env has been created
localrules: setup_env
rule setup_env:
    input:
        os.path.join(conda_env_path,"{envname}.env")
    output:
        "config/conda_envs/{envname}.touch"
    resources:
        conda_lock=1 #prevent multiple conda installs happening at once
    shell:
        "conda create -y -n {os.environ[MUT_PREFIX]}{wildcards.envname} "
        "    $(cat {os.environ[MUT_DIR]}/Pipelines/conda_envs/channels) "
        "    --file {input} && "
        "touch {output}"
