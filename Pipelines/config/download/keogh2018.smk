module keogh2018:
    snakefile: "_keogh2018.smk"
    config: config
use rule * from keogh2018 as keogh2018_*