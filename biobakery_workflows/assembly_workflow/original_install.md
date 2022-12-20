Install everything manually
```
conda create --name biobakery_assembly -c biobakery biobakery_workflows
conda activate biobakery_assembly
git clone --recursive https://github.com/chklovski/checkm2.git && cd checkm2 && conda env update --file checkm2.yml
conda activate biobakery_assembly
conda install -c bioconda megahit metabat2 assembly-stats fastani
conda install -c r r
conda install -c conda-forge python-leveldb r-stringi multiprocess
cd ..
checkm2/bin/checkm2 database --download --path /custom/checkm/path/
export CHECKM_DATA_PATH=/custom/checkm/path/
mkdir tmp_to_delete && mkdir -p /custom/phylophlan/path/ && touch tmp_to_delete/tmp.fa && phylophlan_metagenomic -d SGB.Jul20 --database_folder /custom/phylophlan/path/ -i tmp_to_delete/; rm -r tmp_to_delete/
export PHYLOPHLAN_PATH=/custom/phylophlan/path/
```

# Install necessary R packages
```
R
install.packages(c("docopt", "dplyr", "data.table", "stringr", "doParallel", "textshape"))
q()
```