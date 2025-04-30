We first recommend to clone the repository.

=== "git"

    ``` bash
    git clone https://github.com/dienerlab/pipelines
    ```

=== "GitHub CLI"

    ``` bash
    gh repo clone dienerlab/pipelines
    ```

It's also a good idea to install nextflow into your conda base environment.

``` bash
conda install -c conda-forge -c bioconda nextflow
```

## Install environments

Most of the pipelines contain a conda environment file that can be used. For instance,
to set up the metagenomics environment:

``` bash
cd /path/to/cloned/repository
conda env create -f metagenomics/conda.yml
```

No additional setup is required for pipelines using docker or apptainer/singularity.

## Configure Nextflow

If you are part of our lab pleae use the configuration [from our wiki](https://github.com/dienerlab/internal/wiki/configs).

Otherwise you can use the following basic configuration which sets up profiles for local
execution or SLURM:

``` groovy title="nextflow.config"
--8<-- "nextflow.config"
```

Make sure to make any adjustments beforehand. For instance, for SLURM you probably want to use
the correct name for the queue in your HPC system.If you use Seqera Cloud (formerly Tower) also
add in your token and woskspace ID and set it to enabled.

Either copy this file into your project directory or create the file `~/.nextflow/config` in
your home directory and copy the contents there.

``` bash
mkdir ~/.nextflow
```

After that edit and copy the config:

``` bash
cp /path/to/pipelines/nextlow.config ~/.nextflow/config
```

You can switch between the profiles by using
the Nextflow `-profile` flag.

``` bash
nextflow run [PIPELINE] -profile slurm ...
```

!!! NOTE

    For the lab config the default is already using SLURM. You only need to specify the
    profile if you want to force local execution:

    ``` bash
    nextflow run [PIPELINE] -profile local ...
    ```

## Run the pipelines

There are generally two options to run the pipelines:

1. Directly call the pipeline from the cloned directory
2. Create a symlink to the specific `[pipeline].nf` file

=== "direct"

    ```bash
    export DLP=/path/to/cloned/repository
    nextflow run $DLP/metagenomics ...
    ```

=== "link"

    ```bash
    ln -s /path/to/cloned/repository/metagenomics/main.nf
    nextflow run main.nf ...
    ```

If the pipeline is using a conda environment it should be run using the `-with-conda` flag
indicating the directory where the environment was created, usually one in `~/miniforge3/envs`.
For the docker based ones usually it is run using apptainer/singularity.

Please see the individual pipelines for instructions.

For the metagenomics pipeline this looks like:

``` bash
nexflow run $DLP/metagenomics/main.nf -resume -with-conda=$HOME/miniforge3/envs/metagenomics ...
```