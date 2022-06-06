# Running the pipeline
The scripts here were written to automate the running of the `MAKER` and `BRAKER` pipelines in the University of Florida High Powered Computing System (UFRC). UFRC uses the `slurm` job management system. These scripts are written to create standardized directory structures, configuration files, and job files for each step of the pipeline that can be submitted to `slurm`. Job files and shell scripts are labeled in the order in which they are supposed to be run.

All scripts require `python 3` to run. They depend on a configuration file, `template_pl_vars.txt`, which allows the user to control the names of files, the reference sequences used, the location of configuration files for programs like `AUGUSTUS` and `BUSCO`, and the versions of programs used. All scripts have been written with the UFRC's modular program paradigm in mind and are not guaranteed to work on other computer systems as-is. I have had to work around dependency clashes in several programs, so the pipeline takes several steps which may not be necessary on other systems. 

## Set-up
The file `template_pl_vars.txt` should be edited with relevant values for the run, including the genome assembly, run name, and references to be used. Values which I did not change across runs, such as the reference proteomes used and versions/locations of programs, have been pre-filled in the provided template file.

Programs like `BUSCO` and `AUGUSTUS` need a local writeable copy of their configuration directory. In UFRC, you can copy the config directory over like so:

```bash
module load augustus/3.4.0
CONFIG_DIR_LOC="/blue/soltis/kasey.pham/nitfix/augustus_config"
cp -r /apps/augustus/3.4.0/config "$CONFIG_DIR_LOC"
```

This local configuration directory location should be provided in `template_pl_vars.txt` as `AUGUSTUS_CONFIG`. You will also need to make a separate directory for each set of annotations if you are doing multiple and specify the location in `template_pl_vars.txt` as `BASE_DIR`.

Keep the files `parse_snap_errors.py` and `config.ini` in the directory specified by `PIPELINE_DIR` in `template_pl_vars.txt`.

## Usage
After filling out values in `template_pl_vars.txt`, scripts can be run as follows:

```bash
module load python/3.8

python preprocessing.py template_pl_vars.txt
python maker.py template_pl_vars.txt
python braker.py template_pl_vars.txt
```

After doing this, you should see a set of subdirectories, `repeat_lib`, `maker`, and `braker`, in the directory specified  in `template_pl_vars.txt`. Jobs in `repeat_lib` must be run first to mask the genome assembly for repeats before annotation, but `maker` and `braker` pipelines can be run in any order or concurrently.