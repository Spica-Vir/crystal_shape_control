# LAMMPS job submitter - Imperial cluster version

[LAMMPS](https://www.lammps.org/) job submitter for [Imperial Cluster](https://www.imperial.ac.uk/admin-services/ict/self-service/research-support/rcs/), PBS job scheduler. 

## Install

1. Use `scp -r` to upload this folder to any sub-directory of `${HOME}` on the cluster.  
2. Enter the folder and execute the script `setup.sh`.  
3. Following the instructions, specify the directory of job submitter and the directory of executables, or module loading command. 
4. Type `source ~/.bashrc` to implement user-defined commands. 

**Note**

1. All the scripts should be placed in the same directory.  
2. By default, job submitter scripts will be stored in `${HOME}/runLAMMPS/`.  
3. By default, module 'lammps/19Mar2020' will be set as the executable. MPI: 'mpi/intel-2019' and 'intel-suite/2019.4'  

## Command list & Usage

### Commands and options

The script adopts the command-line options to launch parallel jobs, which are in the similar fashion to the native LAMMPS. 

Here are the user defined commands: 

1. `Plmp` - executing parallel LAMMPS calculations  

``` bash
Plmp -in <input> -wt <walltime> -nd <node> -- <opts>
```

`-in`       : str, main input file, must include '.in'  
`-wt`       : str, walltime, hh:mm format  
`-nd`       : int, number of nodes  
`-- <opts>` : str, optional, other LAMMPS options, in `-opt value` pairs  

The sequence of `-in` `-wt` `-nd` is arbitrary, but `-- <opts>` should always be placed at the end. The sequence within `<opts>` is arbitrary. 

Examples:

``` bash
> Plmp -in input.in -wt 12:00 -nd 2
> Plmp -in=input.in -wt=12:00 -nd=2 
> Plmp -in input.in -wt 12:00 -nd 2 -- -screen screen.dump
```

Submit files:

``` bash
> qsub input.qsub
```

2. `setlmp` - print the file `settings`. No input required.

### Output and job executing information

1. **\<jobname\>.o\<jobid\>**  
Information of file transferring and node synchronizing before the job is actually launched is printed out in this file. Scripts executed in the temporary directory are listed here. Check this file if computation / post-processing is not executed. 

2. **\<jobname\>.log**  
Message on screen, including basic information of the current job (header), mpi message (comments), LAMMPS output (body), and post-processing message (ending). Check this file for the progress of current job, or problems with post-processing.   

3. **\<jobname\>.out**  
LAMMPS output file, originally named as 'log.lammps'. Same as the 'body' section of '\<jobname\>.log'. 

If a '.out' or a '.log' file with the same name as the job to be submitted exists in the same directory, that job won't be executed before output is either transferred to another folder or removed. An error message will be given when generating the '.qsub' file. 

4. **\<jobname\>.restart/**  
A folder including all restart files. Restart information generated by LAMMPS using command `restart`. According to the manual, the name of restart file is arbitrary. To normalize the naming scheme, by default, names of restart files should include '\*.restart\*', i.e., at least named as '.restart'. The keyword identifying restart files can be changed in 'settings' file. 

5. **\<jobname\>.dump/**   
A folder including all dump files. Dump files are generated by LAMMPS using command `dump`. Their naming scheme is also arbitrary. By default, names of dump files should include '\*.dump\*'. The keyword identifying dump files can be changed in 'settings' file.

### Temporary directory

By default, the temporary directory is set as '/rds/general/ephemeral/user/${USER}/ephemeral', where the contents will be automatically removed after 30 days. The folder for current job is named as '\<jobname\>\_\<jobid\>'

* If the job is terminated due to exceeding walltime, temporary files will be saved in the output directory. The temporary directory will be removed.

* If the job is terminated due to improper settings of calculation parameters, temporary files will be saved in the output directory. The temporary directory will be removed.

* If the job is killed before 'timeout' (usually by user), temporary will be saved in the temporary directory with temporary names. The temporary directory will not be removed. Refer to '.log' file or '.o\<jobid\>' file for the path. 

## Script list

**setup.sh**          : set up the 'settings' file and create job submission commands  
**settings**          : store all parameters needed for LAMMPS jobs. see the 'Keyword list' below.  
**settings_template** : empty 'settings' file, will be used to cover 'settings' file when installing/re-installing the job submitter  
**gen_sub**           : generate submission file ('.qsub')  
**runPlmp**           : execute LAMMPS calculations in parallel   
**postlmp**           : Copy & save files from temporary directory to the output directory.  

**NOTE**

1. The name of file 'settings', 'gen_sub', 'settings_template' shouldn't be changed.  
2. Names of 'runPlmp', 'postlmp' can be changed, but should corresponds to the values in 'settings' - not recommended.  

## Keyword list
Keywords used for the script 'settings' are listed in the table below. Modify the values in the same file to change the parameters used during computation.

| KEYWORD                 | DEFAULT VALUE   | DEFINITION |
|:------------------------|:---------------:|:-----------|
| SUBMISSION_EXT          | .qsub           | extension of job submission script |
| NCPU_PER_NODE           | 24              | Number of processors per node |
| MEM_PER_NODE            | 50              | Unit: GB. Allocated memory per node |
| BUDGET_CODE             | -               | Budget code of a research project, for ARCHER2|
| QOS                     | -               | Quality of service, for ARCHER2 |
| PARTITION               | -               | Partition of jobs, for ARCHER2 |
| TIME_OUT                | 3               | Unit: min. Time spared for post processing |
| LMP_SCRIPT              | runPlmp         | Script for parallel LAMMPS |
| POST_PROCESSING_SCRIPT  | postlmp         | Post processing script |
| JOB_TMPDIR              | /rds/general/ephemeral/user/${USER}/ephemeral | Temporary directory for calculations |
| EXEDIR                  | module load  lammps/19Mar2020 | Directory of executables / Available module |
| EXE_PLMP                | lmp_mpi         | Executable for parallel LAMMPS |
| EXE_LMP                 | -               | Executable for serial LAMMPS, for workstation |
| PRE_CALC                | \[Table\]       | Saved names, temporary names, and definitions of input files |
| POST_CALC               | \[Table\]       | Saved names, temporary names, and definitions of output files |
| JOB_SUBMISSION_TEMPLATE | \[script\]      | Template for job submission files |

**NOTE**

1. Keyword `JOB_SUBMISSION_TEMPLATE` should be the last keyword, but the sequences of other keywords are allowed to change.  
2. Empty lines between keywords and their values are forbidden.  
3. All listed keywords have been included in the scripts. Undefined keywords are left blank.  
3. Multiple-line input for keywords other than `PRE_CALC`, `POST_CALC`, and `JOB_SUBMISSION_TEMPLATE` is forbidden.  
4. Dashed lines for `PRE_CALC`, `POST_CALC`, and `JOB_SUBMISSION_TEMPLATE` are used to define input blocks and are not allowed to be modified. Minimum length: '------------------'