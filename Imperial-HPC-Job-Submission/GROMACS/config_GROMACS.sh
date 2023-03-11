#!/bin/bash

function welcome_msg {
    core_version=`grep 'core' ${CTRLDIR}/version_control.txt | awk '{printf("%s", substr($0,22,11))}' | awk '{sub(/^ */, ""); sub(/ *$/, "")}1'`
    core_date=`grep 'core' ${CTRLDIR}/version_control.txt | awk '{printf("%s", substr($0,33,21))}' | awk '{sub(/^ */, ""); sub(/ *$/, "")}1'`
    core_author=`grep 'core' ${CTRLDIR}/version_control.txt | awk '{printf("%s", substr($0,54,21))}' | awk '{sub(/^ */, ""); sub(/ *$/, "")}1'`
    core_contact=`grep 'core' ${CTRLDIR}/version_control.txt | awk '{printf("%s", substr($0,75,31))}' | awk '{sub(/^ */, ""); sub(/ *$/, "")}1'`
    core_acknolg=`grep 'core' ${CTRLDIR}/version_control.txt | awk '{printf("%s", substr($0,106,length($0)))}' | awk '{sub(/^ */, ""); sub(/ *$/, "")}1'`
    code_version=`grep 'GROMACS' ${CTRLDIR}/version_control.txt | awk '{printf("%s", substr($0,22,11))}' | awk '{sub(/^ */, ""); sub(/ *$/, "")}1'`
    code_date=`grep 'GROMACS' ${CTRLDIR}/version_control.txt | awk '{printf("%s", substr($0,33,21))}' | awk '{sub(/^ */, ""); sub(/ *$/, "")}1'`
    code_author=`grep 'GROMACS' ${CTRLDIR}/version_control.txt | awk '{printf("%s", substr($0,54,21))}' | awk '{sub(/^ */, ""); sub(/ *$/, "")}1'`
    code_contact=`grep 'GROMACS' ${CTRLDIR}/version_control.txt | awk '{printf("%s", substr($0,75,31))}' | awk '{sub(/^ */, ""); sub(/ *$/, "")}1'`
    code_acknolg=`grep 'GROMACS' ${CTRLDIR}/version_control.txt | awk '{printf("%s", substr($0,106,length($0)))}' | awk '{sub(/^ */, ""); sub(/ *$/, "")}1'`
    cat << EOF
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
    _   _   _   _______   _          ___        ___       _______    _______   
   | | | | | | |  _____| | |       / ___ \    / ___ \    / _   _ \  |  _____|  
   | | | | | | | |       | |      / /   \_\  / /   \ \  | / \ / \ | | |        
   | | | | | | | |____   | |     | |        | |     | | | | | | | | | |____    
   | | | | | | | |____|  | |     | |        | |     | | | | | | | | | |____|   
   | | | | | | | |       | |     | |     __ | |     | | | | | | | | | |        
   | \_/ \_/ | | |_____  | |_____ \ \___/ /  \ \___/ /  | | |_| | | | |_____   
    \__/\___/  |_______| |_______| \ ___ /    \ ___ /   |_|     |_| |_______|  

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
GROMACS job submission script for Imperial HPC - Setting up

Job submission script installed date : `date`
Batch system                         : PBS
Job submission script version        : ${code_version} (${code_date})
Job submission script author         : ${code_author} (${code_contact})
Core script version                  : ${core_version} (${core_date})
Job submission script author         : ${core_author} (${core_contact})

${code_acknolg}
${core_acknolg}

EOF
}


function get_scriptdir {
    cat << EOF
================================================================================
    Note: all scripts should be placed into the same directory!
    Please specify your installation path.

    Default Option:
    ${HOME}/etc/runGROMACS/

EOF

    read -p " " SCRIPTDIR
    SCRIPTDIR=`realpath $(echo ${SCRIPTDIR})`

    if [[ -z ${SCRIPTDIR} ]]; then
        SCRIPTDIR=${HOME}/etc/runGROMACS/
    fi

    if [[ ${SCRIPTDIR: -1} == '/' ]]; then
        SCRIPTDIR=${SCRIPTDIR%/*}
    fi

    source_dir=`realpath $(dirname $0)`
    if [[ ${source_dir} == ${SCRIPTDIR} ]]; then
        cat << EOF
--------------------------------------------------------------------------------
    ERROR: You cannot specify source directory as your working directory. 
    Your option:  ${SCRIPTDIR}

EOF
        exit
    else
        ls ${SCRIPTDIR} > /dev/null 2>&1
        if [[ $? == 0 ]]; then
            cat << EOF
--------------------------------------------------------------------------------
    Warning: Directory exists - current folder will be removed.

EOF
            rm -r ${SCRIPTDIR}
        fi
    fi
}

function set_exe {
    cat << EOF
================================================================================
    Please specify the directory of GROMACS exectuables, 
    or the command to load GROMACS modules 

    Default Option:
    module load  gromacs/2021.3-mpi

EOF
    
    read -p " " EXEDIR
    EXEDIR=`echo ${EXEDIR}`

    if [[ -z ${EXEDIR} ]]; then
        EXEDIR='module load  gromacs/2021.3-mpi'
    fi

    if [[ ! -d ${EXEDIR} && (${EXEDIR} != *'module load'*) ]]; then
        cat << EOF
--------------------------------------------------------------------------------
    Error: Directory or command does not exist. Exiting current job. 

EOF
        exit
    fi

    if [[ ${EXEDIR} == *'module load'* ]]; then
        ${EXEDIR} > /dev/null 2>&1
        if [[ $? != 0 ]]; then
                    cat << EOF
--------------------------------------------------------------------------------
    Error: Module specified not available. Check your input: ${EXEDIR}

EOF
            exit
        fi
    fi
}

function set_mpi {
    cat << EOF
================================================================================
    Please specify the directory of MPI executables or mpi modules

    Default Option
    module load mpi/intel-2019

EOF
    
    read -p " " MPIDIR
    MPIDIR=`echo ${MPIDIR}`

    if [[ -z ${MPIDIR} ]]; then
        MPIDIR='module load mpi/intel-2019'
    fi

    if [[ ! -d ${EXEDIR} && (${EXEDIR} != *'module load'*) ]]; then
        cat << EOF
--------------------------------------------------------------------------------
    Error: Directory or command does not exist. Check your input: ${MPIDIR}

EOF
        exit
    fi

    if [[ ${MPIDIR} == *'module load'* ]]; then
        ${MPIDIR} > /dev/null 2>&1
        if [[ $? != 0 ]]; then
            cat << EOF
--------------------------------------------------------------------------------
    Error: Module specified not available. Check your input: ${MPIDIR}

EOF
            exit
        fi
    fi
}

function copy_scripts {
    mkdir -p ${SCRIPTDIR}
    cp ${CTRLDIR}/settings_template ${SCRIPTDIR}/settings

    cat << EOF
================================================================================
    modified scripts at ${SCRIPTDIR}/
EOF
}

# Configure settings file

function set_settings {

    # Values for keywords
    SETFILE=${SCRIPTDIR}/settings
    sed -i "/SUBMISSION_EXT/a .qsub" ${SETFILE}
    sed -i "/NCPU_PER_NODE/a 24" ${SETFILE}
    sed -i "/MEM_PER_NODE/a 50" ${SETFILE}
    sed -i "/N_THREAD/a 1" ${SETFILE}
    sed -i "/NGPU_PER_NODE/a 0" ${SETFILE}
    sed -i "/GPU_TYPE/a RTX6000" ${SETFILE}
    sed -i "/TIME_OUT/a 3" ${SETFILE}
    sed -i "/JOB_TMPDIR/a nodir" ${SETFILE}
    sed -i "/EXEDIR/a ${EXEDIR}" ${SETFILE}
    sed -i "/MPIDIR/a ${MPIDIR}" ${SETFILE}

    # Executable tabel

    LINE_EXE=`grep -nw 'EXE_TABLE' ${SETFILE}`
    LINE_EXE=`echo "scale=0;${LINE_EXE%:*}+3" | bc`
    sed -i "${LINE_EXE}a\mdrun      mpiexec              gmx_mpi mdrun -s [job].tpr     Parallel GROMACS MD calculation" ${SETFILE}

    # # Input file table - calculation performed in current directory

    # LINE_PRE=`grep -nw 'PRE_CALC' ${SETFILE}`
    # LINE_PRE=`echo "scale=0;${LINE_PRE%:*}+3" | bc`
    # sed -i "${LINE_PRE}a\[jobname].tpr        [jobname].tpr        GROMACS portable binary run input file" ${SETFILE}

    # # Reference file table - calculation performed in current directory

    # LINE_REF=`grep -nw 'REF_FILE' ${SETFILE}`
    # LINE_REF=`echo "scale=0;${LINE_REF%:*}+3" | bc`
    # sed -i "${LINE_REF}a\[refname].anything   anything             Anything" ${SETFILE}

    # # Post-processing file table

    # LINE_POST=`grep -nw 'POST_CALC' ${SETFILE}`
    # LINE_POST=`echo "scale=0;${LINE_POST%:*}+3" | bc`
    # sed -i "${LINE_POST}a\[jobname].CheckPoint\/*.cpt                Checkpoint files" ${SETFILE}
    # sed -i "${LINE_POST}a\[jobname].gro        *.gro                Gromacs geometry file" ${SETFILE}
    # sed -i "${LINE_POST}a\[jobname].trr        *.trr                Trajectory file" ${SETFILE}
    # sed -i "${LINE_POST}a\[jobname].trj        *.trj                Trajectory file" ${SETFILE}
    # sed -i "${LINE_POST}a\[jobname].log        *.log                MD output file" ${SETFILE}

    # Job submission file template - should be placed at the end of file

    cat << EOF >> ${SETFILE}
-----------------------------------------------------------------------------------
#!/bin/bash  --login
#PBS -N \${V_JOBNAME}
#PBS -l select=\${V_ND}:ncpus=\${V_NCPU}:mem=\${V_MEM}:mpiprocs=\${V_PROC}:ompthreads=\${V_TRED}\${V_NGPU}\${V_TGPU}:avx2=true
#PBS -l walltime=\${V_TWT}

echo "PBS Job Report"
echo "--------------------------------------------"
echo "  Start Date : \$(date)"
echo "  PBS Job ID : \${PBS_JOBID}"
echo "  Status"
qstat -f \${PBS_JOBID}
echo "--------------------------------------------"
echo ""

# number of cores per node used
export NCORES=\${V_NCPU}
# number of processes
export NPROCESSES=\${V_TPROC}

# Make sure any symbolic links are resolved to absolute path
export PBS_O_WORKDIR=\$(readlink -f \${PBS_O_WORKDIR})

# Set the number of threads
export OMP_NUM_THREADS=\${V_TRED}

# to sync nodes
cd \${PBS_O_WORKDIR}

# start calculation: command added below by gen_sub
-----------------------------------------------------------------------------------

EOF
    cat << EOF
================================================================================
    Paramters specified in ${SETFILE}.

EOF
}

function set_commands {
    bgline=`grep -nw "# >>> begin GROMACS job submitter settings >>>" ${HOME}/.bashrc`
    edline=`grep -nw "# <<< finish GROMACS job submitter settings <<<" ${HOME}/.bashrc`

    if [[ ! -z ${bgline} && ! -z ${edline} ]]; then
        bgline=${bgline%%:*}
        edline=${edline%%:*}
        sed -i "${bgline},${edline}d" ${HOME}/.bashrc
    fi

    echo "# >>> begin GROMACS job submitter settings >>>" >> ${HOME}/.bashrc
    echo "alias Pmdrun='${CTRLDIR}/gen_sub -x mdrun -set ${SCRIPTDIR}/settings'" >> ${HOME}/.bashrc
    echo "alias Xgmx='${CTRLDIR}/gen_sub -set ${SCRIPTDIR}/settings'" >> ${HOME}/.bashrc
    echo "alias SETgmx='cat ${SCRIPTDIR}/settings'" >> ${HOME}/.bashrc
    echo "alias HELPgmx='source ${CONFIGDIR}/run_help; print_ALIAS_HOWTO_; print_GENSUB_HOWTO_'" >> ${HOME}/.bashrc
    # echo "chmod 777 ${SCRIPTDIR}/gen_sub" >> ${HOME}/.bashrc
    # echo "chmod 777 ${SCRIPTDIR}/run_exec" >> ${HOME}/.bashrc
    # echo "chmod 777 ${SCRIPTDIR}/post_proc" >> ${HOME}/.bashrc 
    echo "# <<< finish GROMACS job submitter settings <<<" >> ${HOME}/.bashrc

    source ${CONFIGDIR}/run_help; print_ALIAS_HOWTO_
}

# Main I/O function
## Disambiguation : Here is a historical problem
## Variables and functions with 'script' in configure script refer to the user's local settings file and its directory
## In the current implementation, ${SCRIPTDIR} only has 1 file, i.e., user-defined settings file
## Executable scripts are now centralized and shared in ${CTRLDIR}
## For executable scripts, ${SCRIPTDIR} refer to their own directory. ${SETTINGS} refers to local settings file. 
CONFIGDIR=`realpath $(dirname $0)`
CTRLDIR=`realpath ${CONFIGDIR}/../`

welcome_msg
get_scriptdir
copy_scripts
set_exe
set_mpi
set_settings
set_commands
