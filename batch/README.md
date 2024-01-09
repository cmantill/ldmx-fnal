# SLURM Batch System

To submit jobs one can use SLAC's SDF SLURM system.
More general instructions can be found on SDF's [README](https://github.com/slaclab/sdf-docs/blob/master/batch-compute.md#slurmexample)

An example is shown in template.slurm where you would replace ```<config file>``` with the config file of your choice. Note for the SLURM batch system, you must pass the entire environment through the job submission script since it does not pass either the environment you are currently working in or your bashrc. Since you are working inside the container, you can simply source the ldmx-env.sh script in you job submission script.

```bash
source "$LDMX_BASE/ldmx-sw/scripts/ldmx-env.sh"
```

In addition, make sure to source the ldmx-env.sh script before job submission in order to properly set the ```LDMX_BASE``` path. In order to submit a job, simply use the ```sbatch``` command.

```bash
sbatch <slurm file>
```

## Other Useful Commands

To check the status of your jobs:
```bash
squeue -u <computing id>
```

To cancel a job:
```bash
scancel <job id>
```