#!/bin/sh
# embedded options to bsub - start with #BSUB
### -- set the job Name AND the job array --
#BSUB -J isoSD_sim_phi2[1-2]
### -- specify queue --
#BSUB -q hpc
### -- ask for number of cores (default: 1) --
#BSUB -n 1
### -- set walltime limit: hh:mm --
#BSUB -W 02:00
### -- specify that we need 2GB of memory per core/slot --
#BSUB -R "rusage[mem=2GB]"
### -- set the email address --
# please uncomment the following line and put in your e-mail address,
# if you want to receive e-mail notifications on a non-default address
##BSUB -u "your email here"
### -- send notification at start --
#BSUB -B
### -- send notification at completion--
#BSUB -N
### -- Specify the output and error file. %J is the job-id %I is the job-array index --
### -- -o and -e mean append, -oo and -eo mean overwrite --
#BSUB -o ./isoSD_sim_phi2/Output_%J_%I.out
#BSUB -e ./isoSD_sim_phi2/Error_%J_%I.err

# here follow the commands you want to execute
# Program_name_and_options

matlab -nodisplay -r 'isoSD_sim_phi2($LSB_JOBINDEX)' -logfile testOut
