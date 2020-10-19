#!/bin/sh
# embedded options to bsub - start with #BSUB
### -- set the job Name AND the job array --
#BSUB -J My_array[1-25]
### –- specify queue --
#BSUB -q hpc
### -- ask for number of cores (default: 1) --
#BSUB -n 4
### -- set walltime limit: hh:mm --
#BSUB -W 03:00
### -- specify that we need 2GB of memory per core/slot --
#BSUB -R "rusage[mem=2GB]"
### -- set the email address --
# please uncomment the following line and put in your e-mail address,
# if you want to receive e-mail notifications on a non-default address
##BSUB -u your_email_address
### -- send notification at start --
#BSUB -B
### -- send notification at completion--
#BSUB -N
### -- Specify the output and error file. %J is the job-id %I is the job-array index --
### -- -o and -e mean append, -oo and -eo mean overwrite --
#BSUB -o Output_%J_%I.out
#BSUB -e Error_%J_%I.err

# here follow the commands you want to execute
# Program_name_and_options

matlab -nodisplay -r 'test($LSB_JOBINDEX)' -logfile testOut
