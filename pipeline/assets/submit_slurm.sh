#!/bin/bash

# Dynamic messages about standard error log files
function err_msg() {
  err_nfiles=`find $1 -name "e_*.out" -type f ! -empty | wc -l`
  if [[ $err_nfiles > 0 ]]; then
    err_nlines=`find $1 -name "*e_*.out" -exec cat {} + | wc -l`
    echo -e "$RED_CROSS ${err_nlines} lines written to ${err_nfiles} standard error log files"
  else
    echo -e "$GREEN_TICK 0 lines written to standard error log files"
  fi
}

# Queue status for user without the column headers
function squ_pr() {
  squeue -u $1 -h -t pending,running -r --start
}


# Current number of pending/running jobs in queue for user
function njobs() {
  squ_pr $1 | wc -l
}

# Elapsed time since first job submitted for user.
# Total time when last job has completed.
# Reference: https://www.linuxquestions.org/questions/linux-newbie-8/time-difference-calculation-4175459414/
function elapsed() {
  if [[ `njobs $1` > 0 ]]; then
    first=`squ_pr $1 | awk '{print $6}' | head -n 1 | sed 's/^.*T//'`
  else
    first=$2
  fi
  now=$(date +"%T")
  SEC1=`date +%s -d ${first}`
  SEC2=`date +%s -d ${now}`
  DIFFSEC=`expr ${SEC2} - ${SEC1}`
  date +%H:%M:%S -ud @${DIFFSEC}
}

# Print number of total jobs to submit
jtotal=${#file_to_submit[@]}
echo -e "$BLUE_BULLET Total jobs: $jtotal"

# Store job ids in a job array
jobarray=()

# Store start time
start=$(date +"%T")

# Step 1: Loop over array of shell scripts and submit batch jobs
for shname in "${file_to_submit[@]}"; do
	# Check if data set was specified
	if [ "$shname" = "" ]; then
		  echo "Shell script must be specified."
		  exit 1
	fi

	# Check current max job limit
  mkdir -p $logDir
  startQLength=`njobs $user`
  curr_max_submit=$(($maxQueueLength - $startQLength))

  # Job submission
  if command -v sbatch &>/dev/null; then
      # Ensure no jobs are submitted if current max job limit is reached
      while [[ "$currQLength" -ge "$curr_max_submit" ]]; do
        currQLength=`njobs $user`
        echo -e "$BLUE_BULLET Waiting on: ${currQLength} jobs to complete before submitting job"
        sleep 30s
      done

      # Submit shell script to queue without extra verbosity using variable qq
		  qcmd=$(sbatch -D $logDir $shname)
      echo -e "$GREEN_TICK Submitted to queue: $shname"
  fi
done
echo -e "$GREEN_TICK Job submission finished. Check status with \"squeue -u ${user}\""

# Step 2: wait until all jobs are complete before proceeding

# Initialize queue length at starting number of jobs
if command -v sbatch &>/dev/null; then
  currQLength=`njobs $user`
else
  currQLength=0
fi

sleep 5s
while [[ $currQLength > 0 && shouldWait ]]; do
    currQLength=`njobs $user`

    if [[ $currQLength > 0 ]]; then
      echo -e "$BLUE_BULLET Time elapsed: `elapsed ${user} ${start}` | Jobs remaining: ${currQLength}/${jtotal}"
      sleep 30s
    else
      echo -e "$GREEN_TICK All jobs completed in: `elapsed ${user} ${start}`"
    fi
done

# Step 3: Verbose printouts
err_msg $logDir
echo -e "$GREEN_BULLET Logs written to \"${logDir}\""
if [[ -z "$outputDir" ]]; then
  echo -e "$GREEN_BULLET No outputs written"
else
  echo -e "$GREEN_BULLET Outputs written to \"${outputDir}\""
fi
