#!/bin/bash

# Specify the script directory. If empty, logout.
if [ "$scriptDir" = "" ]; then
  echo "Script directory must be specified"
	exit 1
fi

# Specify the output directory. If empty, logout.
if [ "$outputDir" = "" ]; then
  echo "Output directory must be specified"
	exit 1
fi

# Specify the log directory. If empty, logout.
if [ "$logDir" = "" ]; then
  echo "Log directory must be specified"
	exit 1
fi

# Specify the normalization method. If empty, logout.
if [ "$norm_by" = "" ]; then
  echo "Normalization method must be specified"
  exit 1
fi

# Specify the minimum variance. If negative, logout.
if ! [[ $min_var -ge 0 ]]; then
  echo "Minimum variance must be greater than or equal to 0"
	exit 1
fi

# Specify the normalization type. If empty, logout.
if [ "$norm_type" = "" ]; then
  echo "Normalization type must be specified"
	exit 1
fi

# Specify R path to place in $PATH. If empty, logout.
if [ "$RPath" = "" ]; then
  echo "Path to R directory must be specified"
  exit 1
fi
