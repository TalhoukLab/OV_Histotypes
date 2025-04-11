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

# Specify R path to place in $PATH. If empty, logout.
if [ "$RPath" = "" ]; then
  echo "Path to R directory must be specified"
  exit 1
fi
