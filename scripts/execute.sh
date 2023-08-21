#! /usr/bin/env bash

echo "Plotting datasets..."
do_save=false;
save_type="pdf";
sampling_rate=0.01;
fanout=20000;
while getopts t:r:f:s arg; do
  case $arg in
    s) do_save=true;;
    t) save_type=${OPTARG};;
    r) sampling_rate=${OPTARG};;
    f) fanout=${OPTARG};;
  esac
done

BUILD=build/build
if [ ! -f ${BUILD} ]; then
  echo "plot binary does not exist"
  exit
fi

if [ ${do_save} = true ]; then
  mkdir -p ./plot
fi

mkdir -p ./parameters

PYTHON=python
PYPLOT=plot.py
if [ ! -f ${PYPLOT} ]; then
  echo "plot python code does not exist"
  exit
fi

for dataset in $(cat scripts/datasets_for_plotting.txt); do
  echo "Building LCDE on ${dataset}"
  ${BUILD} "${dataset}" -r ${sampling_rate} -f ${fanout}
  echo "Plotting dataset ${dataset}"
  if  [ ${do_save} = true ]; then
    ${PYTHON} ${PYPLOT} -d "${dataset}" -s True -e ${save_type}
    echo "Saved plot to ./plot/${dataset}.${save_type}"
  else
    ${PYTHON} ${PYPLOT} -d "${dataset}"
  fi
done
