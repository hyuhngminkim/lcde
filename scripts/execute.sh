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

PLOT=build/plot
if [ ! -f ${PLOT} ]; then
  echo "plot binary does not exist"
  exit
fi

if [ ${do_save} = true ]; then
  mkdir -p ./plot
fi

for dataset in $(cat scripts/datasets_for_plotting.txt); do
  echo "Plotting dataset ${dataset}"
  if [ ${do_save} = true ]; then
    ${PLOT} "${dataset}" -s "${save}" -e "${save_type}" -r ${sampling_rate} -f ${fanout}
    echo "Saved plot to ./plot/${dataset}.${save_type}"
  else 
    ${PLOT} ${dataset} -r ${sampling_rate} -f ${fanout}
  fi
done
