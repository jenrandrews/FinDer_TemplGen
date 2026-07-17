#!/bin/bash

idir='event_inis'
if [ ! -d $idir ]; then
  echo 'No ini files available in '$idir
  exit
fi

odir='event_outputs'
if [ ! -d $odir ]; then
  mkdir $odir
fi

avodir='averaged_gmms'
if [ ! -d $avodir ]; then
  mkdir $avodir
fi

script_dir=$(dirname $BASH_SOURCE)

for fpath in ${idir}/*.ini; do
  ini=$(basename $fpath)
  source_id="${ini%.*}"
  av_outputs="${avodir}/${source_id}.csv"

  if ! test -f ${av_outputs}.gz; then
    echo "${source_id} NOT PROCESSED"
  fi

  if test -f ${av_outputs}.gz; then
    echo "Skipping ${source_id}: already processed"
    continue
  fi

  if ! test -f ${av_outputs}; then
    echo $fpath

    ## Run and create outputs
    oq engine --run $fpath
    mkdir -p "${odir}/${source_id}"
    oq engine --export-outputs -1 "${odir}/${source_id}"
    calc_id=$(python "${script_dir}"/get_job_id.py "${source_id}")
    oq engine --lo $calc_id

    ## Using the AvgGMPE method
    cp ${odir}/${source_id}/gmf-data*.csv ${av_outputs}
   
    ## Clean up
    oq engine --dc ${calc_id} -y
    rm -rf "${odir}/${source_id}"

  else
    echo "Skipping ${source_id}: already processed"
  fi
done
