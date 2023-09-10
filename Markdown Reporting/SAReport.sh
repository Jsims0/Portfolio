#!/usr/bin/env bash
############################################################
## Created :  2022-04-28
## Author  :  smacarthur
## Version :  1.0.0
## Updated :  2022-04-28
############################################################

## SGE Settings
#$ -pe threaded 20
#$ -l 'excl'
## #$ -l cpuflags='*:avx2:*' ## Mainly for RTA
#$ -cwd
#$ -V
#$ -M $USER@illumina.com
#$ -m es

## BASH Settings
set -o errexit
set -o pipefail
set -o noclobber
set -o nounset
#set -o xtrace

## Command Line Parameters
if [ "$#" -lt 1 ]; then
    echo "
***********************
      SAReport
smacarthur@illumina.com
     2022-04-28
***********************

USAGE: SAReport runFolder

runFolder : Path to a run folder containing an InterOp directory

***********************
"
    exit 1;
fi

runfolder=${1}
outdir=${runfolder}

if [[ -d ${runfolder} ]];
then
    echo "Runfolder : ${runfolder}"
else
    echo "${runfolder} does not exist or doesn't contain InterOp data";
    exit 1;
fi



R_VERSION="/illumina/sync/software/unofficial/corona/bin/R"
PANDOC_LOCATION="/usr/lib/rstudio-server/bin/pandoc"

RMD_FILE=$(readlink -f $(dirname ${0})//SAReport.Rmd)
echo $RMD_FILE

ARGS="--args ${runfolder}"
CMD="Sys.setenv(RSTUDIO_PANDOC=\"${PANDOC_LOCATION}\");rmarkdown::render(\"${RMD_FILE}\", intermediates_dir = \"${outdir}\", knit_root_dir = \"${outdir}\", output_file = \"${outdir}/$(basename ${runfolder})_SAReport.html\")"
R_PARAMS="--no-init-file --quiet --no-save"
echo ${CMD} | ${R_VERSION} ${R_PARAMS} ${ARGS}


echo "************************************************************"
echo "DONE"
echo "************************************************************"



