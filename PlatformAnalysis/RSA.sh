#!/usr/bin/env bash
############################################################
## Created :  2022-02-08 
## Author  :  jsims
## Version :  0.1
## Updated :  2022-02-08
############################################################

## BASH Settings
set -o errexit
set -o pipefail
set -o noclobber
set -o nounset
#set -o xtrace

if [ "$#" -ne 6 ]; then
	if [ "$#" -gt 0 ]; then
		echo "ERROR: Incorrect number of arguments."
	fi
	
    echo "
**************************
  RegionSpecificAnalysis
**************************
USAGE: interopSummary.sh runFolder runFolder2 bed tile1 tile.run2
runFolder : Path to a run folder containing instrumentrd RTA output with h5s
runFolder2: Path to a run folder containing instrumentrd RTA output with h5s
bed : Path to a bed file containing the genommic locations of features of interest
tile: tile number (e.g 1_1011) for the tile of runFolder you want to analyse
tile.run2: tile number (e.g 1_1011) for the tile of runFolder2 you want to analyse 
output: output path and filename
***********************
"
	exit 1;
fi

runfolder=${1}
runfolder2=${2}
bed=${3}
tile=${4}
tileRun2=${5}
output=${6}

if [[ -d ${runfolder} && -d ${runfolder}/Data/InstrumentedData/ ]];
then
	echo "Runfolder : ${runfolder}"
	instrumentedRun=${runfolder}/Data/InstrumentedData/
else
    echo "${runfolder} does not exist or doesn't contain Intrumented RTA data";
    exit 1;
fi

if [[ -d ${runfolder2} && -d ${runfolder2}/Data/InstrumentedData/ ]];
then
	echo "Runfolder 2 : ${runfolder2}"
	instrumentedRun2=${runfolder2}/Data/InstrumentedData/
else
    echo "${runfolder2} does not exist or doesn't contain InterOp data";
    exit 1;
fi

if [[ -f ${bed} ]]; then
	echo "Bed file: ${bed}"
else
	echo "${bed} dose not exist."
	exit 1;
fi
bamMetrics=${instrumentedRun}/$(basename ${runfolder} | awk -v FS='_' '{print $NF}' | cut -c2-)_${tile}.h5
if [[ -f ${bamMetrics} ]]; then	
	echo "Bam metrics h5: ${bamMetrics}"
else
	echo "${bamMetrics} dose not exist"
	exit 1;
fi


bamMetrics=${instrumentedRun2}/$(basename ${runfolder2} | awk -v FS='_' '{print $NF}' | cut -c2-)_${tile}.h5
if [[ -f ${bamMetrics} ]]; then
        echo "Bam metrics h5: ${bamMetrics}"
else    
        echo "${bamMetrics} dose not exist"
        exit 1;
fi
tileMetrics=${instrumentedRun}/L00$(echo ${tile} | cut -c1)/s_${tile}.h5
if [[ -f ${tileMetrics} ]]; then
	echo "Tile metrics for run 1: ${tileMetrics}"
else
	echo "${tileMetrics} dose not exist"
	exit 1;
fi

tileMetrics=${instrumentedRun2}/L00$(echo ${tileRun2} | cut -c1)/s_${tileRun2}.h5
if [[ -f ${tileMetrics} ]]; then
        echo "Tile metrics for run 2: ${tileMetrics}"
else    
        echo "${tileMetrics} dose not exist"
        exit 1;
fi

R_VERSION="XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX/R"
#R_VERSION="XXXXXXXXXXXXXXXXXXXXXXXXX/R"
PANDOC_LOCATION="/usr/lib/rstudio-server/bin/pandoc"


if [[ ! -z ${JOB_ID+x} ]]; then
	echo "Running as an SGE job"
	## Need to hard code this as we loose the ability to do a dirname $0 inside a qsub
	QMD_FILE="XXXXXXXXXXXXXXXXXXXXX/Users/jsims/gquad/PlatformAnalysis/PlatformAnalysis/RegionSpecificAnalysis.qmd"
else
	echo "Do not running locally, please use qsub"
	exit 1;
fi
echo "${QMD_FILE}"
ARGS="list(runFolder = '${runfolder}', runFolder2 = '${runfolder2}', bed = '${bed}', tile = '$(echo ${tile} | sed 's/_/-/g')', tile.run2 = '$(echo ${tileRun2} | sed 's/_/-/g')')"
echo "ARGS: ${ARGS}"
CMD="Sys.setenv(RSTUDIO_PANDOC=\"${PANDOC_LOCATION}\");quarto::quarto_render('XXXXXXXXXXXXXXXXXXXXX/Users/jsims/gquad/PlatformAnalysis/PlatformAnalysis/RegionSpecificAnalysis.qmd', output_format = 'html', execute_params = ${ARGS}, output_file = '${output}')"

echo ${CMD} | ${R_VERSION}

