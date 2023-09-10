#!/usr/bin/env nextflow

nextflow.enable.dsl=2
IPA_SCRIPT_HOME=params.IPA_SCRIPT_HOME
process RTA {

input:
	val file
        val tileSet
        val rta
        val interop

output:
   	val result

script:

	file = file.normalize().toString()
	if (params.output != ''){
                output = params.output
        }else{
        	runName = file.split('/').last()
        	output = file + '/rtaCluster/' + runName
        }
	result=output
	
"""
echo $result
"""
}

process findClusterNeighbours{
    label 'bigTask'

input:
    path outRunFolder

output:
    path outRunFolder

    """
    echo 'FCN: $outRunFolder'
    $IPA_SCRIPT_HOME/scripts/postprocessing/findClusterNeighbours.R --outputFile \
    $outRunFolder/Data/InstrumentedData/ClusterNeighbours.h5 $outRunFolder
    """
}

process InteropSummary {
label 'bigTask'

input:
    path outRunFolder

output:
    path outRunFolder

    """
    echo IS: $outRunFolder
    """

}

process fcToError {
label 'bigTask'

input:
    path outRunFolder

    output:
    path outRunFolder

    """
    echo 'fcToError: $outRunFolder'
    \$IPA_SCRIPT_HOME/scripts/postprocessing/fcToErrorH5.sh  -r $outRunFolder  -o $outRunFolder/Data/InstrumentedData --bwa
    """

}

process OIA {
label 'bigTask'

input:
    val outRunFolder
    val inRunFolder
output:

    script:
    print("OIA: $inRunFolder $outRunFolder")
    
    """
    { echo '/launchOIA.sh $inRunFolder $outRunFolder'; } 2> /dev/null   
    """
}

process cache {
label 'bigTask'

input:
    path outRunFolder
    val fromfcToError
output:
    path outRunFolder

    """
    echo $outRunFolder
    echo $fromfcToError
    \$IPA_SCRIPT_HOME/scripts/postprocessing/cacheH5Results.R $outRunFolder
    """

}

process mergeOIA {
label 'bigTask'

input:
    val outRunFolder
    val OIA
output:
    path outRunFolder

    """
    echo mergeOIA: $outRunFolder
    \$IPA_SCRIPT_HOME/../mergeOIAandInteropSummary.R $outRunFolder
    """

}

process clean {
label 'bigTask'

input:
    path outrunFolder

output:


    '''
    '''

}

process primeta {
label 'bigTask'

input:
    path outRunFolder

output:
    path outRunFolder

    """
    echo primeta: $outRunFolder
    \${IPA_SCRIPT_HOME}/../git/quality_models/runXXXXXXX/.sh -- --data_path '\$(ls $outRunFolder/*WithIQS.csv )' --save_path '${outRunFolder}s/Primeta/'
    """
}

process SAReport {
label 'bigTask'

input:
    path outRunFolder

output:
    path outRunFolder

    """
    $IPA_SCRIPT_HOME/../reporting/SAReport.sh $outRunFolder
    """
}

workflow LAUNCH_RTA {
take:
	file
	tileSet
	rta
	interop
main:
//launch RTA
	rtaOutput = RTA(file, tileSet, rta, interop).view()

//post RTA
// cluster + fcToError | cacheResults
clusterData = findClusterNeighbours(rtaOutput)
fcError = fcToError(clusterData)
results = cache(clusterData, fcError)

//clusters = cluster()
//tiles = fcToError()
//cacheResults(clustered, tiled).view()

// interopsummery + OIA | mergeOIA | primeta
interops=InteropSummary(rtaOutput)
images=OIA(file, interops)
//mergedInterops=mergeOIA(rtaOutput, images)
//PRI=primeta(mergedInterops)
}

