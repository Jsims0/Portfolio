#!/usr/bin/env nextflow

//defaults
nextflow.enable.dsl=2
inputPath = "XXXXXXXXXXXXXXXXXXXXX/Users/jsims/"

//import sub-workflow
include { LAUNCH_RTA } from '/home/jsims/RTA-C-automation/nextflowScripts/main.nf'

process checkFlags{
//detect if run has been processed before
input:
    path dir

output:
    env passRun
    
    shell:
    """
    if [[  -f $inputPath/$dir/RTA_PP_Launched.txt ]]; then
        passRun=''
    else
        passRun='$inputPath/$dir/'
    fi
        

    """
}

process getRTA{
//get rta version from run info

input:
    path dir

output:
    env rta_version

    shell:
    """
    if [[ -f $inputPath/$dir/Logs/info_00000.log ]]; then
        rta_version=\$(cat $inputPath/$dir/Logs/info_00000.log | head -n 1 |awk -v FS=' ' '{print \$5}' | sed 's/v//g')
	else
        rta_version='WARNING: No $inputPath/$dir/Logs/info.0000 to retrive rta verion from. Skipping'
    fi
    """
}

process getInterop{
//get interop version from run info

input:
	path dir

output:
    env interop_version

    shell:
    """
    if [[ -f $inputPath/$dir/Logs/info_00000.log ]]; then
        echo "interop: Logs found"
        interop_version=\$(cat $inputPath/$dir/Logs/info_00000.log | grep InterOp: | head -n 1 |awk -v FS=' ' '{print \$5}')
        if ! [[ \$interop_version =~ ^[0-9]+.[0-9]+.[0-9]+\$ ]]; then
            echo "Interop: Incompatable versions extracted from logs"
        fi
    else
            interop_version="WARNING: No Logs/info_00000.log to retrive rta verion from. Attempting with default"
        fi
    """

}


process selectTiles{
//get tiles with images in run folder

input:
	path dir

output:
	stdout
shell:
"""
$params.appDir/selectTile.sh $inputPath $dir
"""
}

process gatherLanes{
//get lane numbers from runfolder

input:
	path dir

output:
	stdout

shell:
"""
$params.appDir/gatherLanes.sh $inputPath $dir
"""

}


process buildTileSet{

input:
	val lanes
	val tiles

output:
	stdout

script:
lanesList = lanes.split()
tilesList = tiles.split()
combinedList = [lanesList, tilesList].combinations()

result = combinedList.collect { sublist -> sublist.join('_') }.take(8)

"""
echo $result
"""
}

process buildCommand{

input:
	val file
	val tileSet
	val rta
	val interop

output:
	tuple val(file)

script:
	cmd = ["runFolder":"${file}", "tileSet":tileSet, "rta":rta, "interop":interop]

}



// setup a chennel from runfolders inside inputdir
file_ch = Channel
        .fromPath( "${inputPath}/**" )
	.filter{ it.name =~ /(RunInfo.xml)$/ }
	.map { it.parent }


//file_ch.choice( positive_ch, negative_ch ) { it -> fromPath(it).collect()  }

//runInData_ch = tmpFolders.join(tmpTiles)
workflow {

//gather run information
path=checkFlags(file_ch).view()
rta=getRTA(file_ch).view()
interop=getInterop(file_ch).view()
tiles=selectTiles(file_ch).view()
lanes=gatherLanes(file_ch).view()

//subset lane/tile combinations above a pf limit
tile_set=buildTileSet(lanes, tiles)
//build rta-c command
//cmd=buildCommand(file_ch, tile_set, rta, interop)

//submit rta-c to other script
LAUNCH_RTA(file_ch, tile_set, rta, interop)
}
