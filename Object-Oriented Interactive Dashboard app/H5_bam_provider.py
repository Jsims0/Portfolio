import re
import h5py
import numpy as np
import pandas as pd
import codecs
from urllib.request import urlopen
import io
import sys, os
from random import sample
from os import path, listdir
from os.path import basename
from XXXXXXXXXXXXXXXXXXXXXXX import instrumentation_output


class H5LocalProvider(object):
    def __init__(self):
        pass

    @staticmethod
    def find_h5_file(runpath:str):
        """find H5 file"""
        H5_files = []
        H5_files.extend(
            [f for f in listdir(runpath+'/Data/InstrumentedData') if f.endswith('h5')])
        return H5_files

    @staticmethod
    def select_bam_metric_h5(files, run):
        """select the bam metric file from H5 files"""
        run_name = basename(run)
        bam_h5 = [file for file in files if file.split('_')[0] in run_name]

        return bam_h5[0]

    @staticmethod
    def list_bam_h5s(files, run):
        """list the bam metric file from H5 files"""
        run_name = basename(run)
        bam_h5s = []
        bam_h5s = [file for file in files if file.split('_')[0] in run_name]

        return bam_h5s

    @staticmethod
    def select_cluster_h5(files, run):
        """select the cluster metric file from H5 files"""
        if path.exists(run+'/Data/InstrumentedData/ClusterMetrics.h5'):
            cluster.h5 = run+'/Data/InstrumentedData/ClusterMetrics.h5'
            return cluster.h5
        else:
            return None


    @staticmethod
    def bam_data(runpath, file):
        h5 = h5py.File(path.join(runpath, 'Data/InstrumentedData/'+file))
        return np.array(h5.get('ID')).tolist()


    @staticmethod
    def select_tileH5(runpath, lane, tile):
        tile_h5_path = path.join(runpath, '/Data/InstrumentedData/', lane, tile, 'h5')
        if path.exists(tile_h5_path):
            return tile_h5_path
        else:
            return None

    @staticmethod
    def list_lanes(runpath):
        """list lanes in run directory"""
        inst_path = instrumentation_output(runpath)
        lanes = list(filter(lambda x: x.startswith("L"), listdir(inst_path)))
        return lanes

    @staticmethod
    def list_tiles(runpath, lanes):
        """list the tile h5 files"""
        if isinstance(lanes, str) == True:
            lanes = [lanes]

        inst_path = instrumentation_output(runpath)
        h5_files = []
        for lane in lanes:
            h5_files.extend(
                [f for f in listdir(path.join(inst_path, lane)) if f.endswith('h5')]
            )
        return h5_files

    @staticmethod
    def get_tile(runpath, lane, file):
        """get tile name from file name and select from list of h5 files"""
        tile_name = file.split('_')[2]
        inst_path = instrumentation_output(runpath)


        if isinstance(lane, str) == True:
            lanes = [lane]

        inst_path = instrumentation_output(runpath)
        h5_files = []
        for lane in lanes:
            h5_files.extend(
                [f for f in listdir(path.join(inst_path, lane)) if f.endswith('h5')]
            )


        tile = [i for i in h5_files if tile_name in i]
        return tile[0]


    @staticmethod
    def get_lane(run, file):
        """list lanes in run directory"""

        lane_num = file.split('_')[1]
        inst_path = instrumentation_output(run)
        lanes = list(filter(lambda x: x.startswith("L"), listdir(inst_path)))

        lane = [i for i in lanes if i.endswith(lane_num)]
        return str(lane[0])



###################################################################

    def load_files(self, lane, tile_file, bam_file, run):
        """ceate h5 file handles"""

        inst_path = instrumentation_output(run)
        self.bam_file = h5py.File(path.join(inst_path, bam_file), 'r')
        self.tile_file = h5py.File(path.join(inst_path, lane, tile_file), 'r')
        self.cluster_file = h5py.File(path.join(inst_path, 'ClusterMetrics.h5'), 'r')

    def get_h5_dimension(self):
        """get the number of cycels in bam h5 file"""
        return self.bam_file.get('cycleMisMatch').shape[0]

    def box_plot_data(self, run, n_cycle, n_samples):
        """generate data frame for dot plot"""
        warn = [None, None]
        try:
            bam_file = self.bam_file
            tile_file = self.tile_file
        except Exception as e:
            return ({'data':[]}, ['HTF5 file not found', 'danger'])

        pf = np.where(np.array(tile_file.get('PassingFilter')) == 1)[1]
        mismatches = np.isin(np.array(bam_file.get('cycleMisMatch')[n_cycle]), np.array([65, 67, 71, 84]))

        idx_mismatches = np.where(mismatches == True)[0]
        if n_samples > idx_mismatches.shape[0]:
            warn = ['Sampling Number greater then number of mismatches', 'warning']
            n_samples = idx_mismatches.shape[0]
        else:
            pass
        idx_mismatches = idx_mismatches[np.random.choice(np.arange(len(idx_mismatches)), n_samples, replace=False)]

        pf_idx_mismatches = pf[idx_mismatches]

        refcall = np.array(bam_file.get('cycleMisMatch')[n_cycle][idx_mismatches])
        basecall = np.array(tile_file.get('Basecall')[n_cycle][pf_idx_mismatches])
        basecall_sample = basecall[:1000]

        intensity = np.array(tile_file.get('FullyCorrectedIntensity')[n_cycle][pf_idx_mismatches])

        chastity = np.array(tile_file.get('Chastity')[n_cycle][pf_idx_mismatches])
        mapQ = np.array(bam_file.get('mapq')[0][idx_mismatches])

        df = pd.DataFrame({'ref': [codecs.decode(i) for i in refcall],
                           'basecall': [codecs.decode(i) for i in basecall],
                           'ch1': [i[0] for i in intensity],
                           'ch2': [i[1] for i in intensity],
                           'Chastity': chastity,
                           'mapQ': mapQ})

        return (df, warn[0], warn[1])
