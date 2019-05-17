"""
Quantifying miRNA expression.
"""

import os
import re
import subprocess
from collections import defaultdict
from multiprocessing import Pool

from .utils import dir_check


class Quantifier:
    def __init__(self, config):
        self.config = config
        self.output = self.config.get('OUTPUT')
        self.report = os.path.join(self.config.get('REPORT'), 'read_count')
        dir_check(self.report)
        self.noprecursor = ['hsa-miR-12119','hsa-miR-12120','hsa-miR-12121',
                'hsa-miR-12122','hsa-miR-12123','hsa-miR-12124','hsa-miR-12125',
                'hsa-miR-12126','hsa-miR-12127','hsa-miR-12128','hsa-miR-12129',
                'hsa-miR-12130','hsa-miR-12131','hsa-miR-12132','hsa-miR-12133',
                'hsa-miR-12135','hsa-miR-12136']

    def record_mature_pre(self):
        self.mat2pre = defaultdict(set)
        mature_to_pre = self.config['REF']['M2P']

        with  open(mature_to_pre) as fh:
            for line in fh:
                larr = line.split('\t')
                m_id = larr[0].split()[0]
                pre_id = larr[2]
                self.mat2pre[m_id].add(pre_id)

    def reads_map_time_count(self, sample):
        reads_map_time = defaultdict(int)
        map_res = os.path.join(self.output, sample, f'{sample}.mature.bwt')

        with open(map_res) as fh:
            for line in fh:
                larr = line.split('\t')
                readid = larr[0].split()[0]
                reads_map_time[readid] += 1
        return reads_map_time

    def reads_map_hairpin(self, sample):
        reads2hairpin = defaultdict(set)
        map_res = os.path.join(self.output, sample, f'{sample}.hairpin.bwt')
        
        with open(map_res) as fh:
            for line in fh:
                larr = line.split('\t')
                readid = larr[0].split()[0]
                reads2hairpin[readid].add(larr[2])
        return reads2hairpin

    def read_count(self, sample):
        map_res = os.path.join(self.output, sample, f'{sample}.mature.bwt')
        reads_map_time = self.reads_map_time_count(sample)
        reads2hairpin = self.reads_map_hairpin(sample)
        counter = defaultdict(int)

        output = os.path.join(self.report, f'{sample}.count')

        with open(map_res) as fh:
            for line in fh:
                larr = line.split('\t')
                readid = larr[0].split()[0]
                if reads_map_time.get(readid) > 1:
                    continue
                mat_id = larr[2]
                mapped_hairpin = reads2hairpin.get(readid)
                hairpin_of_mat = self.mat2pre.get(mat_id)

                # reads mapped to mature miRNA and related hairpin simultaneously will be
                # count as a correct mapping, multi alignment reads will be ignored
                if mat_id in self.noprecursor:
                    counter[mat_id] += 1
                    continue

                if mapped_hairpin & hairpin_of_mat and (mapped_hairpin | hairpin_of_mat == hairpin_of_mat):
                    counter[mat_id] += 1

        with open(output, 'w') as fh:
            for mir in counter:
                fh.write(f'{mir}\t{counter.get(mir)}\n')
        #return counter

    def process(self):
        self.record_mature_pre()
        with Pool(8) as p:
            p.map(self.read_count, self.config.get('SAMPLES'))











