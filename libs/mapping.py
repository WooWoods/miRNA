"""
    Small non-coding RNA filtering, know miRNA idenfication.
"""

import os
import subprocess
from multiprocessing import Pool

from .utils import dir_check


class NCFilter:
    """mapping against Rfam database to filter non-coding RNAs"""

    def __init__(self, config):
        self.config = config
        self.output = self.config.get('OUTPUT')
        self.samples = self.config.get('SAMPLES')
        self.rfam = self.config['REF']['RFAM']
        self.thread = self.config['THREADS']['MAPPING']

    def db_align(self, sample):
        out_dir = os.path.join(self.output, sample)

        try:
            subprocess.run([
                'bowtie',
                '--threads', '3',
                '--seedmms', '1',
                '--seedlen', '18',
                '-a',
                '--best',
                '--strata',
                '--sam',
                self.rfam,
                '--al', os.path.join(out_dir, f'{sample}.rfam.mapped.fq'),
                '--un', os.path.join(out_dir, f'{sample}.rfam.unmapped.fq'),
                os.path.join(out_dir, f'{sample}.trim.fq.gz'),
                os.path.join(out_dir, f'{sample}.rfam.sam')
                ])
        except subprocess.CalledProcessError:
            pass

    def stats(self):
        pass

    def process(self):
        with Pool(self.thread) as p:
            p.map(self.db_align, self.samples)
    

class RefMapping:
    """idenfy known miRNAs by mapping against miRBase"""

    def __init__(self, config):
        self.config = config
        self.output = self.config.get('OUTPUT')
        self.samples = self.config.get('SAMPLES')
        self.mature = self.config['REF']['MATURE']
        self.hairpin = self.config['REF']['HAIRPIN']
        self.thread = self.config['THREADS']['MAPPING']

    def hairpin_mapping(self, sample):
        out_dir = os.path.join(self.output, sample)

        try:
            subprocess.run([
                'bowtie',
                '--threads', '3',
                '--seedmms', '1',
                '--seedlen', '18',
                '-a',
                '-v', '0',
                '-m', '5',
                '--best',
                '--strata',
                '--norc',
                self.hairpin,
                '--al', os.path.join(out_dir, f'{sample}.hairpin.mapped.fq'),
                '--un', os.path.join(out_dir, f'{sample}.hairpin.unmapped.fq'),
                os.path.join(out_dir, f'{sample}.rfam.unmapped.fq'),
                os.path.join(out_dir, f'{sample}.hairpin.bwt')
                ])
        except subprocess.CalledProcessError:
            pass

    def mature_mapping(self, sample):
        out_dir = os.path.join(self.output, sample)

        try:
            subprocess.run([
                'bowtie',
                '--threads', '3',
                '--seedmms', '1',
                '--seedlen', '18',
                '-a',
                '-v', '0',
                '-m', '5',
                '--best',
                '--strata',
                '--norc',
                self.mature,
                os.path.join(out_dir, f'{sample}.rfam.unmapped.fq'),
                os.path.join(out_dir, f'{sample}.mature.bwt')
                ])
        except subprocess.CalledProcessError:
            pass
    
    def mapping(self, sample):
        self.hairpin_mapping(sample)
        self.mature_mapping(sample)

    def process(self):
        with Pool(self.thread) as p:
            p.map(self.mapping, self.samples)

    
