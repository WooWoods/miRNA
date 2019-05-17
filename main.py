import os
import sys

from libs.config import Config
from libs.fastqc import PreProcessing
from libs.mapping import NCFilter, RefMapping
from libs.quantifier import Quantifier



if __name__ == '__main__':
    cfg = sys.argv[1]
    config = Config(cfg)

    pp = PreProcessing(config)
    pp.process()

    ncf = NCFilter(config)
    ncf.process()
    
    ref = RefMapping(config)
    ref.process()
    
    quanti = Quantifier(config)
    quanti.process()

