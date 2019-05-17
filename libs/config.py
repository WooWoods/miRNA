"""
    config module
"""

import os
import types
import errno
import yaml
from collections import UserDict


class Config(UserDict):
    """Instantiate a configuration object."""

    def __init__(self, cfgfile):
        self.data = {}

        if os.path.isfile(cfgfile):
            self.cfgfile = cfgfile
            self.from_yaml()
        

    def from_yaml(self):
        d = yaml.load(open(self.cfgfile), Loader=yaml.Loader)
        self.data.update(d)

    def __setitem__(self, key, value):
        self.data[key] = value

    def __getitem__(self, key):
        return self.data[key]
