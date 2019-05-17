"""
    Implement commonly used funtions
"""
import os
import sys


def file_check(f):
    return os.path.exists(f)

def dir_check(d):
    if os.path.exists(d):
        return None
    os.makedirs(d)

def fqfile_check(fq):
    if os.path.is_file(fq) and os.path.getsize(fq):
        return True
    return False

