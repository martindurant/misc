# -*- coding: utf-8 -*-
"""
Code to "Open" a file, i.e., call the os's function, which in turn
should find the correct program for the given file
"""

import sys,os,subprocess

def launch(filename):
    """Pass filename to the OS, which should select the correct
    program to open it.
    """
    if sys.platform=='win32':
        os.startfile(filename)
    elif sys.platform=='darwin':
        subprocess.Popen(['open', filename])
    else:
        try:
            subprocess.Popen(['xdg-open', filename])
        except OSError:
            print 'File open failed: '+filename