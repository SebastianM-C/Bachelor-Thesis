#!/usr/bin/env python

import os


def exists(b, d, n):
    """Check if the folder specified by the given parameters exists"""
    return os.path.isdir("../Output/B" +
                         str(b) + " D" + str(d) + " N" + str(n))


def cd(b, d, n):
    """ Try to cd to the given path. In case of an error go back to ../../Scripts
    and try again (maybe the last run had an error or
    the script did not reach the end)"""

    try:
        os.chdir("../Output/B" + str(b) + " D" + str(d) + " N" + str(n))
    except FileNotFoundError as fnf:
        try:
            print(fnf.strerror + ':' + fnf.filename + '\nTrying again')
            os.chdir("../../Scripts")
            os.chdir("../Output/B" + str(b) +
                     " D" + str(d) + " N" + str(n))
            print('Succes\n' + 'Now in: ' + os.getcwd())
        except FileNotFoundError:
            print('Check the filename\n' + 'Now in: ' + os.getcwd())


def format_float(n):
    """Converts the number to int if possible"""
    return ('{:n}' if n == int(n) else '{:.8g}').format(n)