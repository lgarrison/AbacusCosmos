'''
Created on Feb 2, 2012

@author: dferrer
Updates by Lehman Garrison

Wrapper for input files. Files consist of two parts: a parameters section that essentially looks like python code, and a binary section. When this class is instantiated
for a file, it attempts to 'register' all of the parameters by simply executing this portion of of the file as slightly modified python code. This is a sketchy way to
do it, but it avoids the difficulty of having to rewrite ParseHeader in python. These parameters can then be accessed as class variables. Right now, full functionality
is only available for scalar data. Any parameter type that isn't understood is turned into a string, which can be parsed elsewhere.

e.g if infile.txt has the lines

a = 1
b = 2 2 2
c = 35v kls a

then if infile = InputFile('infile.txt')

infile.a returns 1
infile.b returns [2,2,2]
infile.c returns "35v kls a"

EXERCISE EXTREME CAUTION: THIS CLASS ALLOWS ARBITRARY EXECUTION OF TEXT FILES AS CODE
'''

import sys
import os.path as path
from cStringIO import StringIO

class InputFile:
    def __init__(self, fn=None, str_source=None):
        """
        Construct an InputFile from filename `fn`, or a string containing the file contents `str_source`.
        """
        
        # Read the input as either a file or string
        if fn and str_source:
            raise RuntimeError('Cannot specify both `fn` = "{}" and `str_source` = "{}"!'.format(fn, str_source))
        if not fn and not str_source:
            raise RuntimeError('Must specify one of `fn` and `str_source`!')
            
        if fn:
            param = open(fn, "rb")
        elif str_source:
            param = StringIO(str_source)
        else:
            raise RuntimeError('Invalid state. Should never be able to get here!')
            
        self.code = param.readlines()
        param.close()
        
        for line in self.code:
            line = line.strip()
            comment = line.find('#')
            if comment > -1:  # Trim everything after a comment
                line = line[:comment]
            equals = line.find('=')
            if equals == -1:
                continue
            try:
                exec 'self.'+line  # valid python as-is?
            except (SyntaxError,NameError):
                try:
                    lhs = line[equals+1:]
                    items = lhs.split()
                    vec  = ','.join(items)
                    exec 'self.'+line[0:equals+1] + "("+vec+")"  # valid as a vector, if we replace spaces with commas?
                except:
                    try:
                        exec 'self.'+line[0:equals+1] + "'"+line[equals+1:] +"'"  # valid as a string if we wrap the whole RHS in quotes?
                    except:
                        raise RuntimeError('Error: Could not parse line: "{:s}"'.format(line))
        
        # A common pathology is that BoxSize is interpreted as an int by ParseHeader
        # We should endeavor to always write "50." instead of "50" in the .par files
        for field in ['BoxSize', 'InitialRedshift', 'ZD_PLT_target_z', 'wa']:
            if field in self:
                setattr(self, field, float(self[field]))

    # Allows access as "params[key]"
    def __getitem__(self, key):
        if key not in vars(self):
            raise KeyError("InputFile has no field {}".format(repr(key)))
        return getattr(self, key)
    def __setitem__(self, key, value):
        return setattr(self, key, value)
    
    # Allows testing "key in params"
    def __contains__(self, key):
        return hasattr(self, key)
    
    # Don't return 'code' when asking for a string representation
    def __repr__(self):
        selfvars = vars(self)
        try:
            del selfvars['code']
        except KeyError:
            pass
        return str(selfvars)
    