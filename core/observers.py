# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 10:45:28 2016

@author: Johannes Hartung
"""

import numpy as np

class AbstractObserver(object):
    __init__(self):
        """ 
        Get actualized from to be observed object
        """
        super(AbstractObserver, self).__init__()
    
    setValues(self, vals):
        raise NotImplemented()
    

