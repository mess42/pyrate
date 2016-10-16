# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 10:45:28 2016

@author: Johannes Hartung
"""

import numpy as np

class AbstractObserver(object):
    def __init__(self):
        """ 
        Get actualized from to be observed object.
        Most simple case: Observed object calls informUpdate()
        """
        super(AbstractObserver, self).__init__()
    
    def informAboutUpdate(self):
        raise NotImplemented()
    

