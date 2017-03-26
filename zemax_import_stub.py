# -*- coding: utf-8 -*-
"""
"""

def zmx2dict(filename):
    """
    Reads a Zemax file and converts it to a python dictionary.
    
    @param: filename (str)
            filename of the .zmx file. May be utf8 or ansi.
            
    @return: zmxdict (dict)
             Dictionary keys are the first 4 characters of each line 
             in the zemax file.
             Dictionary values are lists of lines beginning with that key.
             Indented lines are treated as if they were 
             the same line as the previous one.
    """
    f = open(filename,"r")
    zmx = f.read() 
    if zmx.startswith("\xff\xfe"): # utf8
        zmx = zmx[2::2].replace('\"', "") # dirty way to convert to ascii
    f.close
    
    modeblocks = zmx.split("\r\nMODE ")
    # Before the first occurence of "MODE", 
    # there is the Zemax version and dongle number of the file creator.
    # Then either "MODE SEQ" or "MODE NSC" occurs, 
    # followed by a sequential or non-sequential lens decription. 
    
    for current_modeblock in modeblocks[1:]:
        if current_modeblock.startswith("SEQ"):
            # Current block describes a lens in Zemax sequential raytracing mode
            
            # Split at \r\n windows newline commands,
            # except if a line begins with an indent 
            lines = current_modeblock.split("\r\n")[1:]
            data = []
            for l in lines:
                if l.startswith("  "):
                    data[-1] += l + "\r\n"
                else:
                    data.append(l)
                    
            # Group lines that begin with the same 4 characters
            zmxdict = {}
            for d in data:
                if d[0:4] in zmxdict:
                    zmxdict[d[0:4]] += [d[5:]]
                else:
                    zmxdict[d[0:4]]  = [d[5:]]
        else:
            raise Exception("Nonsequential mode import not supported yet.")
    return zmxdict
    
zmxdict = zmx2dict("lenssystem.ZMX")