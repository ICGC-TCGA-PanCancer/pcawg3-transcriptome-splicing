import os
import sys
import pandas as pd
import glob

__conf_ext = '.conf'
conf_dir = os.path.join(os.path.dirname(__file__), 'confs')
conf_files = glob.glob(os.path.join(conf_dir, '*' + __conf_ext))
scheme_names = [os.path.basename(conf).lower().replace(__conf_ext, '') for conf in conf_files]

def _load_scheme(name):
    name = name.lower()
    if not name in scheme_names:
        msg = 'Scheme %s does not exist.\nScheme names: %s' %(name, ', '.join(scheme_names))
        raise ValueError
    conf = os.path.join(conf_dir, name + __conf_ext)
    scheme = dict(pd.read_csv(conf, sep='\t').values)
    for key in scheme.keys():
        try:
            int_key = int(key)
        except ValueError:
            continue
        else:
            scheme[int_key] = scheme[key]
    return scheme

def get_color_scheme(scheme_name, labels=None, return_scheme=False):
    '''Given a data type, loads a corresponding color scheme specified by ICGC.

    Returns a color pallete matching terms to colors speficied by the scheme.
    Color palletes are implemented as a dictionary matching terms to hex colors
    Color scheme configuration files are located in ./confs

    scheme_name (str): specifies which pallete to use, (from ./confs/scheme_name.conf)
    labels (array-like): terms corresponding to colors. If specified, a list of hex
        colors is returned.
    return_scheme (bool): if labels is specified, returns the pallete dictionary
        in addition to the color list
    '''
    scheme = _load_scheme(scheme_name)
    if labels is None:
        return scheme
    colors = list()
    for lab in labels:
        try:
            colors.append(scheme[lab])
        except KeyError:
            rep = ', '.join([str(key) for key in scheme])
            raise ValueError('Label %s not in scheme\nLabels: %s'%(lab, rep))
    if return_scheme:
        return colors, scheme
    else:
        return colors
