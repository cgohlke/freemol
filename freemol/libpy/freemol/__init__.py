# master FreeMOL module file

import sys

def quote_exe(exe):
   if ' ' in exe:
      if sys.platform == 'win32':
         exe = '"'+exe+'"'
   return exe


def find_exe(name):
    import os
    import itertools

    path = []
    if 'FREEMOL' in os.environ:
        path.append(os.path.join(os.environ['FREEMOL'], 'bin'))
    if 'PATH' in os.environ:
        path.extend(os.environ['PATH'].split(os.pathsep))
    d = os.path.dirname(sys.executable)
    if d not in path:
        path.append(d)

    if sys.platform.startswith('win'):
        suffixes = ('.exe',)
    else:
        suffixes = ('.exe', '')

    for (p, suffix) in itertools.product(path, suffixes):
        exe = os.path.join(p, name + suffix)
        if os.path.isfile(exe):
            return exe

    return None
