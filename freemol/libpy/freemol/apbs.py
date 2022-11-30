# apbs module for the FreeMOL library

# executable path
import os

_apbs_exe = None
_psize_py = None

def validate():
    '''
    make sure that the FREEMOL environment variable is defined and
    that it contains a valid directory path
    
    returns 1 if valid, 0 if not
    '''
    import itertools
    import sys

    global _apbs_exe
    if _apbs_exe == None:
        path = []
        if 'FREEMOL' in os.environ:
            path.append(os.path.join(os.environ['FREEMOL'], 'bin'))
            lib_path = os.path.join(os.environ['FREEMOL'], 'lib')
            LLP = "DYLD_LIBRARY_PATH" if sys.platform == "darwin" else "LD_LIBRARY_PATH"
            os.environ[LLP] = lib_path + ("" if LLP not in os.environ else
                    os.pathsep + os.environ[LLP])

        if 'PATH' in os.environ:
            path.extend(os.environ['PATH'].split(os.pathsep))
        for (p, suffix) in itertools.product(path, ('.exe', '')):
            test_path = os.path.join(p, 'apbs' + suffix)
            if os.path.isfile(test_path):
                _apbs_exe = test_path
                break
        else:
            return 0

    global _psize_py
    if _psize_py == None:
        path = []
        if 'FREEMOL' in os.environ:
            path.append(os.path.join(os.environ['FREEMOL'], 'share', 'apbs', 'tools', 'manip'))
            path.append(os.path.join(os.environ['FREEMOL'], 'share', 'pdb2pqr', 'src'))
        for p in path:
            test_path = os.path.join(p, 'psize.py')
            if os.path.isfile(test_path):
                _psize_py = test_path
                break
        else:
            return 0

    return 1

def get_exe_path():
    if _apbs_exe == None:
        validate()
    return _apbs_exe

def get_psize_path():
    if _psize_py == None:
        validate()
    return _psize_py

