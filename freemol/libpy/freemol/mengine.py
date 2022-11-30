# mengine module for the FreeMOL library

# executable path
import traceback
import os

_mengine_exe = None

def validate():
    '''
    make sure that the FREEMOL environment variable is defined and
    that it contains a valid directory path
    
    returns 1 if valid, 0 if not
    '''
    ok = 0
    global _mengine_exe
    if _mengine_exe != None: # _mengine_exe assigned, so presume valid
        ok = 1
    else:
        from . import find_exe
        test_path = find_exe('mengine')
        if test_path:
            _mengine_exe = test_path
            ok = 1
    return ok

def run(input):
    from subprocess import Popen, PIPE
    result = None
    global _mengine_exe
    if _mengine_exe == None:
        validate()
    if _mengine_exe != None:

        if not isinstance(input, bytes):
            input = input.encode()

        try:
            p = Popen([_mengine_exe], stdin=PIPE, stdout=PIPE, stderr=PIPE)
            result = p.communicate(input)
        except:
            print("Error: mengine did not run")
            traceback.print_exc()
    return result
    
