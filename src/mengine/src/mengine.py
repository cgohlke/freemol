import os

cmd = os.curdir + os.sep + "mengine.exe -dxi"
# mengine will read multi-molblock sdf files,
#  but calling through popen2 with some (large?) files hangs on windows,
#  even though the same file runs fine at command line.
#  So, split into separate files.
for molblock_in in (file("drugs.sdf").read()).split('$$$$\n'):
#for molblock_in in (file("caffeine.sdf").read()).split('$$$$\n'):
#for molblock_in in (file("acaffeine.sdf").read()).split('$$$$\n'):
 if len(molblock_in) > 0:
  mengine_in, mengine_out, mengine_err = os.popen3(cmd)
  mengine_in.write(molblock_in)
  mengine_in.close()
  molblock_out = mengine_out.read()
  print molblock_out,
  print mengine_err.read(),
