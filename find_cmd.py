#!/usr/bin/env python3

'''
This function implements the (misnamed but very useful) ksh utility
type which, given an application name and the user's PATH, returns
the application's pathname or ""
The list of directories given by PATH has been extended to
include PYTHON_PATH (via sys.path).

For human users with a real PATH this works quite well; for
web apps it is not so useful so the basic search has been supplemented
with some likely places
'''

LIKELY_LOCS = ["/usr/bin", "/usr/local/bin", "/sw/bin", "/opt/bin"]

import os, sys

def find_cmd(name) :
  if "PATH" in os.environ :
    pathlist = os.environ["PATH"].split(":") + ["."]
  else :
    pathlist = ["."]
  pathlist = pathlist + sys.path + LIKELY_LOCS
  # This covers the case where an absolute path has been given!
  if os.path.isfile(name) :
    return(name)
  for p in pathlist :
    if os.path.isfile(p + "/" + name) :
      return(p + "/" + name)
  return(None)


if __name__ == "__main__" :
  if len(sys.argv) != 2 :
    sys.stderr.write("Usage: %s <commmand>\n" % sys.argv[0])
    sys.exit(0)

  full_path = find_cmd(sys.argv[1])
  if full_path != None:
    print(full_path)
  else:
    print(sys.argv[1] + " not found")
