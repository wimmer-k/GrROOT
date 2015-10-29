#!/usr/bin/env python

import os
from sys import argv

command = 'tar zcvf - %s | ssh %s "cat > %s"'

from_loc = argv[1]
dest_machine = argv[2]
to_loc = argv[3]

command = command % (from_loc, dest_machine, to_loc)
print command
os.system(command)
