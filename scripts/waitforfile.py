#!/usr/bin/env python

import os
from sys import argv
from time import sleep

n = int(argv[1])

while not os.path.exists("Run%04d/Global.dat"%n):
    sleep(0.25)
    
print "File now exists"
        
while os.stat("Run%04d/Global.dat"%n).st_size==0:
    sleep(0.25)
            
print "File filled"
