#!/usr/bin/env python

import os
import sys

starting_folder = os.getcwd()

def RunOneRun(n,settings):
    runStr = "%04d" % n
    print "----------------------------------------------------"
    print "----------------------------------------------------"
    print "            Now starting run %s" % runStr
    print "----------------------------------------------------"
    print "----------------------------------------------------"
    os.chdir("Run%s" % runStr)
    os.system("GEB_HFC Global.dat")
    os.move("HFC.dat","HFC_%s" % runStr)
    #os.system("GrROOT_changing -i HFC_%s.dat -o raw.root -s %s -rt" % (runStr,settings))
    #os.system("	Histos -i raw.root -o rawmode2_histos.root")
    #os.system("Calculate -i raw.root -o cal.root -s %s" % settings)
    #os.chdir("../kathrin")
    #os.system("ln -s ../Run%s/cal.root cal_run%s" % (runStr,runStr))
    os.chdir(starting_folder)

def main(argv):
    if len(argv)==3:
        start = int(argv[1])
        stop = int(argv[1])
        settings = argv[2]
    elif len(argv)==4:
        start = int(argv[1])
        stop = int(argv[2])
        settings = argv[3]
    else:
        print "Unknown arguments"
        return 1

    for n in range(start,stop+1):
        RunOneRun(n)

if __name__=="__main__":
    sys.exit(main(sys.argv))
