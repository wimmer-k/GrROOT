#!/usr/bin/env python

import os
import sys

def FileInCWD(filename):
    return filename in os.listdir(os.getcwd())

def PreserveCWD(f):
    def Inner(*args,**kwargs):
        starting_folder = os.getcwd()
        f(*args,**kwargs)
        os.chdir(starting_folder)
    return Inner

@PreserveCWD
def RunOneRun(n):
    runStr = "%04d" % n
    print "----------------------------------------------------"
    print "----------------------------------------------------"
    print "            Now starting run %s" % runStr
    print "----------------------------------------------------"
    print "----------------------------------------------------"
    #Bail out if the folder is not there.
    if not FileInCWD("Run%s"%runStr):
        return
    os.chdir("Run%s" % runStr)
    
    #Bail out if the file is not there.
    if not FileInCWD("Global.dat") or os.stat("Global.dat").st_size==0:
        print "Nothing in Global.dat, or it is missing"
        return

    #Unpacking with GrROOT
    rawfile = "raw"
    if FileInCWD("%s.root"%rawfile):
        print "unpacked file already exists, skipping"
    else:
        command = "GrROOT_devbuild -i Global.dat -o %s.root -s ~/kathrin/settings/inbeam2.dat -rt" % rawfile
        print command
        os.system(command)

    #Find any files into which root has overflowed.
    rawfilelist = [rawfile]
    for i in range(1,10):
        name = "%s_%d"%(rawfile,i)
        if FileInCWD(name+".root"):
            rawfilelist.append(name)

    #Histogramming with Histos
    #if FileInCWD("rawhists.root"):
    #    print "histograms already exist, skipping"
    #else:
    #    command = "Histos_devbuild -i %s -o rawhists.root -s ~/kathrin/settings/inbeam2.dat" % " ".join(s+".root" for s in rawfilelist)
    #    print command
    #    os.system(command)

    #Calibrating
    calfile = "cal"
    calfilelist = []
    for name in rawfilelist:
        calname = name.replace(rawfile,calfile)
        calfilelist.append(calname)
        if FileInCWD(calname+".root"):
            print "Calibrated file %s.root exists, skipping" % calname
        else:
            command = "Calculate_devbuild -i %s.root -o %s.root -s ~/kathrin/settings/inbeam2.dat" % (name,calname)
            print command
            os.system(command)

    #Make the calibrated histograms
    #if FileInCWD("calhists.root"):
    #    print "calibrated histograms already exist, skipping"
    #else:
    #    command = "Cal_histos_devbuild -i %s -o calhists.root -s ~/kathrin/settings/inbeam2.dat" % " ".join(s+".root" for s in calfilelist)
    #    print command
    #    os.system(command)

    #for name in calfilelist:
    #    command = "ln -s %.root ../cal_links"
    #    print command
    #    os.system(command)

def main(argv):
    if len(argv)==2:
        start = int(argv[1])
        stop = int(argv[1])
    elif len(argv)==3:
        start = int(argv[1])
        stop = int(argv[2])
    else:
        print "Unknown arguments"
        return 1

    for n in range(start,stop+1):
        RunOneRun(n)

if __name__=="__main__":
    sys.exit(main(sys.argv))
