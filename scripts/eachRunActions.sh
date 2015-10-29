#!/bin/sh

folder=$1
settingsfile=$2
addtocumsum=$3
cutsfile=$4
if [ $# > 3 ]
 then
    cutsfile="-c $cutsfile"
fi
cd $folder
cd ..
basedir=`pwd`
cd $folder

echo "------------------------------------------------------------\n" >> log.txt
echo "------------------------------------------------------------\n" >> log.txt
echo "                 new running of eachRunActions.sh\n" >> log.txt
echo "------------------------------------------------------------\n" >> log.txt
echo "------------------------------------------------------------\n" >> log.txt

#For safety reasons, removing write privileges from accounts other than gretina
ssh gretina@a2 "cd $basedir; cd $folder; chmod 744 Global.dat GlobalRaw.dat"


#Make each of the root tree and root histogram files.
#Skip any that are already made.
if [ -f raw.root ]
then
    echo "raw.root already exists, skipping"
else
    GrROOT -i Global.dat -o raw.root -s $settingsfile -rt | tee -a log.txt
fi

if [ -f rawhists.root ]
then
    echo "rawhists.root already exists, skipping"
else
    Histos -i raw.root -o rawhists.root -s $settingsfile | tee -a log.txt
fi

if [ -f cal.root ]
then
    echo "cal.root already exists, skipping"
else
    Calculate -i raw.root -o cal.root -s $settingsfile | tee -a log.txt
fi

if [ -f calhists.root ]
then
    echo "calhists.root already exists, skipping"
else
    Cal_histos -i cal.root -o calhists.root -s $settingsfile $cutsfile -t 1 | tee -a log.txt
fi
#rm raw.root cal.root


#Only modify the cumsum if told to.
if [ "$addtocumsum" == "1" ]
    then
    #Make the runs_added file modifiable.
    fullname=`pwd`/calhists.root
    if [ -f ../runs_added.txt ]
	then
	chmod 666 ../runs_added.txt
    fi
    if [ -f ../current_cumsum.root ]
	then
	#current_cumsum exists, add to it
	if grep -Fxq "$fullname" ../runs_added.txt
	    then
	    echo -e "\e[00;31mError, added that file before.\e[00m"
	else
	    hadd -f calhists_cumsum.root ../current_cumsum.root calhists.root
	    cp ../runs_added.txt .
	    echo "$fullname" >> runs_added.txt
	    cd ..
	    ln -sf $folder/calhists_cumsum.root current_cumsum.root
	    ln -sf $folder/runs_added.txt runs_added.txt
	fi
    else
	#current_cumsum doesn't exist, make it
	cd ..
	ln -sf $folder/calhists.root current_cumsum.root
	echo "$fullname" > $folder/runs_added.txt
	ln -sf $folder/runs_added.txt runs_added.txt
    fi
    cd $basedir
    chmod 444 runs_added.txt
fi