#!/bin/bash

if [ -z "$ROOTSYS" ]; then
    echo "Please set your ROOTSYS variable to the location of your root installation."
    echo "Preferably, place this in your .bashrc"
    exit 1
fi

if [ ! -d bin ]; then
    echo "Making bin directory"
    mkdir -p bin/bin32
    mkdir -p bin/bin64
fi

if [ ! -d lib ]; then
    echo "Making lib directory"
    mkdir -p lib/lib32
    mkdir -p lib/lib64
fi

#There are some files which will be modified frequently on local machines,
# but should not be updated in the repository.
#Copy the repository version into to working folder version.
cp src/RawEventLoopBase.cc src/RawEventLoop.cc
cp src/CalEventLoopBase.cc src/CalEventLoop.cc

GRROOT_DIR=`pwd`

echo "Making env.sh"
echo "#!/bin/bash" > env.sh
echo "" >> env.sh 

#Standard setup variables.
echo "export GRROOT_DIR=$GRROOT_DIR" >> env.sh
echo "" >> env.sh

#System-dependent setup variables.
echo "if [ \`file /bin/bash | awk '{print \$3}'\` = \"64-bit\" ]; then" >>  env.sh
echo "    export GRROOT_BINDIR=\$GRROOT_DIR/bin/bin64" >> env.sh
echo "    export GRROOT_LIBDIR=\$GRROOT_DIR/lib/lib64" >> env.sh
echo "else" >> env.sh
echo "    export GRROOT_BINDIR=\$GRROOT_DIR/bin/bin32" >> env.sh
echo "    export GRROOT_LIBDIR=\$GRROOT_DIR/lib/lib32" >> env.sh
echo "fi" >> env.sh
echo "export PATH=\$PATH:\$GRROOT_BINDIR" >> env.sh
echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$GRROOT_LIBDIR" >> env.sh

chmod +x env.sh

echo "In addition, please make the following changes manually."
echo "  In ~/.rootrc, add \$GRROOT_LIBDIR to the Unix.*.Root.DynamicPath"
echo "  In your rootlogon.C, add the following lines."
echo "     gSystem->Load(\"libGretina\");"
echo "     gSystem->Load(\"libS800\");"
echo "     gSystem->Load(\"libScaler\");"
echo "     gSystem->Load(\"libSettings\");"
echo "     gSystem->Load(\"libRunInfo\");"
echo "The program will run without these changes, but you will be unable"
echo "to interact with the objects in the root interpreter beyond histogramming"
