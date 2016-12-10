#!/bin/bash

# for mac osx: multiple thread make
#ncpu=`sysctl hw.ncpu | awk '{print $2}'`
#make -j $ncpu

# If eigen is not installed but still want to use most of the tools, run make without -j options!
# In this case SNPipeline.abg2bin wouldn't get compiled but everything else would.
make

if [ -d ~/bin ]
then
    cd ~/bin
    [ -d ~/SNPipeline/bin ] && ln -sf ~/SNPipeline/bin/* .
    cd -
fi