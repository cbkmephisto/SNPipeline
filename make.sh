#!/bin/bash

# for mac osx: multiple thread make
#ncpu=`sysctl hw.ncpu | awk '{print $2}'`
#make -j $ncpu

# If eigen is not installed but still want to use most of the tools, run make without -j options!
# In this case SNPipeline.abg2bin wouldn't get compiled but everything else would.
make

old_dir=$(pwd)
if [ -d ~/bin ]
then
    cd ~/bin
    [ -d "$old_dir/bin" ] && ln -sf "$old_dir"/bin/* .
    cd -
fi

