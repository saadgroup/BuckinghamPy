#!/bin/bash
dir=doc
if [ -d $dir ]
then
    pdoc --html buckinghampy/buckinghampi.py --html-dir ./doc --overwrite
else
    mkdir $dir
    pdoc --html buckinghampy/buckinghampi.py --html-dir ./doc
fi