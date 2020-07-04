#!/bin/bash
dir=doc
if [ -d $dir ]
then
    pdoc --html src/buckinghampi.py --html-dir ./doc --overwrite
else
    mkdir $dir
    pdoc --html src/buckinghampi.py --html-dir ./doc
fi