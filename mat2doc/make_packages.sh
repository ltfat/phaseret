#!/bin/bash
cd ..
mat2doc . mat --tgz --unix --script=release.py

mat2doc . mat --zip --dos --addon=phaseret_win64 --packagename=phaseret-%s-win64 --script=release.py


