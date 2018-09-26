#!/bin/bash
cd ..
make -C mex -f Makefile_crossmingw
/home/susnak/dev/mat2doc/mat2doc.py . mat --tgz --unix --script=release.py --packagename=phaseret-%s-src
/home/susnak/dev/mat2doc/mat2doc.py . mat --zip --dos --addon=phaseret_win64 --packagename=phaseret-%s-win64 --script=release.py
/home/susnak/dev/mat2doc/mat2doc.py . html


