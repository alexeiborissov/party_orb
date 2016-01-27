#!/bin/bash
for file in multipar_fmt_t*.dat; do sed -e '1d' $file > `basename $file .dat`.nohead; done