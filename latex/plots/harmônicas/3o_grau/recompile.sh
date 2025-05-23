#!/bin/bash

set -x

#for fn in "laço_histerese" "indução" "campo" "componentes"
for fn in "laço_histerese" "indução" "campo" "componentes"
do
gnuplot $fn.gpi
done

pdflatex teste.tex

rm *aux* *log *out
