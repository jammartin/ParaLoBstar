#!/bin/bash

TEMPLATE="standalonePlot.tex"
OUT_DIR="../plots"

for pgfFile in ../pgf/*.tex
do
    sed -i "~" "s,PGF_FILE,${pgfFile%.tex},g" $TEMPLATE
    pdflatex -output-directory $OUT_DIR $TEMPLATE
    base=$(basename $pgfFile)
    mv $OUT_DIR/${TEMPLATE%.tex}.pdf $OUT_DIR/${base%.tex}.pdf
    rm -f $OUT_DIR/${TEMPLATE%.tex}*
    mv $TEMPLATE~ $TEMPLATE
done
