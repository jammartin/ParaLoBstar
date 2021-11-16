#!/bin/bash

THETA_LBL="th0_5"
lbIntervals=("1" "5" "10" "25" "50")
numProcs=("2" "4" "8" "12" "16" "20" "28" "36")

if [ -z "$1" ]
then
    echo "Please provide path to to simulation directory and job name base as command line argument. Aborting."
    exit
fi

JOB_BASE_PATH=$1

for lbInterval in "${lbIntervals[@]}"
do
    
outFile="pgf/$(basename ${JOB_BASE_PATH})lb${lbInterval}${THETA_LBL}.tex"

echo "Creating $outFile ..."

#serial
cat <<EOF > $outFile
\begin{tikzpicture}
  \begin{axis}[height=9cm, width=9cm, legend pos=north west,
               grid=major, xlabel=numProcs, ylabel=t]
    \addplot coordinates {
      (1,$(./main.py -s -f ${JOB_BASE_PATH}Se.h5 -o))
      (28,$(./main.py -s -f ${JOB_BASE_PATH}Se.h5 -o))      
    };
    \addlegendentry{serial}
    \addplot coordinates {
EOF

# Lebesgue
for numProc in "${numProcs[@]}"
do
    echo "      ($numProc, $(./main.py -l $lbInterval -f ${JOB_BASE_PATH}lb${lbInterval}np${numProc}Le.h5 -o))" >> $outFile
done

cat <<EOF >> $outFile
    };
    \addlegendentry{Lebesgue}
    \addplot coordinates {
EOF

# Hilbert
for numProc in "${numProcs[@]}"
do
    echo "      ($numProc, $(./main.py -l $lbInterval -f ${JOB_BASE_PATH}lb${lbInterval}np${numProc}Hi.h5 -o))" >> $outFile
done

cat <<EOF >> $outFile
    };
    \addlegendentry{Hilbert}
  \end{axis}
\end{tikzpicture}
EOF

echo "... done."
done



