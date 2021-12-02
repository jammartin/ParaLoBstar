#!/bin/bash

# BEGIN CONFIGURATION 
lbIntervals=("1" "10" "50")
numProcs=("2" "4" "8" "16" "28" "40" "56" "76" "94" "112" "140")
# END CONFIGURATION

if [ -z "$1" ]
then
    echo "Please provide path to to simulation directory and job name base as command line argument. Aborting."
    exit
fi

JOB_BASE_PATH=$1
PATH_REGEX=".*th([0-1]_[0-9]).*N([0-9]+)$"

if [[ $JOB_BASE_PATH =~ $PATH_REGEX  ]]
then
    THETA_LBL=${BASH_REMATCH[1]}
    NUM_PARTICLES=${BASH_REMATCH[2]}

theta=${THETA_LBL/_/.}

echo "Inferred from given job base path: N = $NUM_PARTICLES, theta=$theta"

for lbInterval in "${lbIntervals[@]}"
do

    outFile="pgf/$(basename ${JOB_BASE_PATH})lb${lbInterval}th${THETA_LBL}.tex"

    echo "Creating $outFile ..."

    #serial
    cat <<EOF > $outFile
\begin{tikzpicture}
  \begin{semilogyaxis}[height=8cm, width=10cm, legend pos=north east, title={\$N=\num{${NUM_PARTICLES}} \text{, } \theta=\num{${theta}}$},
               grid=major, xlabel=Number of cores, ylabel=Execution time per step in $\si{\milli\second}$]
    \addplot[only marks, mark=o] coordinates {
      (1,$(./main.py -s -f ${JOB_BASE_PATH}Se.h5 -o))
    };
    \addlegendentry{serial}
    \addplot[only marks, mark=triangle] coordinates {
EOF

    # Lebesgue
    for numProc in "${numProcs[@]}"
    do
	echo "      ($numProc, $(./main.py -l $lbInterval -f ${JOB_BASE_PATH}lb${lbInterval}np${numProc}Le.h5 -o))" >> $outFile
    done

    cat <<EOF >> $outFile
    };
    \addlegendentry{Lebesgue}
    \addplot[only marks, mark=square] coordinates {
EOF

    # Hilbert
    for numProc in "${numProcs[@]}"
    do
	echo "      ($numProc, $(./main.py -l $lbInterval -f ${JOB_BASE_PATH}lb${lbInterval}np${numProc}Hi.h5 -o))" >> $outFile
    done

    cat <<EOF >> $outFile
    };
    \addlegendentry{Hilbert}
  \end{semilogyaxis}
\end{tikzpicture}
EOF

    echo "... done."
done

else
    echo "Error inferring number of particles and theta from filename. - Exiting."
    echo "Please provide path to h5 profiling files in the format '/path/*th[01]_[0-9]/*N[0-9]+'."
fi




