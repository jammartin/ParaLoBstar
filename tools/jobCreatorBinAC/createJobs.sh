#!/bin/bash

# BEGIN CONFIGURATION
SIM_DIR="../../simulations"
modes=("Le" "Hi")
lbIntervals=("1" "10" "50")
numProcs=("76" "94")
thetas=("0.5" "0.8")
# END CONFIGURATION

if [ -z "$1" ]
then
    echo "Please provide initial distribution file stem as command line argument. Aborting."
    exit
fi

INIT_FILE_STEM=$1
WORKDIR=$(pwd)

for theta in "${thetas[@]}"
do

    # create neccessary directories
    mkdir -p $SIM_DIR/${INIT_FILE_STEM}th${theta//./_}
    cd $SIM_DIR/${INIT_FILE_STEM}th${theta//./_}
    mkdir -p config jobs

    echo "Copying initial distribution file from $HOME if available."
    cp ~/$INIT_FILE_STEM*.h5 $INIT_FILE_STEM.h5

    echo "Creating configs and scripts in '$SIM_DIR/${INIT_FILE_STEM}th${theta//./_}' ..."

    # creating submit script
    cp $WORKDIR/templates/submit.sh .
    echo "echo '... submitting jobs ...'" >> submit.sh

    # main loop
    for mode in "${modes[@]}"
    do
	if [ "$mode" = "Se" ]; then
	    
	    JOBNAME=$INIT_FILE_STEM$mode
	    echo "Creating serial job '$JOBNAME'"
	    
	    # create configuration file
	    cfgPath=config/$JOBNAME.info

	    cp $WORKDIR/templates/config.info $cfgPath
	    
	    sed -i "s/THETA/$theta/g" $cfgPath
	    sed -i "s/INITFILE/$INIT_FILE_STEM.h5/g" $cfgPath
	    sed -i "s/OUTDIR/$JOBNAME/g" $cfgPath
	    sed -i "s/PARALLEL/false/g" $cfgPath
	    sed -i "s/HILBERT/false/g" $cfgPath
	    sed -i "s/LBINTERVAL/1/g" $cfgPath

	    #create job script
	    jobPath=jobs/$JOBNAME.sh

	    cp $WORKDIR/templates/job.sh $jobPath

	    sed -i "s/_JOBNAME_/$JOBNAME/g" $jobPath
	    sed -i "s/_NODES_/1/g" $jobPath
	    sed -i "s/_PPN_/1/g" $jobPath

	    # adding job to submit
	    echo "qsub jobs/$JOBNAME.sh" >> submit.sh
	else
	    for lbInterval in "${lbIntervals[@]}"
	    do
		for numProc in "${numProcs[@]}"
		do
		    JOBNAME="${INIT_FILE_STEM}lb${lbInterval}np$numProc$mode"
		    echo "Creating parallel job '$JOBNAME'"

		    # create configuration file
		    cfgPath=config/$JOBNAME.info

		    cp $WORKDIR/templates/config.info $cfgPath
		    
		    sed -i "s/THETA/$theta/g" $cfgPath
		    sed -i "s/INITFILE/$INIT_FILE_STEM.h5/g" $cfgPath
		    sed -i "s/OUTDIR/$JOBNAME/g" $cfgPath
		    sed -i "s/PARALLEL/true/g" $cfgPath

		    if [ "$mode" = "Hi" ]; then
			sed -i "s/HILBERT/true/g" $cfgPath
		    else
			sed -i "s/HILBERT/false/g" $cfgPath
		    fi
		    
		    sed -i "s/LBINTERVAL/$lbInterval/g" $cfgPath

		    #create job script
		    jobPath=jobs/$JOBNAME.sh

		    cp $WORKDIR/templates/job.sh $jobPath

		    coresPerNode=28
		    let nodes=numProc/coresPerNode

		    if [ "$nodes" -eq "0" ]; then
			nodes=1 # numProc < coresPerNode
			coresPerNode=$numProc
		    else
			let nodesRemainder=numProc%coresPerNode
			if [ "$nodesRemainder" -gt "0" ]; then
			    let nodes=nodes+1
			    let coresPerNode=numProc/nodes
			fi
			let missingNodes=numProc-nodes*coresPerNode
			if [ "$missingNodes" -gt "0" ]; then
			    coresPerNode+="+1:ppn=$missingNodes"
			fi
		    fi
		    
		    sed -i "s/_JOBNAME_/$JOBNAME/g" $jobPath
		    sed -i "s/_NODES_/$nodes/g" $jobPath
		    sed -i "s/_PPN_/$coresPerNode/g" $jobPath

		    # adding job to submit
		    echo "qsub jobs/$JOBNAME.sh" >> submit.sh
		done
	    done
	fi
    done

    echo "echo '... done.'" >> submit.sh
    chmod u+x submit.sh
    echo "... done."
done   
