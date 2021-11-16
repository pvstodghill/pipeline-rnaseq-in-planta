#! /bin/bash

export LC_ALL=C

# ------------------------------------------------------------------------

THREADS=$(nproc --all)

. config.bash

SAMPLES_INDICES=""
for i in $(seq 0 9) ; do
    if [ "${SAMPLES_R1[$i]}" ] ; then
	SAMPLES_INDICES+=" $i"
    fi
done

if [ "${SAMPLES_R2[0]}" ] ; then
    PE=1
fi

# ------------------------------------------------------------------------

if [ "$PACKAGES_FROM" = conda ] ; then
    if [ -z "$CONDA_EXE" ] ; then
	CONDA_EXE=$(type -p conda)
    fi
fi

# In order to help test portability, I eliminate all of my
# personalizations from the PATH, etc.
if [ "$PVSE" ] ; then
    export PATH=/usr/local/bin:/usr/bin:/bin
    export PERL5LIB=
    export PERL_LOCAL_LIB_ROOT=
    export PERL_MB_OPT=
    export PERL_MM_OPT=
    export PYTHONPATH=
fi

case X"$PACKAGES_FROM"X in
    XcondaX)
	CONDA_PREFIX=$(dirname $(dirname $CONDA_EXE))
	. "${CONDA_PREFIX}/etc/profile.d/conda.sh"
	conda activate $CONDA_ENV

	;;
    XhowtoX|XstubsX)
	export PATH=$(dirname ${BASH_SOURCE[0]})/stubs:"$PATH"
	;;
    XnativeX)
	: nothing
	;;
    XX)
	echo 1>&2 "\$PACKAGES_FROM is not set"
	exit 1
	;;
    X*X)
	echo 1>&2 "\$PACKAGES_FROM is recognized: $PACKAGES_FROM"
	exit 1
	;;
    *)
	echo 1>&2 "Cannot happen"
	exit 1
esac

# ------------------------------------------------------------------------

if [ -z "$PARALLEL_CMD" ] ; then
    PARALLEL_CMD="$(type -p parallel)"
fi

function run_commands {
    if [ "$PARALLEL_CMD" ] ; then
	eval $PARALLEL_CMD -j ${THREADS} -kv
    else
	bash -x
    fi
}

# ------------------------------------------------------------------------

function init_FEATURECOUNTS_ARGS {
    FEATURECOUNTS_ARGS=

    FEATURECOUNTS_ARGS+=" -O" # Assign reads to all their overlapping features

    #FEATURECOUNTS_ARGS+=" -M" # all multi-mapping reads reported alignments will be counted
    ##FEATURECOUNTS_ARGS+=" --fraction" # Assign fractional counts to features
    #FEATURECOUNTS_ARGS+=" --primary" # Count primary alignments only !(0x100)

    case X"$ORIENTATION"X in
	XforwardX)
	    FEATURECOUNTS_ARGS+=" -s 1" # stranded
	    ;;
	XreverseX)
	    FEATURECOUNTS_ARGS+=" -s 2" # reverse-stranded
	    ;;
	X*X) echo 1>&2 cannot happen ; exit 1
    esac

    if [ "$PE" ] ; then

	V="$(featureCounts -v 2>&1 | egrep . | sed -e 's/.* v//')"
	case x"$V"x in
	    x2.0.2x|x2.0.3x)
		FEATURECOUNTS_ARGS+=" -p" # input data contains paired-end reads. (>=v2.0.2)
		FEATURECOUNTS_ARGS+=" --countReadPairs" # Count read pairs (fragments) instead of reads
		;;
	    x2.0.1x)
		FEATURECOUNTS_ARGS+=" -p" # fragments (or pairs) will be counted instead of reads (<v2.0.2)
		;;
	    *)
		echo 1>&2 Unknown featureCounts version: "$V"
		exit 1
	esac

	FEATURECOUNTS_ARGS+=" -B" # Only count read pairs that have both ends aligned.
	FEATURECOUNTS_ARGS+=" -P" # Check validity of paired-end distance
	FEATURECOUNTS_ARGS+=" -C" # Only count concordant reads
    fi

}

# ------------------------------------------------------------------------

set -e
set -o pipefail

