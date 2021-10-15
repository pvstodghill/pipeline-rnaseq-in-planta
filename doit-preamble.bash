#! /bin/bash

set -e
set -o pipefail

# ------------------------------------------------------------------------

THREADS=$(nproc --all)

export LC_ALL=C

# ------------------------------------------------------------------------

. config.bash

# ------------------------------------------------------------------------

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

# ------------------------------------------------------------------------

if [ -z "${DONT_USE_STUBS}" ] ; then
    export PATH=$(dirname ${BASH_SOURCE[0]})/stubs:"$PATH"
fi
    
