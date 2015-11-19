#!/bin/bash
# N. Frohner, P.Wissgott 10/01/10

if [ $# -ne 1 ]; then
    echo "Usage: $0 case"
    echo "E.g.: $0 bands53-64"
    exit 0
elif [ "$1" == "-h" ]; then
    echo "converts all psink files in that directory to xsf format"
    echo "Usage: $0 case"
    echo "E.g.: $0 bands53-64"
    exit 0
fi

SEED=$1
PSINKS=`ls ${SEED}_?*.psink`
WPLOT2XSF=wplot2xsf.py

for PSINK in ${PSINKS}; do
	echo "Init wplot2xsf for $PSINK"
	BASENAME=`basename $PSINK .psink`
	WFIDX=`echo $PSINK | awk -F '[_.]' '{print $2}'`
	$WPLOT2XSF ${SEED} $WFIDX > ${BASENAME}.xsf
	gzip ${BASENAME}.xsf
	echo "Finished wplot2xsf for $PSINK"
done
