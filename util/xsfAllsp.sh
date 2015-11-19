#!/bin/bash
# N. Frohner, P.Wissgott 10/01/10
if [ "$1" == "-h" ]; then
    echo "converts all psink files in that directory to xsf format"
    echo "(spin-polarized version)"
    echo "Usage: $0 case"
    echo "E.g.: $0 bands53-64"
    exit 0
fi

if [ $# -ne 2 ]; then
    echo "Usage: $0 -up/-dn case"
    echo "E.g.: $0 -up bands53-64"
    exit
fi

SEED=$2
if [ "$1" == "-up" ]; then
   SP="up"
elif [ "$1" == "-dn" ]; then
   SP="dn"
elif [ "$1" == "-soup" ]; then
   SP="up"
   cp $SEED.woutso $SEED.woutup
   cp "$SEED""_centres.xyzso" "$SEED""_centres.xyzup"
elif [ "$1" == "-sodn" ]; then
   SP="dn"
   cp $SEED.woutso $SEED.woutdn
   cp "$SEED""_centres.xyzso" "$SEED""_centres.xyzdn"
else
   echo "Usage: $0 -up/-dn <seed>"
   echo "E.g.: $0 -up bands53-64"
   exit
fi



PSINKS=`ls ${SEED}_?*.psink$SP`
WPLOT2XSF=wplot2xsfsp.py
SPI="--"$SP
for PSINK in ${PSINKS}; do
	echo "Init wplot2xsf for $PSINK"
	BASENAME=`basename $PSINK .psink$SP`
	WFIDX=`echo $PSINK | awk -F '[_.]' '{print $2}'`
	$WPLOT2XSF $SPI ${SEED} $WFIDX > ${BASENAME}.xsf$SP
	gzip ${BASENAME}.xsf$SP
	echo "Finished wplot2xsf for $PSINK"
done
