#!/bin/bash
#w2wpath=$(find / -name "wien2wannier" 2>/dev/null)
w2wpath=`pwd`
if [ ! -d $w2wpath ]; then
 echo "$w2wpath is not a directory!"
 exit 0
fi
echo "Found wien2wannier directory: ""$w2wpath"
#w2wpathmain="$w2wpath""/w2w"
#w2wpathutil="$w2wpath""/util"
#w2wpathwplot="$w2wpath""/wplot"

echo -n " trying to modify your .cshrc file [y/n]  "
read answer
if [[ $answer = "n" ]]; then
 exit 0
elif [[ $answer = "y" ]]; then

 echo "modifying .cshrc..."
 echo "#" >>~/.cshrc
 echo "#---- wien2wannier add ----" >>~/.cshrc
 echo "setenv W2WROOT $w2wpath" >>~/.cshrc
 echo "setenv W2Wpathmain \$W2WROOT/w2w" >>~/.cshrc
 echo "setenv W2Wpathutil \$W2WROOT/util" >>~/.cshrc
 echo "setenv W2Wpathplot \$W2WROOT/wplot" >>~/.cshrc
 echo "setenv W2Wpathwoptic \$W2WROOT/woptic" >>~/.cshrc
 echo "set path = (\$W2Wpathmain \$W2Wpathutil \$W2Wpathplot \$W2Wpathwoptic \$path ) ">>~/.cshrc
 echo "#---- wien2wannier add end ----" >>~/.cshrc
 echo "added wien2wannier directories to search path in ~/.cshrc"
else
 echo "y or n?"
 exit 0
fi

