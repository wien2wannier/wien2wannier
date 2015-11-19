#!/bin/bash
#w2wpath=$(find / -name "wien2wannier" 2>/dev/null)
w2wpath=`pwd`
if [ ! -d $w2wpath ]; then
 echo "$w2wpath is not a directory!"
 exit 0
fi
echo "Found wien2wannier directory: ""$w2wpath"
w2wpathmain="$w2wpath""/w2w"
w2wpathutil="$w2wpath""/util"
w2wpathwplot="$w2wpath""/wplot"
w2wpathwoptic="$w2wpath""/woptic"
echo -n " trying to modify your .bashrc file [y/n]  "
read answer
if [[ $answer = "n" ]]; then
 exit 0
elif [[ $answer = "y" ]]; then
 echo "modifying .bashrc..."
 echo "#---- wien2wannier add ----" >>~/.bashrc
 echo "export W2WROOT=$w2wpath" >>~/.bashrc
 echo "export PATH=""$""PATH:$w2wpathmain:$w2wpathutil:$w2wpathwplot:$w2wpathwoptic:" >>~/.bashrc
 echo "#---- wien2wannier add end ----" >>~/.bashrc
 echo "added wien2wannier directories to search path in ~/.bashrc"
else
 echo "y or n?"
 exit 0
fi

