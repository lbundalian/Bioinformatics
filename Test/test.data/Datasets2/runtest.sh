#!/bin/bash
for i in {1..12}
do
 cp test.sh D$i/
 cp FireMSA.py D$i/
 cp maxalign.pl D$i/
 cd D$i
 printf "d$i\n" | bash test.sh 
 cd ..
done
