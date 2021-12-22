#!/bin/bash
for i in {1..12}
do
 cp test.sh D$i/
 cd D$i
 printf "d$i\n" | bash test.sh 
 cd ..
done
