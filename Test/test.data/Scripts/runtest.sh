#!/bin/bash
for i in {1..12}
do
 cd D$i
 printf "d$i\n" | bash test.sh 
 cd ..
done
