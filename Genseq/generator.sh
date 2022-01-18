x=20
for i in {2..5}
do
 a=$((x*i))
 
 python Random_Sequence.py -r d${i}_reference -o d${i}_original -u d${i}_mutated -n $a -m 1800 -e True -p 0.05 -s 0.05 -i 0.05 -l 100 -a 0.05
done
