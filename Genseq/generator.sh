x=500
for i in {1..12}
do
 a=$((x*i))
 
 python Random_Sequence.py -r d${i}_reference -o d${i}_original -u d${i}_mutated -n $a -m 3000 -e True -p 0.05 -s 0.20 -i 0.00 -l 100 -a 0.0
done
