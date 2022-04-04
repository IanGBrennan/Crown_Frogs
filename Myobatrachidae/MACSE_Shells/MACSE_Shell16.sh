cd /Users/Ian/Desktop/Crown_Frogs/Myobatrachidae
java -jar /Applications/macse_v2.03.jar -prog refineAlignment -align /Users/Ian/Desktop/Crown_Frogs/Myobatrachidae/Myobatrachidae_L98.fasta -optim 1 -local_realign_init 0.1 -local_realign_dec 0.1 -fs 10
java -jar /Applications/macse_v2.03.jar -prog refineAlignment -align /Users/Ian/Desktop/Crown_Frogs/Myobatrachidae/Myobatrachidae_L99.fasta -optim 1 -local_realign_init 0.1 -local_realign_dec 0.1 -fs 10
echo finished batch 16 of 16
perl -pi -w -e 's/!/N/g;' *.fasta
exit 0
