#! /bin/bash

###
### chromosight quantification for loops and strips
###

module load python/3.11.2

res=(1000000 500000 250000 100000 50000 40000 25000 20000 10000 5000 2500 1000 800 600 400)

for i in "${res[@]}"
do
        distal_size=$(cooler info distal.${i}_balanced.cool | grep sum | tr -dc '0-9')
	proximal_size=$(cooler info proximal.${i}_balanced.cool | grep sum | tr -dc '0-9')

	if [ ${distal_size} -gt ${proximal_size} ];
	then
		mat_size=${proximal_size}
	else
		mat_size=${distal_size}
	fi
	echo "CC size = ${distal_size}"
	echo "UC size = ${proximal_size}"
	echo "mat subsample size = ${mat_size}"


        chromosight quantify --pattern loops \
        --subsample ${mat_size} \
        --win-fmt npy \
        distal_${i}_loops.tsv \
        distal.${i}_balanced.cool \
        distal_${i}_loops_distalquant

        chromosight quantify --pattern loops \
        --subsample ${mat_size} \
        --win-fmt npy \
        distal_${i}_loops.tsv \
        proximal.${i}_balanced.cool \
        distal_${i}_loops_proximalquant

        chromosight quantify --pattern loops \
        --subsample ${mat_size} \
        --win-fmt npy \
        proximal_${i}_loops.tsv \
        distal.${i}_balanced.cool \
        proximal_${i}_loops_distalquant

        chromosight quantify --pattern loops \
        --subsample ${mat_size} \
        --win-fmt npy \
        proximal_${i}_loops.tsv \
        proximal.${i}_balanced.cool \
        proximal_${i}_loops_proximalquant


  	chromosight quantify --pattern stripes_left \
        --subsample ${mat_size} \
        --win-fmt npy \
        distal_${i}_stripes_left.tsv \
        distal.${i}_balanced.cool \
        distal_${i}_stripes_left_distalquant

        chromosight quantify --pattern stripes_left \
        --subsample ${mat_size} \
        --win-fmt npy \
        distal_${i}_stripes_left.tsv \
        proximal.${i}_balanced.cool \
        distal_${i}_stripes_left_proximalquant

        chromosight quantify --pattern stripes_left \
        --subsample ${mat_size} \
        --win-fmt npy \
        proximal_${i}_stripes_left.tsv \
        distal.${i}_balanced.cool \
        proximal_${i}_stripes_left_distalquant

        chromosight quantify --pattern stripes_left \
        --subsample ${mat_size} \
        --win-fmt npy \
        proximal_${i}_stripes_left.tsv \
        proximal.${i}_balanced.cool \
        proximal_${i}_stripes_left_proximalquant


  	chromosight quantify --pattern stripes_right \
        --subsample ${mat_size} \
        --win-fmt npy \
        distal_${i}_stripes_right.tsv \
        distal.${i}_balanced.cool \
        distal_${i}_stripes_right_distalquant

        chromosight quantify --pattern stripes_right \
        --subsample ${mat_size} \
        --win-fmt npy \
        distal_${i}_stripes_right.tsv \
        proximal.${i}_balanced.cool \
        distal_${i}_stripes_right_proximalquant

        chromosight quantify --pattern stripes_right \
        --subsample ${mat_size} \
        --win-fmt npy \
        proximal_${i}_stripes_right.tsv \
        distal.${i}_balanced.cool \
        proximal_${i}_stripes_right_distalquant

        chromosight quantify --pattern stripes_right \
        --subsample ${mat_size} \
        --win-fmt npy \
        proximal_${i}_stripes_right.tsv \
        proximal.${i}_balanced.cool \
        proximal_${i}_stripes_right_proximalquant
done
