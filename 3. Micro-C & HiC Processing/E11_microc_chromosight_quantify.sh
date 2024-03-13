#! /bin/bash

###
### chromosight quantification for loops and strips
###

module load python/3.11.2

res=(1000000 500000 250000 100000 50000 40000 25000 20000 10000 5000 2500 1000 800 600 400)

for i in "${res[@]}"
do
        clamped_size=$(cooler info clamped.${i}_balanced.cool | grep sum | tr -dc '0-9')
	unclamped_size=$(cooler info unclamped.${i}_balanced.cool | grep sum | tr -dc '0-9')

	if [ ${clamped_size} -gt ${unclamped_size} ];
	then
		mat_size=${unclamped_size}
	else
		mat_size=${clamped_size}
	fi
	echo "CC size = ${clamped_size}"
	echo "UC size = ${unclamped_size}"
	echo "mat subsample size = ${mat_size}"


        chromosight quantify --pattern loops \
        --subsample ${mat_size} \
        --win-fmt npy \
        clamped_${i}_loops.tsv \
        clamped.${i}_balanced.cool \
        clamped_${i}_loops_clampedquant

        chromosight quantify --pattern loops \
        --subsample ${mat_size} \
        --win-fmt npy \
        clamped_${i}_loops.tsv \
        unclamped.${i}_balanced.cool \
        clamped_${i}_loops_unclampedquant

        chromosight quantify --pattern loops \
        --subsample ${mat_size} \
        --win-fmt npy \
        unclamped_${i}_loops.tsv \
        clamped.${i}_balanced.cool \
        unclamped_${i}_loops_clampedquant

        chromosight quantify --pattern loops \
        --subsample ${mat_size} \
        --win-fmt npy \
        unclamped_${i}_loops.tsv \
        unclamped.${i}_balanced.cool \
        unclamped_${i}_loops_unclampedquant


  	chromosight quantify --pattern stripes_left \
        --subsample ${mat_size} \
        --win-fmt npy \
        clamped_${i}_stripes_left.tsv \
        clamped.${i}_balanced.cool \
        clamped_${i}_stripes_left_clampedquant

        chromosight quantify --pattern stripes_left \
        --subsample ${mat_size} \
        --win-fmt npy \
        clamped_${i}_stripes_left.tsv \
        unclamped.${i}_balanced.cool \
        clamped_${i}_stripes_left_unclampedquant

        chromosight quantify --pattern stripes_left \
        --subsample ${mat_size} \
        --win-fmt npy \
        unclamped_${i}_stripes_left.tsv \
        clamped.${i}_balanced.cool \
        unclamped_${i}_stripes_left_clampedquant

        chromosight quantify --pattern stripes_left \
        --subsample ${mat_size} \
        --win-fmt npy \
        unclamped_${i}_stripes_left.tsv \
        unclamped.${i}_balanced.cool \
        unclamped_${i}_stripes_left_unclampedquant


  	chromosight quantify --pattern stripes_right \
        --subsample ${mat_size} \
        --win-fmt npy \
        clamped_${i}_stripes_right.tsv \
        clamped.${i}_balanced.cool \
        clamped_${i}_stripes_right_clampedquant

        chromosight quantify --pattern stripes_right \
        --subsample ${mat_size} \
        --win-fmt npy \
        clamped_${i}_stripes_right.tsv \
        unclamped.${i}_balanced.cool \
        clamped_${i}_stripes_right_unclampedquant

        chromosight quantify --pattern stripes_right \
        --subsample ${mat_size} \
        --win-fmt npy \
        unclamped_${i}_stripes_right.tsv \
        clamped.${i}_balanced.cool \
        unclamped_${i}_stripes_right_clampedquant

        chromosight quantify --pattern stripes_right \
        --subsample ${mat_size} \
        --win-fmt npy \
        unclamped_${i}_stripes_right.tsv \
        unclamped.${i}_balanced.cool \
        unclamped_${i}_stripes_right_unclampedquant
done
