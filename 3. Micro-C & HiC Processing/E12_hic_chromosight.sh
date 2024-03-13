#! /bin/bash

###
# Execution code for E12.5 Hi-C chromosight loop and stripe calling
###

module load python/3.11.2
res=(1000000 500000 250000 100000 50000 40000 25000 20000 10000 5000 2500 1000 800 600 400)

for i in "${res[@]}"
do
	chromosight detect --pattern loops -t16 \
	distal_microc.${i}_balanced.cool \
	distal_${i}_loops
	
	chromosight detect --pattern loops -t16 \
        proximal.${i}_balanced.cool \
        proximal_${i}_loops

        chromosight detect --pattern stripes_left -t16 \
        distal.${i}_balanced.cool \
        distal_${i}_stripes_left

        chromosight detect --pattern stripes_left -t16 \
        proximal.${i}_balanced.cool \
        proximal_${i}_stripes_left

        chromosight detect --pattern stripes_right -t16 \
        distal.${i}_balanced.cool \
        distal_${i}_stripes_right

        chromosight detect --pattern stripes_right -t16 \
        proximal.${i}_balanced.cool \
        proximal_${i}_stripes_right
done
