

#!/bin/bash
rm output.txt
number_of_points=7
for (( n=0 ; n<=$number_of_points ; n++ )); 
do
    ./diver.example_c $n >> "output.txt"
done


#chmod +x diver_looper.sh
