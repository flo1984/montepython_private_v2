for i in {1..30}
do
echo "Iteration: $i \n"
    mpirun -np 8 python montepython/MontePython.py run -o ../output/$1/ -f $3 -b ../output/$1/$1.bestfit -N $2 -c ../output/run$1/$1.covmat --silent > ../output/$1/term.txt

    python montepython/MontePython.py info ../output/$1 --noplot --silent --minimal
done