begin=$(date +"%s")
codepath=/home/meng/Desktop/repositories/KernelCut_ECCV16/
executable=$codepath/kernelcut/main
exppath=$codepath/motionsegmentation/ducks01
outdir=$codepath/motionsegmentation/ducks01/output

# weight of smoothness term
lambda="0.5"
for ((frameid=0; frameid<41; frameid++))
do
    frame=`echo "scale=2;$frameid+300"|bc`
    preframe=`echo "scale=2;$frameid+299"|bc`
    frame="0$frame"
    preframe="0$preframe"
    echo "frame = $frame"
    if [ "$frameid" -eq "0" ]
    then
        userinput="-u seeds"
    else
        userinput="-u fromimage $outdir/ducks01_${preframe}_ncutknnmulti_s$lambda.bmp"
    fi
    argv="-d $exppath -i ducks01_$frame -n 6 -e nomeasure -h off -o $outdir -k 50 $exppath/knn -s $lambda $userinput"
    cmd="$executable $argv"
    #echo $cmd
    mycmd=`$cmd`
    echo -e "$mycmd"
done

termin=$(date +"%s")
difftimelps=$(($termin-$begin))
echo "$(($difftimelps / 60)) minutes and $(($difftimelps % 60)) seconds elapsed for Script Execution."

