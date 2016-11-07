begin=$(date +"%s")
codepath=/home/tang/Desktop/KernelCut_ECCV16/
executable=$codepath/kernelcut/main
exppath=$codepath/motionsegmentation/fencing
outdir=$codepath/motionsegmentation/fencing/output

# weight of smoothness term
lambda="1"
for ((frameid=0; frameid<3; frameid++))
do
    frame=`echo "scale=2;$frameid+40"|bc`
    preframe=`echo "scale=2;$frameid+39"|bc`
    frame="0$frame"
    preframe="0$preframe"
    echo "frame = $frame"
    if [ "$frameid" -eq "0" ]
    then
        userinput="-u seeds"
        hardconstraintsflag="on"
    else
        userinput="-u fromimage $outdir/ft_${preframe}_ncutknnmulti_s$lambda.bmp"
        hardconstraintsflag="off"
    fi
    argv="-d $exppath -i ft_$frame -n 3 -e nomeasure -h $hardconstraintsflag -o $outdir -k 50 $exppath/knn -s $lambda $userinput"
    cmd="$executable $argv"
    #echo $cmd
    mycmd=`$cmd`
    echo -e "$mycmd"
done

termin=$(date +"%s")
difftimelps=$(($termin-$begin))
echo "$(($difftimelps / 60)) minutes and $(($difftimelps % 60)) seconds elapsed for Script Execution."

