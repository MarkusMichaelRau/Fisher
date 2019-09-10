
for dir in dir*clout;
do 
    #echo $dir
	paste $dir/*para_bias* > $dir/compl_bias.dat
	awk 'NR>1{print}' $dir/compl_bias.dat > $dir/compl_bias_red.dat
	mv $dir/compl_bias.dat $dir/compl_bias_old.dat
	mv $dir/compl_bias_red.dat $dir/compl_bias.dat
done

