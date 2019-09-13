for i in {1..38}
do
	##myjob='MergeChr'$i
	##qsub -cwd -V -N ${myjob} -l highp,h_data=4G,time=06:00:00 -pe shared 2 -M eplau -m ea merge_dogs.sh ${i}
	
	 ##myjob='FilterChr'$i
	 ##qsub -cwd -V -N ${myjob} -l highp,h_data=10G,time=00:20:00 -pe shared 2 -M eplau -m a Filter_AN120_biallelic.sh ${i}
	

	sleep 5
done
