#!/bin/csh

#$ -M 	brit004@ucr.edu # Email address for job notification
#$ -m  abe		 # Send mail when job begins, ends and aborts
#$ -q  gpu 	 # Specify queue
#$ -l gpu_card=1
#s -pe smp 4         #specifies threads??? maybe
#$ -N  Test	 # Specify job name
#$ -t 1       #specify number of data input files

set data = ( data_Rec_0.1_3_25_25_25_1.xml )


#( data_Rec_0.013_3_20_20_20_0.05_5.xml data_Rec_0.013_3_20_20_20_0.05_9.xml data_Rec_0.0653_3_20_20_20_0.05_5.xml data_Rec_0.0653_3_20_20_20_0.05_9.xml data_Rec_0.1177_3_20_20_20_0.05_5.xml data_Rec_0.1177_3_20_20_20_0.05_9.xml data_Rec_0.17_3_20_20_20_0.05_5.xml data_Rec_0.17_3_20_20_20_0.05_9.xml )

#( data_Rec_0.013_3_30_30_30_0.05_5.xml data_Rec_0.013_3_30_30_30_0.05_9.xml data_Rec_0.0653_3_30_30_30_0.05_5.xml data_Rec_0.0653_3_30_30_30_0.05_9.xml data_Rec_0.1177_3_30_30_30_0.05_5.xml data_Rec_0.1177_3_30_30_30_0.05_9.xml data_Rec_0.17_3_30_30_30_0.05_5.xml data_Rec_0.17_3_30_30_30_0.05_9.xml )

#( data_Rec_0.013_3_40_40_40_0.05_5.xml data_Rec_0.013_3_40_40_40_0.05_9.xml data_Rec_0.0653_3_40_40_40_0.05_5.xml data_Rec_0.0653_3_40_40_40_0.05_9.xml data_Rec_0.1177_3_40_40_40_0.05_5.xml data_Rec_0.1177_3_40_40_40_0.05_9.xml data_Rec_0.17_3_40_40_40_0.05_5.xml data_Rec_0.17_3_40_40_40_0.05_9.xml )

module load gcc/6.2.0
module load gsl/2.3
module load cuda/8.0
#module load boost/1.61
echo -n "It is currently: ";date
echo -n "I am logged on as ";whoami
echo -n "This computer is called ";hostname
echo -n "I am currently in the directory ";pwd


./bend-model diagram -eps=0.0001 -df=1.0 -dt=0.001 --lagTime=0.1 --pullPercent=0.2 $data[${SGE_TASK_ID}]

#-dataField=3

#./bend-model diagram -eps=0.0001 -df=1 -dt=0.01 --nodecoords data_Hem_1_8_0.15_25_1000_Origin.xml
#./bend-model diagram -eps=0.00001 -df=1 -dt=0.01 --nodecoords calibration.xml
