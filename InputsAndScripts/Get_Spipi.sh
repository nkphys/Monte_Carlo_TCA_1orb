
file="U14.0_spipi.txt"
rm ${file}
echo "#Temp  spipi  mpipi   <s^2>    m^2" > ${file}
for temp in 0.0200 0.0400 0.0600 0.0800 0.1000 0.1200 0.1400 0.1600 0.1800 0.2000 0.2200 0.2400 0.2600 0.2800 0.3000 0.3200 0.3400 0.3600 0.3800 0.4000 0.4200 0.4400 0.4600 0.4800 0.5000
do

val=$(grep "3.14159        3.14159" quantum_momentum_space_corr${temp}.txt | awk '{print $3}')
val2=$(grep "3.14159        3.14159" classical_momentum_space_corr${temp}.txt | awk '{print $3}')

val3=$(grep "0              0 " quantum_real_space_corr${temp}.txt | awk '{print $3}')
val4=$(grep "0              0 " classical_real_space_corr${temp}.txt | awk '{print $3}')

echo "${temp} ${val} ${val2}  ${val3}  ${val4}" >> ${file}
done
