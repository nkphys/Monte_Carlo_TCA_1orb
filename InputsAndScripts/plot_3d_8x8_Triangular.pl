
# diagonal-1
# Lm1=7;arrow_offset=1;for i in {0..7};do arrow=$(echo "${arrow_offset}+${i}" | bc -l);x1=$(echo "${i}" | bc -l);y1=0;x2=$(echo "${x1} - (${Lm1}*0.5)" | bc -l);y2=$(echo "${Lm1}*sqrt(3)*0.5" | bc -l);echo "set arrow ${arrow} from ${x1},${y1} to ${x2},${y2} lw 2 lc "black" nohead" ;done

set arrow 1 from 0,0 to -3.5,6.06217782649107052732 lw 2 lc black nohead
set arrow 2 from 1,0 to -2.5,6.06217782649107052732 lw 2 lc black nohead
set arrow 3 from 2,0 to -1.5,6.06217782649107052732 lw 2 lc black nohead
set arrow 4 from 3,0 to -.5,6.06217782649107052732 lw 2 lc black nohead
set arrow 5 from 4,0 to .5,6.06217782649107052732 lw 2 lc black nohead
set arrow 6 from 5,0 to 1.5,6.06217782649107052732 lw 2 lc black nohead
set arrow 7 from 6,0 to 2.5,6.06217782649107052732 lw 2 lc black nohead
set arrow 8 from 7,0 to 3.5,6.06217782649107052732 lw 2 lc black nohead


# horizontal
#Lm1=7;arrow_offset=$(echo "${Lm1}+2" | bc -l);for i in {0..7};do arrow=$(echo "${arrow_offset}+${i}" | bc -l);x1=$(echo "-0.5*${i}" | bc -l);y1=$(echo "${i}*sqrt(3)*0.5" | bc -l);x2=$(echo "${x1}+${Lm1}" | bc -l);y2=${y1};echo "set arrow ${arrow} from ${x1},${y1} to ${x2},${y2} lw 2 lc "black" nohead" ;done
set arrow 9 from 0,0 to 7,0 lw 2 lc black nohead
set arrow 10 from -.5,.86602540378443864676 to 6.5,.86602540378443864676 lw 2 lc black nohead
set arrow 11 from -1.0,1.73205080756887729352 to 6.0,1.73205080756887729352 lw 2 lc black nohead
set arrow 12 from -1.5,2.59807621135331594028 to 5.5,2.59807621135331594028 lw 2 lc black nohead
set arrow 13 from -2.0,3.46410161513775458704 to 5.0,3.46410161513775458704 lw 2 lc black nohead
set arrow 14 from -2.5,4.33012701892219323380 to 4.5,4.33012701892219323380 lw 2 lc black nohead
set arrow 15 from -3.0,5.19615242270663188056 to 4.0,5.19615242270663188056 lw 2 lc black nohead
set arrow 16 from -3.5,6.06217782649107052732 to 3.5,6.06217782649107052732 lw 2 lc black nohead

#diagonal-2, set 1
#Lm1=7;arrow_offset=$(echo "((${Lm1} + 1)*2) + 1" | bc -l);for i in {0..6};do arrow=$(echo "${arrow_offset}+${i}" | bc -l);x1=$(echo "-${i}*0.5" | bc -l );y1=$(echo "${i}*sqrt(3)*0.5" | bc -l);x2=$(echo "(${Lm1}*0.5) - ${i}*1.0" | bc -l);y2=$( echo "${Lm1}*sqrt(3)*0.5" | bc -l);echo "set arrow ${arrow} from ${x1},${y1} to ${x2},${y2} lw 2 lc "black" nohead" ;done
set arrow 17 from 0,0 to 3.5,6.06217782649107052732 lw 2 lc black nohead
set arrow 18 from -.5,.86602540378443864676 to 2.5,6.06217782649107052732 lw 2 lc black nohead
set arrow 19 from -1.0,1.73205080756887729352 to 1.5,6.06217782649107052732 lw 2 lc black nohead
set arrow 20 from -1.5,2.59807621135331594028 to .5,6.06217782649107052732 lw 2 lc black nohead
set arrow 21 from -2.0,3.46410161513775458704 to -.5,6.06217782649107052732 lw 2 lc black nohead
set arrow 22 from -2.5,4.33012701892219323380 to -1.5,6.06217782649107052732 lw 2 lc black nohead
set arrow 23 from -3.0,5.19615242270663188056 to -2.5,6.06217782649107052732 lw 2 lc black nohead




#diagonal-2, set2
#Lm1=7;arrow_offset=$(echo "(${Lm1}+1)*3" | bc -l);for i in {0..5};do arrow=$(echo "${arrow_offset}+${i}" | bc -l);x1=$(echo "1+${i}" | bc -l );y1=0;x2=$(echo "(${Lm1}*0.5) + 0.5+ ${i}*0.5" | bc -l);y2=$( echo "(${Lm1} -1 - ${i})*(sqrt(3)*0.5)" | bc -l);echo "set arrow ${arrow} from ${x1},${y1} to ${x2},${y2} lw 2 lc "black" nohead" ;done
set arrow 24 from 1,0 to 4.0,5.19615242270663188056 lw 2 lc black nohead
set arrow 25 from 2,0 to 4.5,4.33012701892219323380 lw 2 lc black nohead
set arrow 26 from 3,0 to 5.0,3.46410161513775458704 lw 2 lc black nohead
set arrow 27 from 4,0 to 5.5,2.59807621135331594028 lw 2 lc black nohead
set arrow 28 from 5,0 to 6.0,1.73205080756887729352 lw 2 lc black nohead
set arrow 29 from 6,0 to 6.5,.86602540378443864676 lw 2 lc black nohead


set xr [-4:8]
set yr [-4:8]
set zr [-4:4]

set ticslevel 0


set palette define (0 "red", 1 "blue")

unset cbr
#set cbr [-1:-1]

#sp "/home/nitin/Desktop/Telperion/data/home/n01/Triangular_Lattice_MCMF_Aug17_2020/Random_Theta_and_Phi_A0.0_8x8/U_10.0/ThetaPhi_Temp0.0001000000MicroState0.txt" u ($1-($2*0.5)):($2*(sqrt(3)/2.0)):(0):(2*$5*sin($3)*cos($4) ):(2*$5*sin($3)*sin($4)):((2*$5*cos($3))) w vec head size 0.2,30,60 filled  lw 5 lc "blue" #linecolor palette z


#sp "/home/nitin/Desktop/Telperion/data/home/n01/Triangular_Lattice_MCMF_Aug17_2020/Random_Theta_and_Phi_A0.5_8x8/U_6.0/ThetaPhi_Temp0.0001000000MicroState0.txt" u ($1-($2*0.5)):($2*(sqrt(3)/2.0)):(0):(($5*cos($3))):($5*sin($3)*cos($4) ):($5*sin($3)*sin($4)) w vec head size 0.2,30,60 filled  lw 5 #linecolor palette

#sp "/home/nitin/Desktop/Telperion/data/home/n01/Triangular_Lattice_MCMF_Aug17_2020/Coplanar120_XZPlane_8x8_SIA0.5/U_6.0/ThetaPhi_Temp0.0001000000MicroState0.txt" u ($1-($2*0.5)):($2*(sqrt(3)/2.0)):(0):((2*$5*cos($3))):(2*$5*sin($3)*cos($4) ):(2*$5*sin($3)*sin($4)) w vec head size 0.2,30,60 filled  lw 5 #linecolor palette

sp "/home/nitin/Documents/Codes/Monte_Carlo_TCA_1orb/ThetaPhi_Temp0.0001000000MicroState0.txt" u ($1-($2*0.5)):($2*(sqrt(3)/2.0)):(0):(2*$5*sin($3)*cos($4) ):(2*$5*sin($3)*sin($4)):((2*$5*cos($3))) w vec head size 0.2,30,60 filled lw 2 linecolor palette

#p "/home/nitin/Desktop/Telperion/data/home/n01/Triangular_Lattice_MCMF_Aug17_2020/Random_justPhi_A0.5_8x8/U_8.5/ThetaPhi_Temp0.0001000000MicroState0.txt" u ($1-($2*0.5)):($2*(sqrt(3)/2.0)):($5*cos($4)*1 ):($5*sin($4)*1) w vec lw 5 lc "blue"
