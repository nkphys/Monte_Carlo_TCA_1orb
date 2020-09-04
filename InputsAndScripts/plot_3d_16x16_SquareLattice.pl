
#vertical
#Lm1=15;arrow_offset=$(echo "${Lm1}-14" | bc -l);for i in {0..15};do arrow=$(echo "${arrow_offset}+${i}" | bc -l);x1=$(echo "0.0*${i}" | bc -l);y1=$(echo "${i}" | bc -l);x2=$(echo "${x1}+${Lm1}" | bc -l);y2=${y1};echo "set arrow ${arrow} from ${y1},${x1} to ${y2},${x2} lw 0.75 lc black dt 2 nohead" ;done
set arrow 1 from 0,0 to 0,15 lw 0.75 lc black dt 2 nohead
set arrow 2 from 1,0 to 1,15 lw 0.75 lc black dt 2 nohead
set arrow 3 from 2,0 to 2,15 lw 0.75 lc black dt 2 nohead
set arrow 4 from 3,0 to 3,15 lw 0.75 lc black dt 2 nohead
set arrow 5 from 4,0 to 4,15 lw 0.75 lc black dt 2 nohead
set arrow 6 from 5,0 to 5,15 lw 0.75 lc black dt 2 nohead
set arrow 7 from 6,0 to 6,15 lw 0.75 lc black dt 2 nohead
set arrow 8 from 7,0 to 7,15 lw 0.75 lc black dt 2 nohead
set arrow 9 from 8,0 to 8,15 lw 0.75 lc black dt 2 nohead
set arrow 10 from 9,0 to 9,15 lw 0.75 lc black dt 2 nohead
set arrow 11 from 10,0 to 10,15 lw 0.75 lc black dt 2 nohead
set arrow 12 from 11,0 to 11,15 lw 0.75 lc black dt 2 nohead
set arrow 13 from 12,0 to 12,15 lw 0.75 lc black dt 2 nohead
set arrow 14 from 13,0 to 13,15 lw 0.75 lc black dt 2 nohead
set arrow 15 from 14,0 to 14,15 lw 0.75 lc black dt 2 nohead
set arrow 16 from 15,0 to 15,15 lw 0.75 lc black dt 2 nohead


# horizontal
#Lm1=15;arrow_offset=$(echo "${Lm1}+2" | bc -l);for i in {0..15};do arrow=$(echo "${arrow_offset}+${i}" | bc -l);x1=$(echo "0.0*${i}" | bc -l);y1=$(echo "${i}" | bc -l);x2=$(echo "${x1}+${Lm1}" | bc -l);y2=${y1};echo "set arrow ${arrow} from ${x1},${y1} to ${x2},${y2} lw 0.75 lc black dt 2 nohead" ;done
set arrow 17 from 0,0 to 15,0 lw 0.75 lc black dt 2 nohead
set arrow 18 from 0,1 to 15,1 lw 0.75 lc black dt 2 nohead
set arrow 19 from 0,2 to 15,2 lw 0.75 lc black dt 2 nohead
set arrow 20 from 0,3 to 15,3 lw 0.75 lc black dt 2 nohead
set arrow 21 from 0,4 to 15,4 lw 0.75 lc black dt 2 nohead
set arrow 22 from 0,5 to 15,5 lw 0.75 lc black dt 2 nohead
set arrow 23 from 0,6 to 15,6 lw 0.75 lc black dt 2 nohead
set arrow 24 from 0,7 to 15,7 lw 0.75 lc black dt 2 nohead
set arrow 25 from 0,8 to 15,8 lw 0.75 lc black dt 2 nohead
set arrow 26 from 0,9 to 15,9 lw 0.75 lc black dt 2 nohead
set arrow 27 from 0,10 to 15,10 lw 0.75 lc black dt 2 nohead
set arrow 28 from 0,11 to 15,11 lw 0.75 lc black dt 2 nohead
set arrow 29 from 0,12 to 15,12 lw 0.75 lc black dt 2 nohead
set arrow 30 from 0,13 to 15,13 lw 0.75 lc black dt 2 nohead
set arrow 31 from 0,14 to 15,14 lw 0.75 lc black dt 2 nohead
set arrow 32 from 0,15 to 15,15 lw 0.75 lc black dt 2 nohead


set xr [-2:18]
set yr [-2:18]
set zr [-1:1]

set ticslevel 0


set palette define (0 "red", 1 "blue")
unset cbr
#rgb(r,g,b) = r
#int(r)*65536 + int(g)*256 + int (b)
#set cbr [-1:-1]
#set palette rgbformulae 33,13,10

sp "/home/nitin/Documents/Codes/Monte_Carlo_TCA_1orb/ThetaPhi_Temp0.0001000000MicroState_given.txt" u ($1):($2):(0):($5*sin($3)*cos($4) ):($5*sin($3)*sin($4)):(($5*cos($3))) w vec head size 0.2,30,60 filled lw 2 linecolor palette


#sp "/home/nitin/Desktop/Telperion/data/home/n01/Triangular_Lattice_MCMF_Aug17_2020/ZigZagAlongZ_12x12_SIA0.5/U_8.0/ThetaPhi_Temp0.0001000000MicroState0.txt" u ($1-($2*0.5)):($2*(sqrt(3)/2.0)):(0):(2*$5*sin($3)*cos($4) ):(2*$5*sin($3)*sin($4)):((2*$5*cos($3))) w vec head size 0.2,30,60 filled lw 2 linecolor palette


#sp "/home/nitin/Desktop/Telperion/data/home/n01/Triangular_Lattice_MCMF_Aug17_2020/Random_Ising_theta_and_homogeneous_m_A1.0_12x12_MP/U_10.0/ThetaPhi_Temp0.0001000000MicroState0.txt" u ($1-($2*0.5)):($2*(sqrt(3)/2.0)):(0):(5*$5*sin($3)*cos($4) ):(5*$5*sin($3)*sin($4)):((5*$5*cos($3))) w vec head size 0.2,30,60 filled lw 2 linecolor palette

#sp "/home/nitin/Desktop/Telperion/data/home/n01/Triangular_Lattice_MCMF_Aug17_2020/Random_Theta_and_Phi_A0.5_12x12_Cooling/U_10.0/ThetaPhi_Temp0.0010000000MicroState0.txt" u ($1-($2*0.5)):($2*(sqrt(3)/2.0)):(0):(2*$5*sin($3)*cos($4) ):(2*$5*sin($3)*sin($4)):((2*$5*cos($3))) w vec head size 0.2,30,60 filled lw 2 linecolor palette

#sp "/home/nitin/Desktop/Telperion/data/home/n01/Triangular_Lattice_MCMF_Aug17_2020/Diagonal_ZigZag_Ising_alongZ_A0.5_12x12/U_8.0/ThetaPhi_Temp0.0001000000MicroState0.txt" u ($1-($2*0.5)):($2*(sqrt(3)/2.0)):(0):(2*$5*sin($3)*cos($4) ):(2*$5*sin($3)*sin($4)):((2*$5*cos($3))) w vec head size 0.2,30,60 filled lw 2 linecolor palette

#sp "/home/nitin/Desktop/Telperion/data/home/n01/Triangular_Lattice_MCMF_Aug17_2020/ZigZagAlongZ_12x12_SIA0.5/U_5.0/ThetaPhi_Temp0.0001000000MicroState0.txt" u ($1-($2*0.5)):($2*(sqrt(3)/2.0)):(0):(2*$5*sin($3)*cos($4) ):(2*$5*sin($3)*sin($4)):((2*$5*cos($3))) w vec head size 0.2,30,60 filled lw 2 linecolor palette

#sp "/home/nitin/Desktop/Telperion/data/home/n01/Triangular_Lattice_MCMF_Aug17_2020/Random_Theta_and_Phi_A0.5_8x8/U_6.0/ThetaPhi_Temp0.0001000000MicroState0.txt" u ($1-($2*0.5)):($2*(sqrt(3)/2.0)):(0):(($5*cos($3))):($5*sin($3)*cos($4) ):($5*sin($3)*sin($4)) w vec head size 0.2,30,60 filled  lw 5 #linecolor palette

#sp "/home/nitin/Desktop/Telperion/data/home/n01/Triangular_Lattice_MCMF_Aug17_2020/Coplanar120_XZPlane_8x8_SIA0.5/U_6.0/ThetaPhi_Temp0.0001000000MicroState0.txt" u ($1-($2*0.5)):($2*(sqrt(3)/2.0)):(0):((2*$5*cos($3))):(2*$5*sin($3)*cos($4) ):(2*$5*sin($3)*sin($4)) w vec head size 0.2,30,60 filled  lw 5 #linecolor palette

#p "/home/nitin/Documents/Codes/Monte_Carlo_TCA_1orb/ThetaPhi_Temp0.0001000000MicroState0.txt" u ($1-($2*0.5)):($2*(sqrt(3)/2.0)):($5*cos($3)*1):($5*sin($3)*cos($4)*1) w vec lw 5 lc "blue"

#p "/home/nitin/Desktop/Telperion/data/home/n01/Triangular_Lattice_MCMF_Aug17_2020/Random_justPhi_A0.5_8x8/U_8.5/ThetaPhi_Temp0.0001000000MicroState0.txt" u ($1-($2*0.5)):($2*(sqrt(3)/2.0)):($5*cos($4)*1 ):($5*sin($4)*1) w vec lw 5 lc "blue"
