#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <random>
#include <complex>
#include <cmath>
#include <cassert>
using namespace std;

#include "Matrix.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#include "Hamiltonian.h"
#include "Observables.h"
#include "MCEngine.h"
#include "HF_SC_Engine.h"
#include "random"


int main(int argc, char *argv[]) {
    if (argc<2) { throw std::invalid_argument("USE:: executable inputfile"); }
    string inputfile = argv[1];

    bool check_Non_Int=false;
    bool Momentum_space_calculations=false;



    Parameters Parameters_;
    Parameters_.Initialize(inputfile);

    Coordinates Coordinates_(Parameters_.lx, Parameters_.ly);
    Coordinates CoordinatesCluster_(Parameters_.lx_cluster, Parameters_.ly_cluster);

    mt19937_64 Generator_(Parameters_.RandomSeed); //for random fields
    mt19937_64 Generator2_(Parameters_.RandomDisorderSeed); //for random disorder

    MFParams MFParams_(Parameters_,Coordinates_,Generator_, Generator2_);

    Hamiltonian Hamiltonian_(Parameters_,Coordinates_,CoordinatesCluster_,MFParams_);


    Observables Observables_(Parameters_,Coordinates_,MFParams_,Hamiltonian_);


    Parameters_.Single_Diagonalization=false;


    if(check_Non_Int==true){

        //Parameters_.J_HUND=0.0;
        Hamiltonian_.InteractionsCreate();
        //  Hamiltonian_.Ham_.print();
        // Hamiltonian_.Check_up_down_symmetry();
        //Hamiltonian_.Check_Hermiticity();
        Hamiltonian_.Diagonalize('V');
        int temp=Parameters_.ns*Parameters_.Fill*2.0;
        cout<<"mu for n=4 = "<<0.5*(Hamiltonian_.eigs_[temp-1] + Hamiltonian_.eigs_[temp])<<"   "<<
             Hamiltonian_.eigs_[temp-1]<<"   "<<Hamiltonian_.eigs_[temp]<<endl;
        Parameters_.mus=0.5*(Hamiltonian_.eigs_[temp-1] + Hamiltonian_.eigs_[temp]);
        double Quantum_E=Hamiltonian_.E_QM();
        double Classical_E=Hamiltonian_.GetCLEnergy();
        cout<<setprecision(9);
        cout<<"Total_Energy = "<<Quantum_E+Classical_E<<endl;
        //double mu = chemicalpotential(0.5, temp);
        // Observables_.Get_Non_Interacting_dispersion();
        //Hamiltonian_.Ham_.print();
        //Observables_.Calculate_Akw();
        //Observables_.Calculate_Akw_at_w(mu);
        Observables_.Calculate_Nw();


    }
    else if(Momentum_space_calculations==true){

       //string Magnetic_state="Pi0_triangular_Lattice_Ising";
       string Magnetic_state="ZZ_vertical_triangular_Lattice_Ising";

        Hamiltonian_.K_space_Calculation(Magnetic_state);
        Observables_.Calculate_Nw();


    }

    else if(Parameters_.Perform_HF_SC_calculation==true){

        cout<<setprecision(9);
        HF_SC_Engine HF_SC_Engine_(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_);
        HF_SC_Engine_.RUN_HF_SC_Engine();


    }

    else if(Parameters_.Single_Diagonalization){


        if(Parameters_.Cooling_!=0){
           cout<<"Only Cooling=0 is allowed i.e. specific temperature"<<endl;
           assert(Parameters_.Cooling_==0);
        }

        double temp_ = Parameters_.Temp_values[0];

        Hamiltonian_.InteractionsCreate();
        Hamiltonian_.Diagonalize('V');
        int temp=Parameters_.ns*Parameters_.Fill*2.0;
        double mu_temp = 0.5*(Hamiltonian_.eigs_[temp-1] + Hamiltonian_.eigs_[temp]);
        cout<<"mu for T=0 : "<<mu_temp<<endl;
        Parameters_.mus=Hamiltonian_.chemicalpotential(mu_temp, Parameters_.Fill);;
        double Quantum_E=Hamiltonian_.E_QM();
        double Classical_E=Hamiltonian_.GetCLEnergy();
        cout<<setprecision(9);
        cout<<"Classical E = "<<Classical_E<<endl;
        cout<<"Quantum E = "<<Quantum_E<<endl;
        cout<<"Total_Energy = "<<Quantum_E+Classical_E<<endl;

        Observables_.calculate_local_density();


        int temp_site_;
        char temp_char[50];
        sprintf(temp_char, "%.10f", temp_);

        string File_Out_theta_phi_microState = "ThetaPhi_Temp" + string(temp_char) +
                "MicroState_given" + ".txt";
        ofstream File_Out_Theta_Phi_MicroState(File_Out_theta_phi_microState.c_str());

        File_Out_Theta_Phi_MicroState << "#x" << setw(15) << "y" << setw(15) << "Theta(x,y)" << setw(15) << "Phi(x,y)"
                                      << setw(15) << "Moment_Size(x,y)" << setw(15) << "Local_density(x,y)" << endl;
        for (int ix = 0; ix < Coordinates_.lx_; ix++)
        {
            for (int iy = 0; iy < Coordinates_.ly_; iy++)
            {
                File_Out_Theta_Phi_MicroState << ix << setw(15) << iy << setw(15) << MFParams_.etheta(ix, iy) << setw(15) << MFParams_.ephi(ix, iy)
                                              << setw(15) << MFParams_.Moment_Size(ix, iy) << setw(15) << MFParams_.Local_density(ix, iy) << endl;
            }
        }


        string File_Out_local_den_microState = "Local_den_Temp" + string(temp_char) +
                "MicroState_given" + ".txt";
        ofstream File_Out_local_den_MicroState(File_Out_local_den_microState.c_str());

        File_Out_local_den_MicroState << "#x" << setw(15) << "y" << setw(15)<< "<n_up>" << setw(15)<< "<n_dn>"<<endl;

        for (int ix = 0; ix < (Coordinates_.lx_+1); ix++)
        {
            for (int iy = 0; iy < (Coordinates_.ly_+1); iy++)
            {
                temp_site_ = Coordinates_.Nc(ix%Coordinates_.lx_, iy%Coordinates_.ly_);
                File_Out_local_den_MicroState << ix << setw(15) << iy << setw(15) << Observables_.local_density[temp_site_][0]<<setw(15)
                                              << Observables_.local_density[temp_site_][1] << endl;
            }
            File_Out_local_den_MicroState<<endl;
        }


    }

    else{

     cout<<setprecision(9);
     MCEngine MCEngine_(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_);

     Observables_.Initialize();     // Set All Observables to zero

     MCEngine_.RUN_MC();      // Monte-Carlo Engine

     Observables_.Calculate_Nw();

    }





    cout << "--------THE END--------" << endl;
} // main
