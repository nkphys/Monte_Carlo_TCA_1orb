#include <algorithm>
#include <functional>
#include <math.h>
#include "tensor_type.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#define PI acos(-1.0)

#ifndef Hamiltonian_class
#define Hamiltonian_class

extern "C" void zheev_(char *, char *, int *, std::complex<double> *, int *, double *,
                       std::complex<double> *, int *, double *, int *);

class Hamiltonian
{
public:
    Hamiltonian(Parameters &Parameters__, Coordinates &Coordinates__, Coordinates &CoordinatesCluster__, MFParams &MFParams__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__), CoordinatesCluster_(CoordinatesCluster__), MFParams_(MFParams__)

    {
        Initialize();
        Hoppings();
        HTBCreate();
        HTBClusterCreate();
    }

    void Initialize();                                     //::DONE
    void Hoppings();                                       //::DONE
    double GetCLEnergy();                                  //::DONE
    void InteractionsCreate();                             //::DONE
    void InteractionsClusterCreate(int Center_site);       //::DONE
    void Check_Hermiticity();                              //::DONE
    void Check_up_down_symmetry();                         //::DONE
    void HTBCreate();                                      //::DONE
    void HTBClusterCreate();                               //::DONE
    double chemicalpotential(double muin, double filling); //::DONE

    double chemicalpotentialCluster(double muin, double filling); //::DONE

    double TotalDensity();   //::DONE
    double ClusterDensity(); //::DONE
    double E_QM();           //::DONE

    double E_QMCluster();                 //::DONE
    void Diagonalize(char option);        //::DONE
    void DiagonalizeCluster(char option); //::DONE
    void copy_eigs(int i);                //::DONE
    void copy_eigs_Cluster(int i);        //::DONE

    void clone(Hamiltonian &Hamiltonian_new);

    void K_space_Calculation(string Magnetic_state);
    int partition(vector<double> &a, int s, int e);
    void quicksort(vector<double> &a, int s, int e);


    Parameters &Parameters_;
    Coordinates &Coordinates_;
    Coordinates &CoordinatesCluster_;
    MFParams &MFParams_;
    int lx_, ly_, ns_, orbs_;
    Matrix<complex<double>> HTB_;
    Matrix<complex<double>> HTBCluster_;
    Matrix<complex<double>> Ham_;
    Matrix<complex<double>> HamCluster_;
    Matrix<complex<double>> HamCluster_saved_, Ham_saved_;
    Matrix<double> Tx, Ty, Tpxpy, Tpxmy;
    vector<double> eigs_, eigsCluster_, eigsCluster_saved_, eigs_saved_, sx_, sy_, sz_;

    double HS_factor;
};



int Hamiltonian::partition(vector<double> &a, int s, int e)
{
    double piviot = a[e];
    int pind = s;
    int i;
    double t;

    for (i = s; i < e; i++) {
        if (a[i] <= piviot) {
            t = a[i];
            a[i] = a[pind];
            a[pind] = t;
            pind++;
        }
    }

    t = a[e];
    a[e] = a[pind];
    a[pind] = t;

    return pind;
}

void Hamiltonian::quicksort(vector<double> &a, int s, int e)
{
    if (s < e) {
        int pind = partition(a, s, e);
        quicksort(a, s, pind - 1);
        quicksort(a, pind + 1, e);
    }
}

void Hamiltonian::K_space_Calculation(string Magnetic_state){



    Matrix<complex<double>>  Ham_up, Ham_dn;
    double U_=-0.5*Parameters_.J_Hund;
    double K_;
    K_=  2.0*(U_ + Parameters_.SIA);
    double M_; //Moment
    M_=Parameters_.Kspace_Moment_fixed;
    double kx_, ky_;
    vector<double> eigs_up, eigs_dn;
    vector<double> eigs_all;


    if(Magnetic_state=="Pi0_triangular_Lattice_Ising"){

        string eigs_out_file = "Eigenvalues_"+Magnetic_state+".txt";
        ofstream eigs_out_fl(eigs_out_file.c_str());
        eigs_out_fl<<"#nx   ny    E_up[0]  E_up[1]  E_dn[0]  E_dn[1]"<<endl;

        for(int nx=0;nx<Parameters_.lx;nx++){
            kx_=(2.0*nx*PI)/(Parameters_.lx);
            for(int ny=0;ny<Parameters_.ly;ny++){
                ky_=(2.0*ny*PI)/(Parameters_.ly);

                Ham_up.resize(2,2);
                Ham_up(0,0)= -0.5*K_*M_  + 0.5*U_ + 2.0*Parameters_.t_hopping*cos(ky_);
                Ham_up(1,1)= 0.5*K_*M_ + + 0.5*U_ + 2.0*Parameters_.t_hopping*cos(ky_);
                Ham_up(0,1)= Parameters_.t_hopping*(one_complex + exp(iota_complex*kx_)
                                                    +exp(-1.0*iota_complex*ky_)+exp(iota_complex*(kx_+ky_))
                                                    );
                Ham_up(1,0)=conj(Ham_up(0,1));

                Ham_dn.resize(2,2);
                Ham_dn(0,0)= 0.5*K_*M_ + 0.5*U_  +2.0*Parameters_.t_hopping*cos(ky_);
                Ham_dn(1,1)= -0.5*K_*M_ + 0.5*U_  +2.0*Parameters_.t_hopping*cos(ky_);
                Ham_dn(0,1)= Parameters_.t_hopping*(one_complex + exp(iota_complex*kx_)
                                                    +exp(-1.0*iota_complex*ky_)+exp(iota_complex*(kx_+ky_))
                                                    );
                Ham_dn(1,0)=conj(Ham_dn(0,1));


                Ham_=Ham_up;
                Diagonalize('V');
                eigs_up=eigs_;

                Ham_=Ham_dn;
                Diagonalize('V');
                eigs_dn=eigs_;

                eigs_out_fl<<nx<<"  "<<ny<<"  "<<eigs_up[0]<<"  "<<eigs_up[1]<<"  "<<eigs_dn[0]<<"  "<<eigs_dn[1]<<endl;


                for(int n=0;n<eigs_up.size();n++){
                    eigs_all.push_back(eigs_up[n]);
                    eigs_all.push_back(eigs_dn[n]);
                }


            }
            eigs_out_fl<<endl;
        }
        eigs_=eigs_all;
        quicksort(eigs_, 0, eigs_.size()-1);


        int n_states_occupied_zeroT = Parameters_.ns * Parameters_.Fill * 2.0;
        double initial_mu_guess = 0.5 * (eigs_[n_states_occupied_zeroT - 1] + eigs_[n_states_occupied_zeroT]);
        Parameters_.mus = chemicalpotential(initial_mu_guess, Parameters_.Fill);

        double Energy_QM, Energy_Total;
        Energy_QM=E_QM();

        Energy_Total = Energy_QM + 2.0*Parameters_.lx*Parameters_.ly*(((U_+Parameters_.SIA)*M_*M_)  -  (U_*0.25));

        cout<<setprecision(10);
        cout<<"mu = "<<Parameters_.mus<<endl;
        cout<<"Energy_Total = "<<Energy_Total<<endl;
        cout<<"Energy_Classical = "<<2.0*Parameters_.lx*Parameters_.ly*(((U_+Parameters_.SIA)*M_*M_)  -  (U_*0.25))<<endl;
    }

    if(Magnetic_state=="ZZ_vertical_triangular_Lattice_Ising"){
        string eigs_out_file = "Eigenvalues_"+Magnetic_state+".txt";
        ofstream eigs_out_fl(eigs_out_file.c_str());
        eigs_out_fl<<"#nx   ny    E_up[0]  E_up[1]  E_up[2]  E_up[3]  E_dn[0]  E_dn[1]  E_dn[2]  E_dn[3]"<<endl;

        for(int nx=0;nx<Parameters_.lx;nx++){
            kx_=(2.0*nx*PI)/(Parameters_.lx);
            for(int ny=0;ny<Parameters_.ly;ny++){
                ky_=(2.0*ny*PI)/(Parameters_.ly);

                Ham_up.resize(8,8);
                Ham_dn.resize(8,8);

                Ham_up(0,0)= -0.5*K_*M_ + 0.5*U_;
                Ham_up(1,1)= 0.5*K_*M_ + 0.5*U_ ;
                Ham_up(2,2)= -0.5*K_*M_+ 0.5*U_;
                Ham_up(3,3)= 0.5*K_*M_+ 0.5*U_;
                Ham_up(4,4)= 0.5*K_*M_ + 0.5*U_;
                Ham_up(5,5)= -0.5*K_*M_ + 0.5*U_ ;
                Ham_up(6,6)= 0.5*K_*M_+ 0.5*U_;
                Ham_up(7,7)= -0.5*K_*M_+ 0.5*U_;


                /*
                Ham_up(0,1)= Parameters_.t_hopping*(one_complex + exp(iota_complex*kx_));
                Ham_up(0,2)= Parameters_.t_hopping*(one_complex + exp(iota_complex*ky_));
                Ham_up(0,3)= Parameters_.t_hopping*(one_complex + exp(iota_complex*ky_));
                Ham_up(1,2)= Parameters_.t_hopping*(exp(-1.0*iota_complex*kx_))*(one_complex + exp(iota_complex*ky_));
                Ham_up(1,3)= Parameters_.t_hopping*(one_complex + exp(iota_complex*ky_));
                Ham_up(2,3)= Parameters_.t_hopping*(one_complex + exp(iota_complex*kx_));
                Ham_up(1,0)=conj(Ham_up(0,1));
                Ham_up(2,0)=conj(Ham_up(0,2));
                Ham_up(3,0)=conj(Ham_up(0,3));
                Ham_up(2,1)=conj(Ham_up(1,2));
                Ham_up(3,1)=conj(Ham_up(1,3));
                Ham_up(3,2)=conj(Ham_up(2,3));
                */

                Ham_up(0,1)= Parameters_.t_hopping*(one_complex + exp(iota_complex*kx_));
                Ham_up(0,2)= Parameters_.t_hopping*(one_complex);
                Ham_up(0,3)= Parameters_.t_hopping*(one_complex);
                Ham_up(0,6)= Parameters_.t_hopping*(exp(iota_complex*ky_));
                Ham_up(0,7)= Parameters_.t_hopping*(exp(iota_complex*(ky_+kx_)));
                Ham_up(1,2)= Parameters_.t_hopping*(exp(-1.0*iota_complex*(kx_)));
                Ham_up(1,3)= Parameters_.t_hopping*(one_complex);
                Ham_up(1,6)= Parameters_.t_hopping*(exp(iota_complex*ky_));
                Ham_up(1,7)= Parameters_.t_hopping*(exp(iota_complex*ky_));
                Ham_up(2,3)= Parameters_.t_hopping*(one_complex + exp(iota_complex*kx_));
                Ham_up(2,4)= Parameters_.t_hopping*(one_complex);
                Ham_up(2,5)= Parameters_.t_hopping*(one_complex);
                Ham_up(3,4)= Parameters_.t_hopping*(exp(-1.0*iota_complex*kx_));
                Ham_up(3,5)= Parameters_.t_hopping;
                Ham_up(4,5)= Parameters_.t_hopping*(one_complex + exp(iota_complex*kx_));
                Ham_up(4,6)= Parameters_.t_hopping;
                Ham_up(4,7)= Parameters_.t_hopping;
                Ham_up(5,6)= Parameters_.t_hopping*(exp(-1.0*iota_complex*kx_));
                Ham_up(5,7)= Parameters_.t_hopping;
                Ham_up(6,7)= Parameters_.t_hopping*(one_complex + exp(iota_complex*kx_));
                for(int i=0;i<8;i++){
                    for(int j=i+1;j<8;j++){
                        Ham_up(j,i)=conj(Ham_up(i,j));
                    }
                }


                Ham_dn=Ham_up;
                Ham_dn(0,0)= 0.5*K_*M_ + 0.5*U_;
                Ham_dn(1,1)= -0.5*K_*M_ + 0.5*U_ ;
                Ham_dn(2,2)= 0.5*K_*M_+ 0.5*U_;
                Ham_dn(3,3)= -0.5*K_*M_+ 0.5*U_;
                Ham_dn(4,4)= -0.5*K_*M_ + 0.5*U_;
                Ham_dn(5,5)= 0.5*K_*M_ + 0.5*U_ ;
                Ham_dn(6,6)= -0.5*K_*M_+ 0.5*U_;
                Ham_dn(7,7)= 0.5*K_*M_+ 0.5*U_;



                Ham_=Ham_up;
                Diagonalize('V');
                eigs_up=eigs_;

                Ham_=Ham_dn;
                Diagonalize('V');
                eigs_dn=eigs_;

                eigs_out_fl<<nx<<"  "<<ny<<"  "<<eigs_up[0]<<"  "<<eigs_up[1]<<"  "<<eigs_up[2]<<"  "<<eigs_up[3]<<"  "
                          <<eigs_up[4]<<"  "<<eigs_up[5]<<"  "<<eigs_up[6]<<"  "<<eigs_up[7]<<"  "
                          <<eigs_dn[0]<<"  "<<eigs_dn[1]<<"   "<<eigs_dn[2]<<"  "<<eigs_dn[3]<<"  "
                         <<eigs_dn[4]<<"  "<<eigs_dn[5]<<"   "<<eigs_dn[6]<<"  "<<eigs_dn[7]<<endl;


                for(int n=0;n<eigs_up.size();n++){
                    eigs_all.push_back(eigs_up[n]);
                    eigs_all.push_back(eigs_dn[n]);
                }


            }
            eigs_out_fl<<endl;
        }
        eigs_=eigs_all;
        quicksort(eigs_, 0, eigs_.size()-1);

        int n_states_occupied_zeroT = 8.0*Parameters_.ns * Parameters_.Fill * 2.0;
        double initial_mu_guess = 0.5 * (eigs_[n_states_occupied_zeroT - 1] + eigs_[n_states_occupied_zeroT]);
        Parameters_.mus = chemicalpotential(initial_mu_guess, Parameters_.Fill);

        double Energy_QM, Energy_Total;
        Energy_QM=E_QM();

        Energy_Total = Energy_QM + 8.0*Parameters_.lx*Parameters_.ly*(((U_+Parameters_.SIA)*M_*M_)  -  (U_*0.25));

        cout<<setprecision(10);
        cout<<"mu = "<<Parameters_.mus<<endl;
        cout<<"Energy_Total = "<<Energy_Total<<endl;
        cout<<"Energy_Classical = "<<8.0*Parameters_.lx*Parameters_.ly*(((U_+Parameters_.SIA)*M_*M_)  -  (U_*0.25))<<endl;



    }



}

void Hamiltonian::clone(Hamiltonian &Hamiltonian_new){

    Hamiltonian_new.lx_ = lx_;
    Hamiltonian_new.ly_ = ly_;
    Hamiltonian_new.ns_ = ns_;
    Hamiltonian_new.orbs_ = orbs_;
    Hamiltonian_new.HTB_ = HTB_;
    Hamiltonian_new.Ham_ = Ham_;
    Hamiltonian_new.HamCluster_ = HamCluster_;
    Hamiltonian_new.Tx = Tx;
    Hamiltonian_new.Ty = Ty;
    Hamiltonian_new.Tpxpy = Tpxpy;
    Hamiltonian_new.Tpxmy = Tpxmy;
    Hamiltonian_new.eigs_ = eigs_;
    Hamiltonian_new.eigsCluster_saved_ = eigsCluster_;
    Hamiltonian_new.eigs_saved_ = eigs_saved_;
    Hamiltonian_new.sx_ = sx_;
    Hamiltonian_new.sy_ = sy_;
    Hamiltonian_new.sz_ = sz_;
    Hamiltonian_new.HS_factor = HS_factor;

}

double Hamiltonian::chemicalpotential(double muin, double filling)
{
    double mu_out;
    double n1, N;
    double dMubydN;
    double nstate = eigs_.size();
    dMubydN = 0.05 * (eigs_[nstate - 1] - eigs_[0]) / nstate;
    N = filling * double(eigs_.size());
    //temp=Parameters_.temp;
    mu_out = muin;
    bool converged = false;

    if (!Parameters_.fix_mu)
    {
        assert(!Parameters_.fix_mu);

        if (1 == 2)
        {
            for (int i = 0; i < 100000; i++)
            {
                n1 = 0.0;
                for (int j = 0; j < nstate; j++)
                {
                    n1 += double(1.0 / (exp((eigs_[j] - mu_out) * Parameters_.beta) + 1.0));
                }
                //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
                if (abs(N - n1) < double(0.00001))
                {
                    //cout<<abs(N-n1)<<endl;
                    converged = true;
                    break;
                }
                else
                {
                    mu_out += (N - n1) * dMubydN;
                    //cout<<i<<"    "<<n1<<"    "<<N-n1<<endl;
                }
            }

            if (!converged)
            {
                cout << "mu_not_converged, N = " << n1 << endl;
            }
            else
            {
                //cout<<"mu converged, N = "<<n1<<endl;
            }
        }

        double mu1, mu2;
        double mu_temp = muin;
        //cout<<"mu_input = "<<mu_temp<<endl;
        if (1 == 1)
        {
            mu1 = eigs_[0];
            mu2 = eigs_[nstate - 1];
            for (int i = 0; i < 4000000; i++)
            {
                n1 = 0.0;
                for (int j = 0; j < nstate; j++)
                {
                    n1 += double(1.0 / (exp((eigs_[j] - mu_temp) * Parameters_.beta) + 1.0));
                }
                //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
                if (abs(N - n1) < double(0.00001))
                {
                    //cout<<abs(N-n1)<<endl;
                    converged = true;
                    break;
                }
                else
                {
                    if (n1 > N)
                    {
                        mu2 = mu_temp;
                        mu_temp = 0.5 * (mu1 + mu_temp);
                    }
                    else
                    {
                        mu1 = mu_temp;
                        mu_temp = 0.5 * (mu2 + mu_temp);
                    }
                }
                //cout<<"mu_temp = "<<mu_temp<<"   "<<n1<<endl;
            }

            if (!converged)
            {
                cout << "mu_not_converged, N = " << n1 << endl;
            }
            else
            {
                //cout<<"mu converged, N = "<<n1<<endl;
            }

            mu_out = mu_temp;
        }

        return mu_out;
    }
    else
    {
        assert(Parameters_.fix_mu);
        return Parameters_.fixed_mu_value;
    }
} // ----------

double Hamiltonian::chemicalpotentialCluster(double muin, double filling)
{
    double mu_out;
    double n1, N;
    double dMubydN;
    double nstate = eigsCluster_.size();
    dMubydN = 0.05 * (eigsCluster_[nstate - 1] - eigsCluster_[0]) / nstate;
    N = filling * double(eigsCluster_.size());
    //temp=Parameters_.temp;
    mu_out = muin;
    bool converged = false;

    if (!Parameters_.fix_mu)
    {
        assert(!Parameters_.fix_mu);
        if (1 == 2)
        {
            for (int i = 0; i < 100000; i++)
            {
                n1 = 0.0;
                for (int j = 0; j < nstate; j++)
                {
                    n1 += double(1.0 / (exp((eigsCluster_[j] - mu_out) * Parameters_.beta) + 1.0));
                }
                //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
                if (abs(N - n1) < double(0.00001))
                {
                    //cout<<abs(N-n1)<<endl;
                    converged = true;
                    break;
                }
                else
                {
                    mu_out += (N - n1) * dMubydN;
                    //cout<<i<<"    "<<n1<<"    "<<N-n1<<endl;
                }
            }

            if (!converged)
            {
                cout << "mu_not_converged, N = " << n1 << endl;
            }
            else
            {
                //cout<<"mu converged, N = "<<n1<<endl;
            }
        }

        double mu1, mu2;
        double mu_temp = muin;
        //cout<<"mu_input = "<<mu_temp<<endl;
        if (1 == 1)
        {
            mu1 = eigsCluster_[0];
            mu2 = eigsCluster_[nstate - 1];
            for (int i = 0; i < 4000000; i++)
            {
                n1 = 0.0;
                for (int j = 0; j < nstate; j++)
                {
                    n1 += double(1.0 / (exp((eigsCluster_[j] - mu_temp) * Parameters_.beta) + 1.0));
                }
                //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
                if (abs(N - n1) < double(0.00001))
                {
                    //cout<<abs(N-n1)<<endl;
                    converged = true;
                    break;
                }
                else
                {
                    if (n1 > N)
                    {
                        mu2 = mu_temp;
                        mu_temp = 0.5 * (mu1 + mu_temp);
                    }
                    else
                    {
                        mu1 = mu_temp;
                        mu_temp = 0.5 * (mu2 + mu_temp);
                    }
                }
                //cout<<"mu_temp = "<<mu_temp<<"   "<<n1<<endl;
            }

            if (!converged)
            {
                cout << "mu_not_converged, N = " << n1 << endl;
            }
            else
            {
                //cout<<"mu converged, N = "<<n1<<endl;
            }

            mu_out = mu_temp;
        }

        return mu_out;
    }
    else
    {
        assert(Parameters_.fix_mu);
        return Parameters_.fixed_mu_value;
    }
} // ----------

void Hamiltonian::Initialize()
{

    //For Hubbard Stratonovich transformation i.e. MCMF use HS_factor=1.0;
    //HS_factor = 1.0;

    //else use
    HS_factor=0.0;

    int ns = (Parameters_.lx_cluster) * (Parameters_.ly_cluster);

    ly_ = Parameters_.ly;
    lx_ = Parameters_.lx;
    ns_ = Parameters_.ns;

    int space = 2 * ns_;
    int spaceCluster = 2 * ns;

    HTB_.resize(space, space);
    Ham_.resize(space, space);
    Ham_saved_.resize(spaceCluster, spaceCluster);
    HTBCluster_.resize(spaceCluster, spaceCluster);
    HamCluster_.resize(spaceCluster, spaceCluster);
    HamCluster_saved_.resize(spaceCluster, spaceCluster);
    eigs_.resize(space);
    sx_.resize(space);
    sy_.resize(space);
    sz_.resize(space);
    eigs_saved_.resize(space);
    eigsCluster_.resize(spaceCluster);
    eigsCluster_saved_.resize(spaceCluster);

} // ----------

double Hamiltonian::TotalDensity()
{
    double n1 = 0.0;
    for (int j = 0; j < eigs_.size(); j++)
    {
        n1 += 1.0 / (exp(Parameters_.beta * (eigs_[j] - Parameters_.mus * 1.0)) + 1.0);
    }
    return n1;
} // ----------

double Hamiltonian::ClusterDensity()
{
    double n1 = 0.0;
    for (int j = 0; j < eigsCluster_.size(); j++)
    {
        n1 += 1.0 / (exp(Parameters_.beta * (eigsCluster_[j] - Parameters_.mus_Cluster * 1.0)) + 1.0);
    }
    return n1;
} // ----------

double Hamiltonian::E_QM()
{
    double E = 0.0;
    for (int j = 0; j < eigs_.size(); j++)
    {
        //E +=  (eigs_[j]-Parameters_.mus)/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
        E += (eigs_[j]) / (exp(Parameters_.beta * (eigs_[j] - Parameters_.mus)) + 1.0);
    }
    return E;
} // ----------

double Hamiltonian::E_QMCluster()
{
    double E = 0.0;
    for (int j = 0; j < eigsCluster_.size(); j++)
    {
        //E +=  (eigs_[j]-Parameters_.mus)/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
        E += (eigsCluster_[j]) / (exp(Parameters_.beta * (eigsCluster_[j] - Parameters_.mus_Cluster)) + 1.0);
    }
    return E;
} // ----------

double Hamiltonian::GetCLEnergy()
{

    double EClassical;
    int site;
    double ei, ai;
    double phasex, phasey;

    for (int i = 0; i < lx_; i++)
    {
        for (int j = 0; j < ly_; j++)
        {
            site = Coordinates_.Nc(i, j); //+x
            ei = MFParams_.etheta(i, j);
            ai = MFParams_.ephi(i, j);
            sx_[site] = MFParams_.Moment_Size(i, j) * cos(ai) * sin(ei);
            sy_[site] = MFParams_.Moment_Size(i, j) * sin(ai) * sin(ei);
            sz_[site] = MFParams_.Moment_Size(i, j) * cos(ei);
        }
    }

    // Classical Energy
    EClassical = double(0.0);

    if(HS_factor!=0.0){
        for (int ix = 0; ix < lx_; ix++)
        {
            for (int iy = 0; iy < ly_; iy++)
            {

                EClassical += HS_factor * (-0.5) * Parameters_.J_Hund * ((MFParams_.Moment_Size(ix, iy) * MFParams_.Moment_Size(ix, iy)) - (0.25 * MFParams_.Local_density(ix, iy) * MFParams_.Local_density(ix, iy)));

            }
        }
    }

    for (int ix = 0; ix < lx_; ix++)
    {
        for (int iy = 0; iy < ly_; iy++)
        {
            ei = MFParams_.etheta(ix, iy);
            EClassical +=  Parameters_.SIA * ((MFParams_.Moment_Size(ix, iy) * MFParams_.Moment_Size(ix, iy)) *cos(ei)*cos(ei) );

        }
    }



    int _ix, _iy;
    for (int i = 0; i < ns_; i++)
    {
        phasex = 1.0;
        phasey = 1.0;
        _ix = Coordinates_.indx(i);
        _iy = Coordinates_.indy(i);
        if (_ix == (Coordinates_.lx_ - 1))
        {
            phasex = Parameters_.BoundaryConnection;
            phasey = 1.0;
        }
        if (_iy == (Coordinates_.ly_ - 1))
        {
            phasey = Parameters_.BoundaryConnection;
            phasex = 1.0;
        }

        site = Coordinates_.neigh(i, 0); //+x
        EClassical += 1.0 *phasex* Parameters_.K1x * (sx_[i] * sx_[site] + sy_[i] * sy_[site] + 1.0 * sz_[i] * sz_[site]);


        site = Coordinates_.neigh(i, 2); //+y
        EClassical += 1.0 *phasey*Parameters_.K1y * (sx_[i] * sx_[site] + sy_[i] * sy_[site] + 1.0 * sz_[i] * sz_[site]);

        site = Coordinates_.neigh(i, 8); //+2x
        EClassical += 1.0 * phasex*Parameters_.K2x * (sx_[i] * sx_[site] + sy_[i] * sy_[site] + 1.0 * sz_[i] * sz_[site]);
        site = Coordinates_.neigh(i, 9); //+2y
        EClassical += 1.0 *phasey*Parameters_.K2y * (sx_[i] * sx_[site] + sy_[i] * sy_[site] + 1.0 * sz_[i] * sz_[site]);

        site = Coordinates_.neigh(i, 4); //+x+y
        EClassical += 1.0*phasex*phasey*Parameters_.K1_prime * (sx_[i] * sx_[site] + sy_[i] * sy_[site] + 1.0 * sz_[i] * sz_[site]);

        site = Coordinates_.neigh(i, 10); //+2x+2y
        EClassical += 1.0*phasex*phasey*Parameters_.K2_prime * (sx_[i] * sx_[site] + sy_[i] * sy_[site] + 1.0 * sz_[i] * sz_[site]);


        EClassical += (-1.0)*phasex*phasey*Parameters_.t_hopping * ( pow(MFParams_.u_pX(_ix,_iy)  ,2.0) + pow(MFParams_.u_pY(_ix,_iy)  ,2.0) );
    }

    return EClassical;
} // ----------

void Hamiltonian::InteractionsCreate()
{

    /*
            For Effective Hamiltonian derived from
            Hubbard-Stratonovich transformation, remember:
            JH= -2*U
            HS_factor=1.0;
            K=0;
             */

    enum {PX=0,MX,PY,MY};

    int space = 2 * ns_;
    int a;
    double ei, ai;
    double upx, umx, upy, umy;
    int i_mx, i_my;

    double den;
    double fix_mu_double;
    if (Parameters_.fix_mu)
    {
        fix_mu_double = 1.0;
    }
    else
    {
        fix_mu_double = 0.0;
    }

    Ham_ = HTB_;
    // Ham_.print();

    for (int i = 0; i < ns_; i++)
    { // For each site

        i_mx = Coordinates_.neigh(i, MX);
        i_my = Coordinates_.neigh(i, MY);
        upx = MFParams_.u_pX( Coordinates_.indx(i), Coordinates_.indy(i) );
        umx = MFParams_.u_pX( Coordinates_.indx(i_mx), Coordinates_.indy(i_mx) );
        upy = MFParams_.u_pY( Coordinates_.indx(i), Coordinates_.indy(i) );
        umy = MFParams_.u_pY( Coordinates_.indx(i_my), Coordinates_.indy(i_my) );

        ei = MFParams_.etheta(Coordinates_.indx(i), Coordinates_.indy(i));
        ai = MFParams_.ephi(Coordinates_.indx(i), Coordinates_.indy(i));
        den = MFParams_.Local_density(Coordinates_.indx(i), Coordinates_.indy(i));



        Ham_(i, i) += HS_factor * (-0.25) * Parameters_.J_Hund * (den);
        Ham_(i + ns_, i + ns_) += HS_factor * (-0.25) * Parameters_.J_Hund * (den);


        Ham_(i, i) += (-2.0)*Parameters_.SIA * (cos(ei))*  0.5 * MFParams_.Moment_Size(Coordinates_.indx(i), Coordinates_.indy(i));
        Ham_(i + ns_, i + ns_) += (2.0)*Parameters_.SIA * (cos(ei))*  0.5 * MFParams_.Moment_Size(Coordinates_.indx(i), Coordinates_.indy(i));



        Ham_(i, i) += -1.0*Parameters_.t_hopping*Parameters_.lambda_lattice*
                (umx-upx  +  umy-upy);
        Ham_(i + ns_, i + ns_) += -1.0*Parameters_.t_hopping*Parameters_.lambda_lattice*
                (umx-upx  +  umy-upy);

        Ham_(i, i) += Parameters_.J_Hund * (cos(ei)) * 0.5 * MFParams_.Moment_Size(Coordinates_.indx(i), Coordinates_.indy(i));
        Ham_(i + ns_, i + ns_) += Parameters_.J_Hund * (-cos(ei)) * 0.5 * MFParams_.Moment_Size(Coordinates_.indx(i), Coordinates_.indy(i));
        Ham_(i, i + ns_) += Parameters_.J_Hund * sin(ei) * complex<double>(cos(ai), -sin(ai)) * 0.5 * MFParams_.Moment_Size(Coordinates_.indx(i), Coordinates_.indy(i)); //S-
        Ham_(i + ns_, i) += Parameters_.J_Hund * sin(ei) * complex<double>(cos(ai), sin(ai)) * 0.5 * MFParams_.Moment_Size(Coordinates_.indx(i), Coordinates_.indy(i));  //S+

        // On-Site potential
        for (int spin = 0; spin < 2; spin++)
        {
            a = i + ns_ * spin;
            Ham_(a, a) += complex<double>(1.0, 0.0) * MFParams_.Disorder(Coordinates_.indx(i), Coordinates_.indy(i));
        }
    }

} // ----------

void Hamiltonian::InteractionsClusterCreate(int Center_site)
{

    int ns = (Parameters_.lx_cluster) * (Parameters_.ly_cluster);

    enum {PX=0,MX,PY,MY};
    int space = 2 * ns;
    int x_pos, y_pos;
    double ei, ai, den;
    int a;
    int i_original;
    int i_mx, i_my;
    double upx, umx, upy, umy;

    double fix_mu_double;
    if (Parameters_.fix_mu)
    {
        fix_mu_double = 1.0;
    }
    else
    {
        fix_mu_double = 0.0;
    }


    if(Parameters_.ED_){
        InteractionsCreate();
        HamCluster_ = Ham_;
    }
    else{
        HamCluster_ = HTBCluster_;

        for (int i = 0; i < ns; i++)
        { // For each site in cluster
            x_pos = Coordinates_.indx(Center_site) - int(Parameters_.lx_cluster / 2) + CoordinatesCluster_.indx(i);
            y_pos = Coordinates_.indy(Center_site) - int(Parameters_.ly_cluster / 2) + CoordinatesCluster_.indy(i);
            x_pos = (x_pos + Coordinates_.lx_) % Coordinates_.lx_;
            y_pos = (y_pos + Coordinates_.ly_) % Coordinates_.ly_;


            i_original=Coordinates_.Nc(x_pos, y_pos);
            i_mx = Coordinates_.neigh(i_original, MX);
            i_my = Coordinates_.neigh(i_original, MY);
            upx = MFParams_.u_pX( Coordinates_.indx(i_original), Coordinates_.indy(i_original) );
            umx = MFParams_.u_pX( Coordinates_.indx(i_mx), Coordinates_.indy(i_mx) );
            upy = MFParams_.u_pY( Coordinates_.indx(i_original), Coordinates_.indy(i_original) );
            umy = MFParams_.u_pY( Coordinates_.indx(i_my), Coordinates_.indy(i_my) );


            ei = MFParams_.etheta(x_pos, y_pos);
            ai = MFParams_.ephi(x_pos, y_pos);
            den = MFParams_.Local_density(x_pos, y_pos);


            HamCluster_(i, i) += (-2.0)*Parameters_.SIA * (cos(ei))*  0.5 * MFParams_.Moment_Size(x_pos, y_pos);
            HamCluster_(i + ns, i + ns) += (2.0)*Parameters_.SIA * (cos(ei))*  0.5 * MFParams_.Moment_Size(x_pos, y_pos);

            HamCluster_(i, i) += HS_factor * (-0.25) * Parameters_.J_Hund * (den) ;
            HamCluster_(i + ns, i + ns) += HS_factor * (-0.25) * Parameters_.J_Hund * (den) ;

            HamCluster_(i, i) += -1.0*Parameters_.t_hopping*Parameters_.lambda_lattice*(umx-upx  +  umy-upy);
            HamCluster_(i + ns, i + ns) += -1.0*Parameters_.t_hopping*Parameters_.lambda_lattice*(umx-upx  +  umy-upy);

            HamCluster_(i, i) += Parameters_.J_Hund * (cos(ei)) * 0.5 * MFParams_.Moment_Size(x_pos, y_pos);
            HamCluster_(i + ns, i + ns) += Parameters_.J_Hund * (-cos(ei)) * 0.5 * MFParams_.Moment_Size(x_pos, y_pos);
            HamCluster_(i, i + ns) += Parameters_.J_Hund * sin(ei) * complex<double>(cos(ai), -sin(ai)) * 0.5 * MFParams_.Moment_Size(x_pos, y_pos); //S-
            HamCluster_(i + ns, i) += Parameters_.J_Hund * sin(ei) * complex<double>(cos(ai), sin(ai)) * 0.5 * MFParams_.Moment_Size(x_pos, y_pos);  //S+

            for (int spin = 0; spin < 2; spin++)
            {
                a = i + ns * spin;
                HamCluster_(a, a) += complex<double>(1.0, 0.0) * MFParams_.Disorder(x_pos, y_pos);
            }
        }


    }

} // ----------

void Hamiltonian::Check_up_down_symmetry()

{
    complex<double> temp(0, 0);
    complex<double> temp2;

    for (int i = 0; i < ns_; i++)
    {
        for (int j = 0; j < ns_; j++)
        {
            temp2 = Ham_(i, j) - Ham_(i + ns_, j + ns_); //+ Ham_(i+orbs_*ns_,j) + Ham_(i,j+orbs_*ns_);
            temp += temp2 * conj(temp2);
        }
    }

    cout << "Assymetry in up-down sector: " << temp << endl;
}

void Hamiltonian::Check_Hermiticity()

{
    complex<double> temp(0, 0);
    complex<double> temp2;

    for (int i = 0; i < HamCluster_.n_row(); i++)
    {
        for (int j = 0; j < HamCluster_.n_row(); j++)
        {
            if (HamCluster_(i, j) != conj(HamCluster_(j, i)))
            {
                cout << i << "," << j << endl;
                cout << "i,j = " << HamCluster_(i, j) << ", j,i=" << conj(HamCluster_(j, i)) << endl;
            }
            assert(HamCluster_(i, j) == conj(HamCluster_(j, i))); //+ Ham_(i+orbs_*ns_,j) + Ham_(i,j+orbs_*ns_);
            //temp +=temp2*conj(temp2);
        }
    }

    // cout<<"Hermiticity: "<<temp<<endl;
}

void Hamiltonian::Diagonalize(char option)
{

    //extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
    //                       std::complex<double> *,int *, double *, int *);

    char jobz = option;
    // jobz = 'V';
    char uplo = 'U'; //WHY ONLY 'L' WORKS?
    int n = Ham_.n_row();
    int lda = Ham_.n_col();
    vector<complex<double>> work(3);
    vector<double> rwork(3 * n - 2);
    int info;
    int lwork = -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(), eigs_.end(), 0);
    // query:
    zheev_(&jobz, &uplo, &n, &(Ham_(0, 0)), &lda, &(eigs_[0]), &(work[0]), &lwork, &(rwork[0]), &info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0].real()));
    work.resize(lwork);
    // real work:
    zheev_(&jobz, &uplo, &n, &(Ham_(0, 0)), &lda, &(eigs_[0]), &(work[0]), &lwork, &(rwork[0]), &info);
    if (info != 0)
    {
        std::cerr << "info=" << info << "\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}
}

void Hamiltonian::DiagonalizeCluster(char option)
{

    //extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
    //                       std::complex<double> *,int *, double *, int *);

    char jobz = option;
    // jobz = 'V';
    // cout << option;
    char uplo = 'U'; //WHY ONLY 'L' WORKS?
    int n = HamCluster_.n_row();
    int lda = HamCluster_.n_col();
    vector<complex<double>> work(3);
    vector<double> rwork(3 * n - 2);
    int info;
    int lwork = -1;

    eigsCluster_.resize(HamCluster_.n_row());
    fill(eigsCluster_.begin(), eigsCluster_.end(), 0);
    // query:
    zheev_(&jobz, &uplo, &n, &(HamCluster_(0, 0)), &lda, &(eigsCluster_[0]), &(work[0]), &lwork, &(rwork[0]), &info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0].real()));
    work.resize(lwork);
    // real work:
    zheev_(&jobz, &uplo, &n, &(HamCluster_(0, 0)), &lda, &(eigsCluster_[0]), &(work[0]), &lwork, &(rwork[0]), &info);
    if (info != 0)
    {
        std::cerr << "info=" << info << "\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}
}

void Hamiltonian::HTBCreate()
{

    int mx = Parameters_.TBC_mx;
    int my = Parameters_.TBC_my;
    complex<double> phasex, phasey;
    int l, m, a, b;

    HTB_.fill(0.0);

    for (l = 0; l < ns_; l++)
    {

        // * +x direction Neighbor
        if (Coordinates_.indx(l) == (Coordinates_.lx_ - 1))
        {
            phasex = Parameters_.BoundaryConnection*exp(iota_complex * 2.0 * (1.0 * mx) * PI / (1.0 * Parameters_.TBC_cellsX));
            phasey = one_complex;
        }
        else
        {
            phasex = one_complex;
            phasey = one_complex;
        }
        m = Coordinates_.neigh(l, 0);
        for (int spin = 0; spin < 2; spin++)
        {

            a = l + ns_ * spin;
            b = m + ns_ * spin;
            assert(a != b);
            if (a != b)
            {
                HTB_(b, a) = complex<double>(1.0 * Parameters_.t_hopping, 0.0) * phasex;
                HTB_(a, b) = conj(HTB_(b, a));
            }
        }

        // * +y direction Neighbor
        if (Coordinates_.indy(l) == (Coordinates_.ly_ - 1))
        {
            phasex = one_complex;
            phasey = Parameters_.BoundaryConnection*exp(iota_complex * 2.0 * (1.0 * my) * PI / (1.0 * Parameters_.TBC_cellsY));
        }
        else
        {
            phasex = one_complex;
            phasey = one_complex;
        }
        m = Coordinates_.neigh(l, 2);
        for (int spin = 0; spin < 2; spin++)
        {

            a = l + ns_ * spin;
            b = m + ns_ * spin;
            assert(a != b);
            if (a != b)
            {

                HTB_(b, a) = complex<double>(1.0 * Parameters_.t_hopping, 0.0) * phasey;

                HTB_(a, b) = conj(HTB_(b, a));
            }
        }


        if(Parameters_.Geometry=="Triangular"){

            // * +x+y direction Neighbor
            if (Coordinates_.indy(l) == (Coordinates_.ly_ - 1))
            {
                phasey = Parameters_.BoundaryConnection*exp(iota_complex * 2.0 * (1.0 * my) * PI / (1.0 * Parameters_.TBC_cellsY));
            }
            else
            {
                phasey = one_complex;
            }
            if (Coordinates_.indx(l) == (Coordinates_.lx_ - 1))
            {
                phasex = Parameters_.BoundaryConnection*exp(iota_complex * 2.0 * (1.0 * mx) * PI / (1.0 * Parameters_.TBC_cellsX));
            }
            else
            {
                phasex = one_complex;
            }

            m = Coordinates_.neigh(l, 4);
            for (int spin = 0; spin < 2; spin++)
            {

                a = l + ns_ * spin;
                b = m + ns_ * spin;
                assert(a != b);
                if (a != b)
                {

                    HTB_(b, a) = complex<double>(1.0 * Parameters_.t_hopping, 0.0) * phasey*phasex;

                    HTB_(a, b) = conj(HTB_(b, a));
                }
            }

        }


    }

} // ----------

void Hamiltonian::HTBClusterCreate()
{

    if(Parameters_.ED_==false){
        int ns = (Parameters_.lx_cluster) * (Parameters_.ly_cluster);

        complex<double> phasex, phasey;
        int l, m, a, b;

        HTBCluster_.fill(0.0);

        for (l = 0; l < ns; l++)
        {

            // * +x direction Neighbor
            if (CoordinatesCluster_.indx(l) == (CoordinatesCluster_.lx_ - 1))
            {
                phasex = one_complex;
                phasey = one_complex;
            }
            else
            {
                phasex = one_complex;
                phasey = one_complex;
            }
            m = CoordinatesCluster_.neigh(l, 0);
            for (int spin = 0; spin < 2; spin++)
            {

                a = l + ns * spin;
                b = m + ns * spin;
                assert(a != b);
                if (a != b)
                {

                    HTBCluster_(a, b) = complex<double>(1.0 * Parameters_.t_hopping, 0.0) * phasex;
                    HTBCluster_(b, a) = conj(HTBCluster_(a, b));
                }
            }

            // * +y direction Neighbor
            if (CoordinatesCluster_.indy(l) == (CoordinatesCluster_.ly_ - 1))
            {
                phasex = one_complex;
                phasey = one_complex;
            }
            else
            {
                phasex = one_complex;
                phasey = one_complex;
            }
            m = CoordinatesCluster_.neigh(l, 2);
            for (int spin = 0; spin < 2; spin++)
            {

                a = l + ns * spin;
                b = m + ns * spin;
                assert(a != b);
                if (a != b)
                {

                    HTBCluster_(a, b) = complex<double>(1.0 * Parameters_.t_hopping, 0.0) * phasey;
                    HTBCluster_(b, a) = conj(HTBCluster_(a, b));
                }
            }


            if(Parameters_.Geometry=="Triangular"){

                // * +x+y direction Neighbor
                m = CoordinatesCluster_.neigh(l, 4);
                for (int spin = 0; spin < 2; spin++)
                {
                    a = l + ns * spin;
                    b = m + ns * spin;
                    assert(a != b);
                    if (a != b)
                    {
                        HTBCluster_(a, b) = complex<double>(1.0 * Parameters_.t_hopping, 0.0);
                        HTBCluster_(b, a) = conj(HTBCluster_(a, b));
                    }
                }

            }


        }

    }
    else{
        HTBCluster_=HTB_;
    }

    // HTBCluster_.print();

} // ----------

void Hamiltonian::Hoppings()
{
    //DOES SOMETHING EXACT i.e NOTHING :)

} // ----------

void Hamiltonian::copy_eigs(int i)
{

    int space = 2 * ns_;

    if (i == 0)
    {
        for (int j = 0; j < space; j++)
        {
            eigs_[j] = eigs_saved_[j];
        }

        //Ham_ = Ham_saved_;
    }
    else
    {
        for (int j = 0; j < space; j++)
        {
            eigs_saved_[j] = eigs_[j];
        }

        //Ham_saved_ = Ham_;
    }
}

void Hamiltonian::copy_eigs_Cluster(int i)
{

    int ns = (Parameters_.lx_cluster) * (Parameters_.ly_cluster);
    int space = 2 * ns;

    if (i == 0)
    {
        for (int j = 0; j < space; j++)
        {
            eigsCluster_[j] = eigsCluster_saved_[j];
        }

        //HamCluster_ = HamCluster_saved_;
    }
    else
    {
        for (int j = 0; j < space; j++)
        {
            eigsCluster_saved_[j] = eigsCluster_[j];
        }

        //HamCluster_saved_ = HamCluster_;
    }
}

#endif
