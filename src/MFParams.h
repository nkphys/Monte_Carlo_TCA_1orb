#include <math.h>
#include "tensor_type.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "random"
#include <stdlib.h>
#define PI acos(-1.0)

#ifndef MFParams_class
#define MFParams_class

class MFParams
{
public:
    // Define Fields
    Matrix<double> etheta, ephi;
    Matrix<double> u_pX, u_pY;
    Matrix<double> Sz, Sx, Sy;
    Matrix<double> etheta_avg, ephi_avg;
    Matrix<double> Moment_Size;
    Matrix<double> Local_density;
    Matrix<double> Disorder;

    // Constructor
    MFParams(Parameters &Parameters__, Coordinates &Coordinates__, mt19937_64 &Generator1__, mt19937_64 &Generator2__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__), Generator1_(Generator1__), Generator2_(Generator2__)
    {
        //setupm_arr();
        initialize();
    }

    double random1();
    double random2();
    void FieldThrow(int site, string mc_dof_type);
    void initialize();
    void Adjust_MCWindow();
    void Calculate_Fields_Avg();
    void Read_classical_DOFs(string filename);

    Parameters &Parameters_;
    Coordinates &Coordinates_;
    mt19937_64 &Generator1_; //for random fields
    mt19937_64 &Generator2_; //for random disorder
    int lx_, ly_, ns_;

    uniform_real_distribution<double> dis1_; //for random fields
    uniform_real_distribution<double> dis2_; //for random disorder

    //mt19937_64 mt_rand(Parameters_.RandomSeed);
};

void MFParams::Adjust_MCWindow()
{
    double ratio;
    ratio = Parameters_.AccCount[0] / (Parameters_.AccCount[0] + Parameters_.AccCount[1]);
    //cout<<"ratio= "<< ratio << "temp= "<<Parameters_.temp << endl;
    Parameters_.AccCount[0] = 0;
    Parameters_.AccCount[1] = 0;
    Parameters_.WindowSize *= abs(1.0 + 1.0 * (ratio - 0.5));
    if(Parameters_.WindowSize>10){
        Parameters_.WindowSize=10.0;
    }
    //Parameters_.WindowSize =0.2;
    cout << "Ratio: " << ratio << "  window size:  " << Parameters_.WindowSize << endl;
    return;
} // ----------

void MFParams::FieldThrow(int site, string mc_dof_type)
{
    int a, b;

    int Pi_multiple;
    double temp_double;

    double Pi = Parameters_.pi;
    double MC_Window = Parameters_.WindowSize;

    a = Coordinates_.indx(site);
    b = Coordinates_.indy(site);

    //ANGLES
    if (mc_dof_type == "phi")
    {
        ephi(a, b) += 2 * Pi * (random1() - 0.5) * MC_Window;

        Pi_multiple = ephi(a, b)/Pi;


        if (ephi(a, b) < 0.0)
        {
            ephi(a, b) = -ephi(a, b);
        }

        ephi(a, b) = fmod(ephi(a, b), 2.0 * Pi);


    }

    if (mc_dof_type == "theta")
    {
        etheta(a, b) += Pi * (random1() - 0.5) * MC_Window;
        if (etheta(a, b) < 0.0)
        {
            etheta(a, b) = -etheta(a, b);
        }

        etheta(a, b) = fmod(etheta(a, b),  Pi);

    }


    if (mc_dof_type == "theta_and_phi")
    {
        //phi
        ephi(a, b) += 2 * Pi * (random1() - 0.5) * MC_Window;

        Pi_multiple = ephi(a, b)/Pi;

        if (ephi(a, b) < 0.0)
        {
            ephi(a, b) = -ephi(a, b);
        }

        ephi(a, b) = fmod(ephi(a, b), 2.0 * Pi);


        //theta
        etheta(a, b) += Pi * (random1() - 0.5) * MC_Window;
        if (etheta(a, b) < 0.0)
        {
            etheta(a, b) = -etheta(a, b);
        }

        etheta(a, b) = fmod(etheta(a, b),  Pi);
    }


    if (mc_dof_type == "Ising_theta")
    {

        if(random1()>0.5){
            etheta(a, b)=Pi;
        }
        else{
            etheta(a, b)=0;
        }
    }

    if (mc_dof_type == "theta_and_phi_and_u")
    {
        //phi
        ephi(a, b) += 2 * Pi * (random1() - 0.5) * MC_Window;

        Pi_multiple = ephi(a, b)/Pi;

        if (ephi(a, b) < 0.0)
        {
            ephi(a, b) = -ephi(a, b);
        }

        ephi(a, b) = fmod(ephi(a, b), 2.0 * Pi);


        //theta
        etheta(a, b) += Pi * (random1() - 0.5) * MC_Window;
        if (etheta(a, b) < 0.0)
        {
            etheta(a, b) = -etheta(a, b);
        }

        etheta(a, b) = fmod(etheta(a, b),  Pi);


        //u's
        u_pX(a, b) += (random1() - 0.5) * MC_Window;
        u_pY(a, b) += (random1() - 0.5) * MC_Window;
    }


    //Moment Size
    if (mc_dof_type == "moment_size")
    {
        if (Parameters_.MC_on_moment_size == true)
        {
            if(Parameters_.Translational_symmetry_imposed){
                temp_double = abs(Moment_Size(a, b) + ((random1() - 0.5) * MC_Window));
                for(int ix_=0;ix_<lx_;ix_++){
                    for(int iy_=0;iy_<ly_;iy_++){
                        Moment_Size(ix_, iy_) = temp_double;
                    }
                }
            }
            else{
                Moment_Size(a, b) = abs(Moment_Size(a, b) + ((random1() - 0.5) * MC_Window));
            }
        }
    }

    // Local Density
    if (mc_dof_type == "local_density")
    {
        if (Parameters_.MC_on_local_density == true)
        {
            Local_density(a, b) = abs(Local_density(a, b) + ((random1() - 0.5) * MC_Window));
        }
    }

} // ----------

double MFParams::random1()
{

    return dis1_(Generator1_);
}

double MFParams::random2()
{

    return dis2_(Generator2_);
}

void MFParams::initialize()
{

    int spin_offset;
    bool ZigZag_Ising_alongZ=false;
    bool ZigZag_Ising_alongX=false;
    bool Coplanar120_in_XYplane=false;
    bool Coplanar120_in_XZplane=false;
    bool Collinear_alongZ_Piby2Pi=false;
    bool Collinear_alongZ_Piby3Pi=false;
    bool Collinear_alongZ_Piby4Pi=false;
    bool Collinear_alongZ_0Pi=false;

    if(Parameters_.Geometry=="Triangular"){

        //ZigZag_Ising_alongZ=true;
        //ZigZag_Ising_alongX=true;
        //Coplanar120_in_XYplane=true;
        //Coplanar120_in_XZplane=true;
        //Collinear_alongZ_Piby2Pi=true;
        //Collinear_alongZ_Piby3Pi=true;
        //Collinear_alongZ_Piby4Pi=true;
        //Collinear_alongZ_0Pi=true;

        if(ZigZag_Ising_alongZ || ZigZag_Ising_alongX ||
                Coplanar120_in_XYplane || Coplanar120_in_XZplane ||
                Collinear_alongZ_Piby2Pi || Collinear_alongZ_Piby3Pi ||
                Collinear_alongZ_Piby4Pi || Collinear_alongZ_0Pi){
            assert (!Parameters_.MC_on_theta_and_phi_and_u  && !Parameters_.MC_on_theta_and_phi  && !Parameters_.MC_on_theta && !Parameters_.MC_on_phi);
        }

    }

    lx_ = Coordinates_.lx_;
    ly_ = Coordinates_.ly_;

    // srand(Parameters_.RandomSeed);

    etheta_avg.resize(lx_, ly_);
    ephi_avg.resize(lx_, ly_);
    Disorder.resize(lx_, ly_);

    etheta.resize(lx_, ly_);
    ephi.resize(lx_, ly_);
    u_pX.resize(lx_,ly_);
    u_pY.resize(lx_,ly_);
    Moment_Size.resize(lx_, ly_);
    Local_density.resize(lx_, ly_);

    Sz.resize(lx_, ly_);
    Sx.resize(lx_, ly_);
    Sy.resize(lx_, ly_);

    ofstream Disorder_conf_file("Disorder_conf_used");
    Disorder_conf_file << "#seed=" << Parameters_.RandomDisorderSeed << " for mt19937_64 Generator is used" << endl;
    Disorder_conf_file << "#ix   iy    Dis[ix,iy]" << endl;

    ofstream Initial_MC_DOF_file("Initial_MC_DOF_values");
    ofstream Initial_HF_OrderParameters_file("Initial_HF_OrderParameters_values");

    if (Parameters_.Perform_HF_SC_calculation == false)
    {

        Initial_MC_DOF_file << "#seed=" << Parameters_.RandomSeed << " for mt19937_64 Generator is used" << endl;
        Initial_MC_DOF_file << "#ix   iy   Theta(x,y)    Phi(x,y)       Moment_Size(x,y)      Local_density(x,y)    u_pX(x,y)   u_pY(x,y)" << endl;
    }
    else if (Parameters_.Perform_HF_SC_calculation == true)
    {

        Initial_HF_OrderParameters_file << "#seed=" << Parameters_.RandomSeed << " for mt19937_64 Generator is used" << endl;
        Initial_HF_OrderParameters_file << "#ix   iy   Sz(x,y)    Sx(x,y)   Sy(x,y)   Local_density(x,y)   u_pX(x,y)   u_pY(x,y)" << endl;
    }

    // File_Out_Theta_Phi_MicroState << "#x" << setw(15) << "y" << setw(15) << "Theta(x,y)" << setw(15) << "Phi(x,y)"
    //                                                       << setw(15) << "Moment_Size(x,y)" << setw(15) << "Local_density(x,y)" << endl;
    //                         for (int ix = 0; ix < lx_; ix++)
    //                         {
    //                             for (int iy = 0; iy < ly_; iy++)
    //                             {
    //                                 File_Out_Theta_Phi_MicroState << ix << setw(15) << iy << setw(15) << MFParams_.etheta(ix, iy) << setw(15) << MFParams_.ephi(ix, iy)
    //                                                               << setw(15) << MFParams_.Moment_Size(ix, iy) << setw(15) << MFParams_.Local_density(ix, iy) << endl;
    //                             }
    //                         }

    string temp_string;
    int ix_, iy_;
    if (Parameters_.Read_Seed_from_file_ == true)
    {
        ifstream Initial_Seed(Parameters_.Seed_file_name_);
        // for (int i = 0; i < 4; i++)
        // {
        //     Initial_Seed >> temp_string;
        //     cout << temp_string << " ";
        // }
        // cout << endl;
        getline(Initial_Seed, temp_string);
        // cout << temp_string << endl;
        for (int ix = 0; ix < lx_; ix++)
        {
            for (int iy = 0; iy < ly_; iy++)
            {

                Initial_Seed >> ix_ >> iy_ >> etheta(ix, iy) >> ephi(ix, iy) >> Moment_Size(ix, iy) >> Local_density(ix, iy) ;//>> u_pX(ix,iy) >> u_pY(ix,iy);
                //                cout << "ix_=" << ix_ << " ix=" << ix << endl;
                //                cout << "iy_=" << iy_ << " iy=" << iy << endl;


                if(ix_!=ix){
                    cout<<"ix_="<<ix_<<"  ix="<<ix<<endl;
                }
                assert(ix_ == ix);
                assert(iy_ == iy);
                // << ix << setw(15) << iy << setw(15) << MFParams_.etheta(ix, iy) << setw(15) << MFParams_.ephi(ix, iy)
                //   << setw(15) << MFParams_.Moment_Size(ix, iy) << setw(15) << MFParams_.Local_density(ix, iy) << endl;
            }
        }
    }

    else
    {
        for (int j = 0; j < ly_; j++)
        {
            for (int i = 0; i < lx_; i++)
            {
                //ephi(i,j)=(0.5+0.5*pow(-1.0f,i))*Parameters_.pi + grnd()*0.2;
                //etheta(i,j)=0.5*Parameters_.pi + grnd()*0.2;

                //q=(pi,pi)
                // ephi(i,j)=0.0; //(0.5+0.5*pow(-1.0f,i))*Parameters_.pi + grnd()*0.2;
                // etheta(i,j)=0.5*(pow(-1.0,j+i)  + 1.0 )*PI ;//+ grnd()*0.2;

                //q=(0,pi)
                //ephi(i,j)=0.0; //(0.5+0.5*pow(-1.0f,i))*Parameters_.pi + grnd()*0.2;
                //etheta(i,j)=0.5*(pow(-1.0,j)  + 1.0 )*PI; //+ grnd()*0.2;

                //q=(0,0)
                // ephi(i,j)=0.0; //(0.5+0.5*pow(-1.0f,i))*Parameters_.pi + grnd()*0.2;
                // etheta(i,j)=0.0; //+ grnd()*0.2;

                if (Parameters_.Perform_HF_SC_calculation == false)
                {
                    //RANDOM fields

                    if (Parameters_.MC_on_theta_and_phi_and_u == true)
                    {
                        ephi(i, j) = 2.0 * random1() * PI;
                        etheta(i, j) = random1() * PI;
                        u_pX(i,j) = ((2.0 * random1()) - 1.0)*0.5;
                        u_pY(i,j) = ((2.0 * random1()) - 1.0)*0.5;

                    }
                    else if (Parameters_.MC_on_theta_and_phi == true)
                    {
                        ephi(i, j) = 2.0 * random1() * PI;
                        etheta(i, j) = random1() * PI;
                        u_pX(i,j)=0.0;
                        u_pX(i,j)=0.0;
                    }
                    else if (Parameters_.MC_on_Ising_theta == true){
                        ephi(i, j)=0.0;

                        if(random1()>0.5){
                            etheta(i, j)=PI;
                        }
                        else{
                            etheta(i, j)=0;
                        }
                    }
                    else
                    {
                        if (Parameters_.MC_on_phi == true)
                        {
                            ephi(i, j) = 2.0 * random1() * PI;
                        }
                        else
                        {
                            if(ZigZag_Ising_alongX){
                                if(j%2==0){
                                    iy_=j/2;
                                }
                                else{
                                    iy_ = (j-1)/2;
                                }
                                if(iy_%2==0){
                                    spin_offset=1;
                                }
                                else{
                                    spin_offset=-1;
                                }
                                ephi(i, j) = (((pow(-1,i)*spin_offset*1.0) + 1.0) *0.5* PI);

                            }
                            else{
                                ephi(i, j) = 0.0;
                            }

                        }

                        if (Parameters_.MC_on_theta == true)
                        {
                            etheta(i, j) = random1() * PI;
                        }
                        else
                        {
                            if(ZigZag_Ising_alongZ){

                                if(j%2==0){
                                    iy_=j/2;
                                }
                                else{
                                    iy_ = (j-1)/2;
                                }
                                if(iy_%2==0){
                                    spin_offset=1;
                                }
                                else{
                                    spin_offset=-1;
                                }
                                etheta(i, j) = ((pow(-1,i)*spin_offset*1.0) + 1.0) *0.5* PI;


                            }
                            else if(Collinear_alongZ_Piby2Pi){
                                if(i%2==0){
                                    ix_=i/2;
                                }
                                else{
                                    ix_ = (i-1)/2;
                                }
                                if(ix_%2==0){
                                    spin_offset=1;
                                }
                                else{
                                    spin_offset=-1;
                                }

                                etheta(i, j) = ((pow(-1,j)*spin_offset*1.0) + 1.0) *0.5* PI;

                            }
                            else if(Collinear_alongZ_Piby3Pi){
                                if(i%3==0){
                                    ix_=i/3;
                                }
                                else{
                                    ix_ = (i  - (i%3) ) /3;
                                }
                                if(ix_%2==0){
                                    spin_offset=1;
                                }
                                else{
                                    spin_offset=-1;
                                }

                                etheta(i, j) = ((pow(-1,j)*spin_offset*1.0) + 1.0) *0.5* PI;

                            }
                            else if(Collinear_alongZ_Piby4Pi){
                                if(i%4==0){
                                    ix_=i/4;
                                }
                                else{
                                    ix_ = (i  - (i%4) ) /4;
                                }
                                if(ix_%2==0){
                                    spin_offset=1;
                                }
                                else{
                                    spin_offset=-1;
                                }

                                etheta(i, j) = ((pow(-1,j)*spin_offset*1.0) + 1.0) *0.5* PI;
                            }
                            else if(Collinear_alongZ_0Pi){
                                spin_offset=1;
                                etheta(i, j) = ((pow(-1,j)*spin_offset*1.0) + 1.0) *0.5* PI;
                            }

                            else{
                                etheta(i, j) = PI*0.5;
                            }

                        }

                        u_pX(i,j)=0.0;
                        u_pY(i,j)=0.0;
                    }

                    if (Parameters_.MC_on_moment_size == true)
                    {
                        if(!Parameters_.Translational_symmetry_imposed){
                            Moment_Size(i, j) = random1();
                        }
                        else{
                            Moment_Size(i,j)=1.0;
                        }

                        //Moment_Size(i,j)=1.0 - (random1() - 0.5)*0.1;
                    }
                    else
                    {
                        Moment_Size(i, j) = 1.0;
                    }

                    if (Parameters_.MC_on_local_density == true)
                    {
                        Local_density(i, j) = random1();
                    }
                    else
                    {
                        Local_density(i, j) = Parameters_.Fill * 2.0;
                    }

                    if(Coplanar120_in_XYplane){
                        ephi(i, j) = (PI*0.5) + ((i*2.0*PI)/3.0) + ((j*2.0*PI)/3.0);
                        ephi(i, j) = fmod(ephi(i, j), 2.0 * PI);
                    }

                    if(Coplanar120_in_XZplane){

                        ephi(i,j)=0.0;
                        etheta(i, j) = (PI*0.5) + ((i*2.0*PI)/3.0) + ((j*2.0*PI)/3.0);
                        etheta(i, j) = fmod(etheta(i, j), 2.0 * PI);

                        if(etheta(i, j) > PI){
                            etheta(i, j) = (2.0*PI) - etheta(i, j);
                            ephi(i,j)=PI;
                        }

                    }



                }
                else if (Parameters_.Perform_HF_SC_calculation == true)
                {
                    Sz(i, j) = random1();
                    Sx(i, j) = random1();
                    Sy(i, j) = random1();
                    Local_density(i, j) = random1();
                    u_pX(i,j) = ((2.0 * random1()) - 1.0)*0.5;
                    u_pY(i,j) = ((2.0 * random1()) - 1.0)*0.5;
                }

                if (Parameters_.Perform_HF_SC_calculation == false)
                {
                    Initial_MC_DOF_file << i << setw(15) << j << setw(15) << etheta(i, j) << setw(15) << ephi(i, j)
                                        << setw(15) << Moment_Size(i, j) << setw(15) << Local_density(i, j) << setw(15) << u_pX(i,j) << setw(15) << u_pY(i,j) << endl;
                }
                else if (Parameters_.Perform_HF_SC_calculation == true)
                {
                    Initial_HF_OrderParameters_file << i << setw(15) << j << setw(15) << Sz(i, j) << setw(15) << Sx(i, j)
                                                    << setw(15) << Sy(i, j) << setw(15) << Local_density(i, j) << setw(15) << u_pX(i,j) << setw(15) << u_pY(i,j) << endl;
                }
            }

            Initial_HF_OrderParameters_file << endl;
        }
    }

    //RANDOM Disorder
    for (int j = 0; j < ly_; j++)
    {
        for (int i = 0; i < lx_; i++)
        {
            Disorder(i, j) = Parameters_.Disorder_Strength * ((2.0 * random2()) - 1.0);
            Disorder_conf_file << i << "  " << j << "  " << Disorder(i, j) << endl;
        }
        Disorder_conf_file << endl;
    }

} // ----------

void MFParams::Calculate_Fields_Avg()
{

    for (int j = 0; j < ly_; j++)
    {
        for (int i = 0; i < lx_; i++)
        {

            ephi_avg(i, j) = ephi_avg(i, j) + ephi(i, j);
            etheta_avg(i, j) = etheta_avg(i, j) + etheta(i, j);
        }
    }

} // ----------

void MFParams::Read_classical_DOFs(string filename)
{

    string tmp_str;
    double tmp_double;
    ifstream fl_in(filename.c_str());
    fl_in >> tmp_str;

    for (int i = 0; i < lx_; i++)
    {
        for (int j = 0; j < ly_; j++)
        {
            fl_in >> tmp_double >> tmp_double >> etheta(i, j) >> ephi(i, j);
        }
    }

} // ----------

#endif
