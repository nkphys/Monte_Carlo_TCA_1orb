#include "tensor_type.h"
#include "ParametersEngine.h"


#ifndef Coordinates_class
#define Coordinates_class
class Coordinates { 

public:
    Coordinates(int lx, int ly)
        :lx_(lx),ly_(ly)
    {
        Numbering();
    }

    enum {PX=0,MX,PY,MY,PXPY,MXPY,MXMY,PXMY};	// not needed - but keep as a "key"
    void Numbering();
    int indx(int i);
    int indy(int i);
    int Nc(int x, int y);
    int NNc(int x, int y);
    int neigh(int site, int wneigh);
    int getneigh(int site,int wneigh);

    int lx_,ly_,ns_;
    Mat_1_int indx_,indy_;
    Matrix<int> Nc_,neigh_;
};

/*  
 * ***********
 *  Functions in Class Coordinates ------
 *  ***********
*/  

int Coordinates::indx(int i){
    if(i>ns_-1){perror("Coordinates.h:x-coordinate of lattice excede limit");}
    return indx_[i];
 // ----------
}


int Coordinates::indy(int i){
    if(i>ns_-1){perror("Coordinates.h:y-coordinate of lattice excede limit");}
    return indy_[i];
 // ----------
}

int Coordinates::Nc(int x, int y){
    if(x>lx_&&y<ly_){perror("Coordinates.h:ith-sitelabel of lattice excede limit");}
    return Nc_(x,y);
 // ----------
}


int Coordinates::neigh(int site, int wneigh){
    if(site>ns_-1 || wneigh>=11){perror("Coordinates.h:getneigh -> ifstatement-1");}
    return neigh_(site,wneigh);
} // ----------


void Coordinates::Numbering(){

    ns_=lx_*ly_;

    indx_.clear(); 	indx_.resize(ns_);
    indy_.clear();	indy_.resize(ns_);
    Nc_.resize(lx_,ly_);
    neigh_.resize(ns_,11);

    // Site labeling
    int icount=0;
    for(int j=0;j<ly_;j++){
        for(int i=0;i<lx_;i++){
            indx_[icount]=i;
            indy_[icount]=j;
            Nc_(i,j)=icount;
            icount++;
        }}

    // Neighbors for each site
    for(int i=0;i<ns_;i++){ 	// ith site
        for(int j=0;j<11;j++) {		// jth neighbor
            neigh_(i,j)=getneigh(i,j);
        }
    }

} // ----------


int Coordinates::getneigh(int site,int wneigh){
    if(site>ns_-1 || wneigh>10){perror("Coordinates.h:getneigh -> ifstatement-1");}
    int nx=indx(site);
    int ny=indy(site);
    int mx=0;
    int my=0;

    // Nearest Neighbours
    if(wneigh==0){ //PX
        mx=(nx+1)%(lx_);
        my=ny;
    }
    if(wneigh==1){ //MX
        mx=(nx+lx_-1)%(lx_);
        my=ny;
    }
    if(wneigh==2){ //PY
        mx=nx;
        my=(ny+1)%(ly_);
    }
    if(wneigh==3){ //MY
        mx=nx;
        my=(ny+ly_-1)%(ly_);
    }


    // Next-Nearest for sqrt(2) distance!
    if(wneigh==4){ //PXPY
        mx=(nx+1)%(lx_);
        my=(ny+1)%(ly_);
    }
    if(wneigh==5){ //MXPY
        mx=(nx+lx_-1)%(lx_);
        my=(ny+1)%(ly_);
    }
    if(wneigh==6){ //MXMY
        mx=(nx+lx_-1)%(lx_);
        my=(ny+ly_-1)%(ly_);
    }
    if(wneigh==7){ //PXMY
        mx=(nx+1)%(lx_);
        my=(ny+ly_-1)%(ly_);
    }


    //Next-Nearest for dis=2
    if(wneigh==8){ //2*PX
        mx=(nx+2)%(lx_);
        my=ny;
    }
    if(wneigh==9){ //2*PY
        mx=nx;
        my=(ny+2)%(ly_);
    }

    if(wneigh==10){ //2*PY + 2*PX
        mx=(nx+2)%(lx_);
        my=(ny+2)%(ly_);
    }



    return Nc_(mx,my); //Nc(mx,my);
} // ----------

#endif
