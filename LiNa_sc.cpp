/*
 *  
 *  Self-consistent calculations
 *
 *  Created by Yulia Shchadilova on 18/05/2014.
 *  Copyright 2014. All rights reserved.
 *
 */
///============================================
#include "LiNa_Lib.h"

int main(int argc, char **argv)
{
	
	
    int Adots=10;
    int Pdots=1;
    int NKdots=1;
    
    
	//part that creates file with the proper name
	std::ostringstream os;
    os << "LiNa_2000_100.dat";
    std::string Filename = os.str();
    ofstream data(Filename.c_str());
    



	for(int ii=1; ii<=NKdots;ii++){
    //grid
    
    Kmax=2000.;
    NK=100*ii;   
    Ntheta=100*ii;  
    chi=1e-5;
    //Iter=2;
    
	Ini();
	for(int na=1; na<=Adots;na++){	

	a=0.1*na;
	aIB = sqrt(a/8./Pi/nBEC);
	gIB = 2.*Pi*(1./mB+1./M)*aIB;
	IniA();
    
    for(int np=1; np<=Pdots;np++){	
    P=0.*M*np;  
    //IniP();
    
	for(int i=0;i<Iter;i++)
	{
	if(i<1)  self_consist_step(1.);  
	else { 
    double E0=E;
    self_consist_step(x);  
    
    
//     if(abs(Xi0-Xi)<Tol && i>1){i=Iter;}
//     else{    Xi=Xi0;
//              dXi=dXi0;};
    Xi=Xi0;
    dXi=dXi0;
    
    if(abs(E0-E)<Tol && i>1){i=Iter;};
    
    double EMinf0=4.*aIB*aIB*nBEC*(1./mB+1./M)*k_vals[NK-1]*mB;
	cout << a << "\t" << P<<"\t"<< 1./NK  << "\t" << EMinf0+E*mB<< "\t" << Xi << "\t" << dXi <<endl;  
	}
	}
	
    double EMinf0=4.*aIB*aIB*nBEC*(1./mB+1./M)*k_vals[NK-1]*mB;
    data << a << "\t" << P<<"\t"<< 1./NK <<"\t" <<EMinf0+E*mB<< "\t" << Xi << "\t" << dXi <<endl;


	
	
	};};}
	
	
// 	save("alpha.dat",alpha);
//  	save("omega.dat",Omega);
// 
// 	save("Fx.dat",Fx);
// 	save("Fy.dat",Fy);
// 	
// 	save("a0xx.dat",A0xx);
// 	save("a1xx.dat",A1xx);
// 	save("Z0xx.dat",Z0yy);
// 	save("Z1xx.dat",Z1yy);
// 	
	return 0;
;}