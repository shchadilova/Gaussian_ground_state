/*
 *  SCParam.h
 *  Self-consistent calculations
 *
 *  Created by Yulia Shchadilova on 18/05/2014.
 *  Copyright 2014. All rights reserved.
 *
 */
///============================================
#pragma once

#include <iomanip>
#include <math.h>
#include <vector>
#include <complex>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>



#define MSG_PREFIX            __FILE__ << ":" << __LINE__ << ": "
#define DEBUG(MSG)            std::cout << MSG_PREFIX << MSG << std::endl;



    typedef std::complex<double> complex_type;
    typedef std::vector<complex_type> vector_type;
    typedef double real_type;
    typedef std::vector<real_type> real_vector_type;

#include "GGS_headers.cpp"

///===========///===========///===========///===========///===========///===========///
double Kmax=200.; //UV cut-off
double chi=1e-5; //IR cut-off

int NK=100;
int Ntheta=150;

int Iter=100000;
double Tol=1e-4;
double xTol=1e-5;
///===========
double xi= 1.;
double mB = 1./sqrt(2.);

double P=0.;
double c=1.;
double M =0.263158*mB; //
double sigma=0.;

double a=0.1;
double aIB =sqrt(a/8./Pi);
double nBEC = 1.;
double gIB = 2.*Pi*(1./mB+1./M)*aIB;
double mu=1./(1./mB+1./M);
///===========
double E=0.;
double E0=0.;
double Xi=0.;
double Xi0=0.;
double dXi=0.;
double dXi0=0.;
double x=0.1;
///===========///===========///===========///===========///===========///===========///


ofstream debug("debug.dat");

real_vector_type k_vals;
real_vector_type om_vals;
real_vector_type dk_vec;
real_vector_type th_vals;
real_vector_type dth_vec;
real_vector_type cos_vals;
real_vector_type sin_vals;
real_vector_type V_vals;

double  ** alpha;
double  ** Fx, ** Fy;
double  ** A0xx, ** A0xy, ** A0yy;
double  ** A1xx, ** A1xy, ** A1yx, ** A1yy;
double  ** A2xx, ** A2xy, ** A2yy;
double  ** Z0xx, ** Z0yy, ** Z0xy;
double  ** Z1xx, ** Z1yy, ** Z1xy, ** Z1yx;
double  ** Z2xx, ** Z2yy, ** Z2xy;
double  ** Volume, ** Omega0, **Omega;

///============================================================

/// grids
double grid_ed(double k, int NK, double Kmax)
{
    return double(k+1.)/NK*Kmax;
}

double grid_pl(double k, int NK, double Kmax)
{	
    return chi+pow(k,2.)/pow(NK-1,2.)*(Kmax-chi);
}

double grid_pl1(double k, int NK, double Kmax)
{	
    return pow(k+1.,2.)/pow(NK+1.,2.)*(Kmax);
}

double grid_exp(double k, int NK, double Kmax)
{	
    return chi+exp(k)/exp(NK-1)*(Kmax-chi);
}


///============================================================

void Ini()
{
// alpha- displacement

	
	k_vals.resize(NK);
	th_vals.resize(Ntheta);
	cos_vals.resize(Ntheta);
	sin_vals.resize(Ntheta);
	
	V_vals.resize(NK);
	om_vals.resize(NK);

	
	for (int i=0; i<NK; i++)
	{
        k_vals[i]=grid_pl(i,NK,Kmax);
        V_vals[i]=sqrt(nBEC)*gIB*sqrt(sqrt(k_vals[i]*k_vals[i]/(2.+k_vals[i]*k_vals[i])))
        /sqrt(2.*Pi)/(2.*Pi)*exp(-k_vals[i]*k_vals[i]*sigma*sigma/2.);
        om_vals[i]=c*k_vals[i]*sqrt(1.+xi*xi*k_vals[i]*k_vals[i]/2.);
        //      cout <<  k_vals[i] <<endl;
	}
	
	for (int th=0; th<Ntheta; th++)
	{
	    th_vals[th]=double(th+1.)/(Ntheta+1.)*Pi;
	    cos_vals[th]=cos(th_vals[th]);
	    sin_vals[th]=sin(th_vals[th]);
	    //	    cout <<  th_vals[th] << "\t"<<Pi-th_vals[th]<<endl;
	}
	
	dk_vec.resize(k_vals.size());
    for(int k=0; k<NK; k++)
    {   
        dk_vec[k]=k_vals[(k+1)% NK]-k_vals[k];
        //if(dk_vec[k]<0) dk_vec[k]+=Kmax;
        if(dk_vec[k]<0) dk_vec[k]=0;
       //cout <<  dk_vec[k] <<endl;
    };
    
    dth_vec.resize(th_vals.size());
    for(int th=0; th<Ntheta; th++)
    {   
        dth_vec[th]=th_vals[(th+1)% Ntheta]-th_vals[th];
        if(dth_vec[th]<0) dth_vec[th]=dth_vec[th-1];
        //cout <<  dth_vec[th] <<endl;
    };
	
	{
	
	alpha=new double *[NK];
	
	Fx=new double *[NK];    Fy=new double *[NK];
	
	A0xx=new double *[NK];	A0xy=new double *[NK];	A0yy=new double *[NK];
	A1xx=new double *[NK];	A1xy=new double *[NK];	
	A1yx=new double *[NK];A1yy=new double *[NK];
	A2xx=new double *[NK];	A2xy=new double *[NK];	A2yy=new double *[NK];	
	
	Z0xx=new double *[NK];	Z0yy=new double *[NK];	Z0xy=new double *[NK];	
	Z1xx=new double *[NK];	Z1yy=new double *[NK];	
	Z1xy=new double *[NK];	Z1yx=new double *[NK];
	Z2xx=new double *[NK];	Z2yy=new double *[NK];	Z2xy=new double *[NK];	
	
	Volume=new double *[NK];
	Omega0=new double *[NK];
	Omega=new double *[NK];
	
	for (int k=0; k<NK; k++)
	{
		alpha[k]=new double [Ntheta];
		
		Fx[k]=new double [Ntheta];  
		Fy[k]=new double [Ntheta];
		
		A0xx[k]=new double [Ntheta]; 
		A0xy[k]=new double [Ntheta];
		A0yy[k]=new double [Ntheta];
        
		A1xx[k]=new double [Ntheta]; 
		A1xy[k]=new double [Ntheta];
		A1yx[k]=new double [Ntheta];
		A1yy[k]=new double [Ntheta];
		
		A2xx[k]=new double [Ntheta]; 
		A2xy[k]=new double [Ntheta];
		A2yy[k]=new double [Ntheta];		
        
        Z0xx[k]=new double [Ntheta];
        Z0yy[k]=new double [Ntheta];
        Z0xy[k]=new double [Ntheta];
        
        Z1xx[k]=new double [Ntheta];
        Z1yy[k]=new double [Ntheta];
        Z1xy[k]=new double [Ntheta];
        Z1yx[k]=new double [Ntheta];

        Z2xx[k]=new double [Ntheta];
        Z2yy[k]=new double [Ntheta];
        Z2xy[k]=new double [Ntheta];
        								
		Volume[k]=new double [Ntheta];
		Omega0[k]=new double [Ntheta];
		Omega[k]=new double [Ntheta];
				
		for (int th=0; th<Ntheta; th++)
		{ 
		    double k_=k_vals[k];
		
		    alpha[k][th]=0.; 
		
		    A0xx[k][th]=0.;
		    A0xy[k][th]=0.;
		    A0yy[k][th]=0.;

		    A1xx[k][th]=0.;
		    A1xy[k][th]=0.;
		    A1yx[k][th]=0.;
		    A1yy[k][th]=0.;

		    A2xx[k][th]=0.;
		    A2xy[k][th]=0.;
		    A2yy[k][th]=0.;
		    		
            Fx[k][th]= 0.;
            Fy[k][th]= 0.;
            
            Z0xx[k][th]= 0.;
            Z0yy[k][th]= 0.;
            Z0xy[k][th]= 0.;
            
            Z1xx[k][th]= 0.;
            Z1yy[k][th]= 0.;
            Z1xy[k][th]= 0.;
            Z1yx[k][th]= 0.; 

            
            Z2xx[k][th]= 0.;
            Z2yy[k][th]= 0.;
            Z2xy[k][th]= 0.;
                    
		    Omega0[k][th]=(om_vals[k]+k_*k_/(2.*M)-k_*cos_vals[th]/M*P);
		    Omega[k][th]=0.;
		
		    Volume[k][th]=2.*Pi*k_vals[k]*k_vals[k]*sin_vals[th]*dth_vec[th]*
		    (dk_vec[k]+dk_vec[k-1])/2.;
		
		}
	;}
	;}
	
	
;}


///============================================================

void IniP()
{

	Omega0=new double *[NK];
	for (int k=0; k<NK; k++)
	{
	    Omega0[k]=new double [Ntheta];
	    for (int th=0; th<Ntheta; th++){ 
		    double k_=k_vals[k];
		    Omega0[k][th]=(om_vals[k]+k_*k_/(2.*M)-k_*cos_vals[th]/M*P);
		}
	;}
;}

///============================================================

void IniA()
{
    	V_vals.resize(NK);
	
	for (int i=0; i<NK; i++)
	{
        V_vals[i]=sqrt(nBEC)*gIB*sqrt(sqrt(k_vals[i]*k_vals[i]/(2.+k_vals[i]*k_vals[i])))
        /sqrt(2.*Pi)/(2.*Pi)*exp(-k_vals[i]*k_vals[i]*sigma*sigma/2.);
	}
;}

///============================================================
///============================================================
