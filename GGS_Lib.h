/*
 *  SCLib.cpp
 *  Self-consistent calculations
 *
 *  Created by Yulia Shchadilova on 18/05/2014.
 *  Copyright 2014. All rights reserved.
 *
 */
///============================================


#pragma once
#include "GGS_Param.h"

inline void mix(double & R, double const & newR, double x)
{
    R=(1.-x)*R+x*newR;
}

void mf_step()
{
    E0=P*P/2./M-Xi*Xi/2./M ;
    Xi0=0;
    
    for (int k1=0; k1<NK; k1++)
    for (int th1=0; th1<Ntheta; th1++)
    {
			double k1_=k_vals[k1];                
			double om1=Omega0[k1][th1]+k1_*cos_vals[th1]/M*Xi;
            double dV1=Volume[k1][th1];
            
            Xi0+=k1_*cos_vals[th1]*norm2(alpha[k1][th1])*dV1;
			
			alpha[k1][th1]=-V_vals[k1]/om1;
			
			//energy 
		    E0+=V_vals[k1]*alpha[k1][th1]*dV1;
		
	;}

}


void self_consist_step(double mix_x)
{

       E=P*P/2./M-Xi*(Xi+dXi)/2./M;
 //    E=P*P/2./M-Xi*(Xi)/2./M;
    
    Xi0=0.;
    dXi0=0.;
    
    for (int k1=0; k1<NK; k1++)
    for (int th1=0; th1<Ntheta; th1++)
    {
            double A0xxupd=0., A0xyupd=0., A0yyupd=0.;
            double A1xxupd=0., A1xyupd=0., A1yxupd=0., A1yyupd=0.;
            double A2xxupd=0., A2xyupd=0., A2yyupd=0.;
            
		    double Z0xxupd=0., Z0xyupd=0., Z0yyupd=0.;
		    double Z1xxupd=0., Z1xyupd=0., Z1yxupd=0., Z1yyupd=0.;
		    double Z2xxupd=0., Z2xyupd=0., Z2yyupd=0.;
            
            
			double k1_=k_vals[k1];    
			double ex1=cos_vals[th1], ey1=sin_vals[th1];            
			
			double om1=Omega0[k1][th1]+k1_*ex1/M*Xi;
		    double dV1=Volume[k1][th1];
		    		    
            
            for (int k2=0; k2<NK; k2++)
	        for (int th2=0; th2<Ntheta; th2++)
	        {           
		          
		        double k2_=k_vals[k2];        
			    double ex2=cos_vals[th2], ey2=sin_vals[th2];
				        
		        double om2=Omega0[k2][th2]+k2_*ex2/M*Xi;
		        double dV2=Volume[k2][th2];       
		        
                //---------------------------------------------------------------------  
                //---------------------------------------------------------------------
                
                double A=om1+om2+k1_*k2_/M*ex1*ex2;
                double B=k1_*k2_/M*ey1*ey2;

                //---------------------------------------------------------------------
                // Integral:  In=(2*Pi)^(-1) \int d\phi cos^n(\phi)/(A+B*cos(\phi))
                //---------------------------------------------------------------------
                double I0=1./sqrt((A-B)/(A+B))/(A+B);
                double I1=-1./(A-B)/B*(B+A*(-1.+sqrt((A-B)/(A+B))));
                double I2=A/B/B/(A-B)*(B+A*(-1.+sqrt((A-B)/(A+B))));
                
                //---------------------------------------------------------------------
                // Integral:  IIn=(2*Pi)^(-1) \int d\phi cos^n(\phi)/(A+B*cos(\phi))^2
                //---------------------------------------------------------------------
                double II0=A*pow(A*A-B*B,-3./2);
                double II1=B*pow(A*A-B*B,-3./2);
                double II2=1./B/B*(1.-A*(A*A-2*B*B)*pow(A*A-B*B,-3./2));
                
                
                //---------------------------------------------------------------------  
                //update of A0xx, A0xy, A0yy;
                
                A0xxupd+=dV2*1./M*alpha[k2][th2]*alpha[k2][th2]*k2_*k2_*ex2*ex2*I0;
                A0xyupd+=dV2*1./M*alpha[k2][th2]*alpha[k2][th2]*k2_*k2_*ex2*ey2*I1;
                A0yyupd+=dV2*1./M*alpha[k2][th2]*alpha[k2][th2]*k2_*k2_*ey2*ey2*I2;
                
                A1xxupd+=dV2*1./M*alpha[k2][th2]*alpha[k2][th2]*Fx[k2][th2]*k2_*ex2*I0;
                A1xyupd+=dV2*1./M*alpha[k2][th2]*alpha[k2][th2]*Fx[k2][th2]*k2_*ey2*I1;
                A1yxupd+=dV2*1./M*alpha[k2][th2]*alpha[k2][th2]*Fy[k2][th2]*k2_*ex2*I1;
                A1yyupd+=dV2*1./M*alpha[k2][th2]*alpha[k2][th2]*Fy[k2][th2]*k2_*ey2*I2;
                
                A2xxupd+=dV2*1./M*alpha[k2][th2]*alpha[k2][th2]*
                    Fx[k2][th2]*Fx[k2][th2]*I0;
                A2xyupd+=dV2*1./M*alpha[k2][th2]*alpha[k2][th2]*
                    Fx[k2][th2]*Fy[k2][th2]*I1;
                A2yyupd+=dV2*1./M*alpha[k2][th2]*alpha[k2][th2]*
                    Fy[k2][th2]*Fy[k2][th2]*I2;
                    
                //---------------------------------------------------------------------  
                //update of Zxx, Zxy, Zyy;
                
                Z0xxupd+=dV2*1./M/M*alpha[k2][th2]*alpha[k2][th2]*k2_*k2_*ex2*ex2*II0;
                Z0xyupd+=dV2*1./M/M*alpha[k2][th2]*alpha[k2][th2]*k2_*k2_*ex2*ey2*II1;
                Z0yyupd+=dV2*1./M/M*alpha[k2][th2]*alpha[k2][th2]*k2_*k2_*ey2*ey2*II2;
                
                Z1xxupd+=dV2*1./M/M*alpha[k2][th2]*alpha[k2][th2]*
                                                                Fx[k2][th2]*k2_*ex2*II0;
                Z1xyupd+=dV2*1./M/M*alpha[k2][th2]*alpha[k2][th2]*
                                                                Fx[k2][th2]*k2_*ey2*II1;
                Z1yxupd+=dV2*1./M/M*alpha[k2][th2]*alpha[k2][th2]*
                                                                Fy[k2][th2]*k2_*ex2*II1;
                Z1yyupd+=dV2*1./M/M*alpha[k2][th2]*alpha[k2][th2]*
                                                                Fy[k2][th2]*k2_*ey2*II2;
                
                
                Z2xxupd+=dV2*1./M/M*alpha[k2][th2]*alpha[k2][th2]*
                                                            Fx[k2][th2]*Fx[k2][th2]*II0;
                Z2xyupd+=dV2*1./M/M*alpha[k2][th2]*alpha[k2][th2]*
                                                            Fy[k2][th2]*Fx[k2][th2]*II1;
                Z2yyupd+=dV2*1./M/M*alpha[k2][th2]*alpha[k2][th2]*
                                                            Fy[k2][th2]*Fy[k2][th2]*II2;
		        
		        
	        ;}
            
            
		    //-------------------------------------------------------------------------\\
		    //new Axx, Axy, Ayy;
		    
		    mix( A0xx[k1][th1], A0xxupd , mix_x);
		    mix( A0xy[k1][th1], A0xyupd , mix_x);
		    mix( A0yy[k1][th1], A0yyupd , mix_x);
	
			mix( A1xx[k1][th1], A1xxupd , mix_x);
		    mix( A1xy[k1][th1], A1xyupd , mix_x);
		    mix( A1yx[k1][th1], A1yxupd , mix_x);
		    mix( A1yy[k1][th1], A1yyupd , mix_x);
		    
			mix( A2xx[k1][th1], A2xxupd , mix_x);
		    mix( A2xy[k1][th1], A2xyupd , mix_x);
		    mix( A2yy[k1][th1], A2yyupd , mix_x);	    

	        //-------------------------------------------------------------------------\\
	        //new Zxx, Zxy, Zyy;
		    
			mix( Z0xx[k1][th1], Z0xxupd , mix_x);
		    mix( Z0xy[k1][th1], Z0xyupd , mix_x);
		    mix( Z0yy[k1][th1], Z0yyupd , mix_x);	    	
		    
		    mix( Z1xx[k1][th1], Z1xxupd , mix_x);
		    mix( Z1xy[k1][th1], Z1xyupd , mix_x);
		    mix( Z1yx[k1][th1], Z1yxupd , mix_x);	
		    mix( Z1yy[k1][th1], Z1yyupd , mix_x);	        
		    
		    mix( Z2xx[k1][th1], Z2xxupd , mix_x);
		    mix( Z2xy[k1][th1], Z2xyupd , mix_x);
		    mix( Z2yy[k1][th1], Z2yyupd , mix_x);	    	    
	}
	
	for (int k1=0; k1<NK; k1++)
    for (int th1=0; th1<Ntheta; th1++)
    {
            
			double k1_=k_vals[k1];    
			double ex1=cos_vals[th1], ey1=sin_vals[th1];            
			
			double om1=Omega0[k1][th1]+k1_*ex1/M*Xi;
		    double dV1=Volume[k1][th1];
	
            
            //-------------------------------------------------------------------------\\  
            //new bosons momentum
            
            Xi0+=dV1*k1_*ex1*alpha[k1][th1]*alpha[k1][th1];
            
            dXi0+=dV1*k1_*ex1*alpha[k1][th1]*alpha[k1][th1]*
                (  
                    (k1_*ex1-Fx[k1][th1])*Z0xx[k1][th1]*(k1_*ex1-Fx[k1][th1])
                    +(k1_*ex1-Fx[k1][th1])*2.*Z0xy[k1][th1]*(k1_*ey1-Fy[k1][th1])
                    +(k1_*ey1-Fy[k1][th1])*Z0yy[k1][th1]*(k1_*ey1-Fy[k1][th1])
                    //strange 2
                    -k1_*ex1*2.*Z1xx[k1][th1]*(k1_*ex1-Fx[k1][th1])
                    -k1_*ex1*2.*Z1xy[k1][th1]*(k1_*ey1-Fy[k1][th1])
                    -k1_*ey1*2.*Z1yx[k1][th1]*(k1_*ex1-Fx[k1][th1])
                    -k1_*ey1*2.*Z1yy[k1][th1]*(k1_*ey1-Fy[k1][th1])
                    
                    +k1_*ex1*Z2xx[k1][th1]*k1_*ex1
                    +k1_*ex1*2.*Z2xy[k1][th1]*k1_*ey1
                    +k1_*ey1*Z2yy[k1][th1]*k1_*ey1
                );
                 
                
            
            //-------------------------------------------------------------------------\\  
            //new dispersion relation  
                    
            double om=Omega0[k1][th1]+
            1./M*(
	    	k1_*ex1*(A1xx[k1][th1]-A0xx[k1][th1])*(k1_*ex1-Fx[k1][th1])+
 		    k1_*ex1*(A1xy[k1][th1]-A0xy[k1][th1])*(k1_*ey1-Fy[k1][th1])+
 		    k1_*ey1*(A1yx[k1][th1]-A0xy[k1][th1])*(k1_*ex1-Fx[k1][th1])+
		    k1_*ey1*(A1yy[k1][th1]-A0yy[k1][th1])*(k1_*ey1-Fy[k1][th1])
		    + k1_*k1_*(
		        A1xx[k1][th1]*ex1*ex1+
		        (A1xy[k1][th1]+A1yx[k1][th1])*ex1*ey1+
		        A1yy[k1][th1]*ey1*ey1)
		    - k1_*k1_*(
		        A2xx[k1][th1]*ex1*ex1+
		        2.*A2xy[k1][th1]*ex1*ey1+
		        A2yy[k1][th1]*ey1*ey1)
		    )
		    +k1_*ex1/M*(Xi+dXi);
		    		    
		    mix(Omega[k1][th1], om ,mix_x);
		    
		    
		    //-------------------------------------------------------------------------\\  
		    //new alpha
            double alphaupd=-V_vals[k1]/Omega[k1][th1];
            mix(alpha[k1][th1], alphaupd ,mix_x);
            
            //-------------------------------------------------------------------------\\
		    //new Fx, Fy;
		    
		    double tFx= (k1_*ex1-Fx[k1][th1])*A0xx[k1][th1]
                    +   (k1_*ey1-Fy[k1][th1])*A0xy[k1][th1] 
                    -   k1_*ex1*A1xx[k1][th1]
                    -   k1_*ey1*A1yx[k1][th1];
            double tFy= (k1_*ex1-Fx[k1][th1])*A0xy[k1][th1]
                    +   (k1_*ey1-Fy[k1][th1])*A0yy[k1][th1] 
                    -   k1_*ex1*A1xy[k1][th1]
                    -   k1_*ey1*A1yy[k1][th1];
        
        	mix( Fx[k1][th1], tFx , mix_x);
		    mix( Fy[k1][th1], tFy , mix_x);
            
            //-------------------------------------------------------------------------\\  
            //new Energy
            
            E+=dV1*V_vals[k1]*alpha[k1][th1];

            
            E+=1/2.*dV1*1./M*alpha[k1][th1]*alpha[k1][th1]*
            (
            (k1_*ex1-Fx[k1][th1])*A0xx[k1][th1]*(k1_*ex1-Fx[k1][th1])+
            (k1_*ex1-Fx[k1][th1])*2.*A0xy[k1][th1]*(k1_*ey1-Fy[k1][th1])+
            (k1_*ey1-Fy[k1][th1])*A0yy[k1][th1]*(k1_*ey1-Fy[k1][th1])
            -k1_*ex1*(2.*A1xx[k1][th1])*(k1_*ex1-Fx[k1][th1])
            -k1_*ex1*(A1xy[k1][th1]+A1xy[k1][th1])*(k1_*ey1-Fy[k1][th1])
            -k1_*ey1*(A1yx[k1][th1]+A1yx[k1][th1])*(k1_*ex1-Fx[k1][th1])
            -k1_*ey1*(2.*A1yy[k1][th1])*(k1_*ey1-Fy[k1][th1])
            +k1_*ex1*A2xx[k1][th1]*k1_*ex1+
            +k1_*ex1*2.*A2xy[k1][th1]*k1_*ey1+
            +k1_*ey1*A2yy[k1][th1]*k1_*ey1
            );
		    

		}
		
}

///========\\\========///========\\\========///========\\\========///========\\\=========
///========\\\========///========\\\========///========\\\========///========\\\=========

void save(string Name, double **B)
{

		ofstream ou_file(Name.c_str());
		
		
		for (int th=0; th<Ntheta; th++){
		for (int k=0; k<NK; k++)
		{
			ou_file<<th_vals[th]<<"\t"<<k_vals[k]<<"\t"
			<< B[k][th];
			ou_file<<"\n"<<flush;
		}
		ou_file<<"\n"<<flush;
		}
}
