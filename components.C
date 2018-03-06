#include <math.h>
#include "TMath.h"

double components(double *x, double *par){
// prepare for drawing

   double m=x[0];
   double mean=par[0];
   double Gamma=par[1];
   double phi=par[2];
   double ratio=par[3];
   double p1 = par[4];
   double p2 = par[5];
   double p3 = par[6];
   double p4 = par[7];
   double p5 = par[8];
   double p6 = par[9];

   double f_sig=par[10];
   double f_bkg=par[11];
   double f_intf=par[12];
   double scale=par[13];

    //PHSP factors
     Double_t m_jpsi=3.097;
     Double_t m_phi=1.019;
     Double_t m_eta=0.548;
     Double_t m_etap=0.958;

   Double_t p=0;
   if ( (m*m-(m_phi+m_etap)*(m_phi+m_etap))>0. ) {
     p=sqrt( (m*m-(m_phi+m_etap)*(m_phi+m_etap))*(m*m-(m_phi-m_etap)*(m_phi-m_etap)) ) / (2*m);
   }
   if (p<0) p=0;

   Double_t q=0;
   if ( (m_jpsi*m_jpsi-(m_eta+m)*(m_eta+m))>0. ) {
     q=sqrt( (m_jpsi*m_jpsi-(m_eta+m)*(m_eta+m))*(m_jpsi*m_jpsi-(m_eta-m)*(m_eta-m)) ) / (2*m_jpsi);
   }
   if (q<0) q=0;

    //PHSP factor=p^{l1+1/2}*q^{l2+1/2}
    Double_t PHSPf=p*sqrt(p)*q*sqrt(q);

   //sig: |BW|^2
   Double_t sig=1./( (m*m-mean*mean)*(m*m-mean*mean)+(mean*Gamma)*(mean*Gamma) )*PHSPf*PHSPf;

   //bkg: poly
   Double_t low=1.96;//
   Double_t up=2.56;//
   Double_t X=-1.0+2.0*( (m-low)/(up-low) );

   Double_t T0=1.0;
   Double_t T1=X;
   Double_t T2=2.0*X*X-1.0;
   Double_t T3=4.0*X*X*X-3.0*X;
   Double_t T4=8.0*X*X*X*X-8.0*X*X+1.0;
   Double_t T5=16.0*X*X*X*X*X-20.0*X*X*X+5.0*X;
   Double_t T6=32.0*X*X*X*X*X*X-48.0*X*X*X*X+18.0*X*X-1.0;

   Double_t bkg=(T0+p1*T1+p2*T2+p3*T3+p4*T4+p5*T5+p6*T6)*ratio*ratio;
   if (bkg<0) bkg=0.0;

   //inter
   Double_t inter=2.0*( cos(phi)*(m*m-mean*mean)-sin(phi)*mean*Gamma ) / ( (m*m-mean*mean)*(m*m-mean*mean)+(mean*Gamma)*(mean*Gamma) ) *sqrt(bkg) * PHSPf;

   //efficiency
   Double_t eff=0;
   if (m>=1.977&&m<=2.549) eff=2.372-3.572*m+1.741*m*m-0.2684*m*m*m;

   return (f_sig*sig + f_bkg*bkg + f_intf*inter)*scale*eff;
} 
