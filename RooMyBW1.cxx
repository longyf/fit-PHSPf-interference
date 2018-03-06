/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "RooMyBW1.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

ClassImp(RooMyBW1) 

 RooMyBW1::RooMyBW1(const char *name, const char *title, 
                        RooAbsReal& _m_phi_etap,
                        RooAbsReal& _mean,
                        RooAbsReal& _Gamma) :
   RooAbsPdf(name,title), 
   m_phi_etap("m_phi_etap","m_phi_etap",this,_m_phi_etap),
   mean("mean","mean",this,_mean),
   Gamma("Gamma","Gamma",this,_Gamma)
 { 
 } 


 RooMyBW1::RooMyBW1(const RooMyBW1& other, const char* name) :  
   RooAbsPdf(other,name), 
   m_phi_etap("m_phi_etap",this,other.m_phi_etap),
   mean("mean",this,other.mean),
   Gamma("Gamma",this,other.Gamma)
 { 
 } 



 Double_t RooMyBW1::evaluate() const 
 { 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE

   Double_t m_jpsi=3.097;
   Double_t m_phi=1.019;
   Double_t m_eta=0.548;
   Double_t m_etap=0.958;

   Double_t BWsq=1./((m_phi_etap*m_phi_etap-mean*mean)*(m_phi_etap*m_phi_etap-mean*mean)+mean*Gamma*mean*Gamma);

   Double_t p=0;
   if ( (m_phi_etap*m_phi_etap-(m_phi+m_etap)*(m_phi+m_etap))>0. ) {
     p=sqrt( (m_phi_etap*m_phi_etap-(m_phi+m_etap)*(m_phi+m_etap))*(m_phi_etap*m_phi_etap-(m_phi-m_etap)*(m_phi-m_etap)) ) / (2*m_phi_etap);
   }
   if (p<0) p=0; 

   Double_t q=0;
   if ( (m_jpsi*m_jpsi-(m_eta+m_phi_etap)*(m_eta+m_phi_etap))>0. ) {
     q=sqrt( (m_jpsi*m_jpsi-(m_eta+m_phi_etap)*(m_eta+m_phi_etap))*(m_jpsi*m_jpsi-(m_eta-m_phi_etap)*(m_eta-m_phi_etap)) ) / (2*m_jpsi);
   }
   if (q<0) q=0;

   //efficiency
   Double_t eff=0;
   if (m_phi_etap>=1.977&&m_phi_etap<=2.549) eff=2.372-3.572*m_phi_etap+1.741*m_phi_etap*m_phi_etap-0.2684*m_phi_etap*m_phi_etap*m_phi_etap;

   //|BW|^2*eff*(p*q)^3
   return BWsq*eff*p*q*p*q*p*q;

   //|BW|^2*eff
//   return BWsq*eff;
 } 



