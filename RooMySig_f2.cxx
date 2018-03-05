/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "RooMySig_f2.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

ClassImp(RooMySig_f2) 

 RooMySig_f2::RooMySig_f2(const char *name, const char *title, 
                        RooAbsReal& _m,
                        RooAbsReal& _mean,
                        RooAbsReal& _Gamma,
                        RooAbsReal& _phi,
                        RooAbsReal& _ratio,
                        RooAbsReal& _p1,
                        RooAbsReal& _p2,
                        RooAbsReal& _p3,
                        RooAbsReal& _p4,
                        RooAbsReal& _p5,
                        RooAbsReal& _p6) :
   RooAbsPdf(name,title), 
   m("m","m",this,_m),
   mean("mean","mean",this,_mean),
   Gamma("Gamma","Gamma",this,_Gamma),
   phi("phi","phi",this,_phi),
   ratio("ratio","ratio",this,_ratio),
   p1("p1","p1",this,_p1),
   p2("p2","p2",this,_p2),
   p3("p3","p3",this,_p3),
   p4("p4","p4",this,_p4),
   p5("p5","p5",this,_p5),
   p6("p6","p6",this,_p6)
 { 
 } 


 RooMySig_f2::RooMySig_f2(const RooMySig_f2& other, const char* name) :  
   RooAbsPdf(other,name), 
   m("m",this,other.m),
   mean("mean",this,other.mean),
   Gamma("Gamma",this,other.Gamma),
   phi("phi",this,other.phi),
   ratio("ratio",this,other.ratio),
   p1("p1",this,other.p1),
   p2("p2",this,other.p2),
   p3("p3",this,other.p3),
   p4("p4",this,other.p4),
   p5("p5",this,other.p5),
   p6("p6",this,other.p6)
 { 
 } 



 Double_t RooMySig_f2::evaluate() const 
 { 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE

     //The range is for 3-body bkg
     Double_t low=1.96;//
     Double_t up=2.56;//

     //sig: |BW|^2
     Double_t sig=(1./( (m*m-mean*mean)*(m*m-mean*mean)+(mean*Gamma)*(mean*Gamma) ));

     //efficiency
     Double_t eff=0;
     if (m>=1.977&&m<=2.549) eff=2.978-4.139*m+1.883*m*m-0.2766*m*m*m;

     //bkg: poly
     Double_t bkg=0.;
     if ( (m>=low)&&(m<=up) ) {
       Double_t X=-1.0+2.0*( (m-low)/(up-low) );

       Double_t T0=1.0;
       Double_t T1=X;
       Double_t T2=2.0*X*X-1.0;
       Double_t T3=4.0*X*X*X-3.0*X;
       Double_t T4=8.0*X*X*X*X-8.0*X*X+1.0;
       Double_t T5=16.0*X*X*X*X*X-20.0*X*X*X+5.0*X;
       Double_t T6=32.0*X*X*X*X*X*X-48.0*X*X*X*X+18.0*X*X-1.0;

       bkg=(T0+p1*T1+p2*T2+p3*T3+p4*T4+p5*T5+p6*T6)*ratio*ratio;
       if (bkg<0.) bkg=0.;
     }//end with if ( (m>low)&&(m<up) )

     //inter
     Double_t inter=2.*( cos(phi)*(m*m-mean*mean)-sin(phi)*mean*Gamma ) / ( (m*m-mean*mean)*(m*m-mean*mean)+(mean*Gamma)*(mean*Gamma) ) *sqrt(bkg);

     return (sig+inter+bkg)*eff; 
 } 



