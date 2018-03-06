void significance(){
//Get the significance from the FCN values.
//From shanw.
//Here FCN is -lnL (root官网上写的是-logL，但实际的意思就是-lnL).
//reference: 39.3.2 (PDG2016)

  //=***************************
  //With interference
  Double_t logLMax = -15597.8;//括号里面的是L。
  //Without interference
  Double_t logLMin = -15596.6;
  //Number of free parameters with signal - that without signal.
  Int_t ndf=1;
  //=***************************

  Double_t deltaL= logLMax - logLMin;
  Double_t chi2 = fabs(2*deltaL);
  Double_t    prob          = TMath::Prob(chi2,ndf)/2;//Chisq distribution
  Double_t    significance  = RooStats::PValueToSignificance(prob);//p: half
  //5.7*10^(-7)/2 ~ 5 sigma

  printf("Dndf = %d\n",ndf);
  printf("D|logL|=%.2f(%.2f - %.2f)\n",deltaL,logLMax,logLMin);
  printf("Prob =%.2G\n",prob);
  printf("Signif =%.2f\n",significance);
  
}
