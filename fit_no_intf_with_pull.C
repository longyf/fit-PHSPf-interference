
using namespace RooFit;

// This fit is much better than the nominal result.
// The difference between this fit and the nominal one is:
// 1. 3-body bkg
//    this fit: 6th-order chebychev polynominal
//    nominal: shape from RooKeysPdf
// 2. ratio (dominant in improving the fit)
//    this fit: not fixed.
//    nominal: nsig1/nsig2 is fixed, nbkg1_1/nbkg1_2 is also fixed. The ratio of BR is the same for signal and PHSP processes.

void fit_no_intf_with_pull(){
// fit to M(phi etap) without considering interference between resonance and 3-body bkg
// 3-body is described by using 6th polynomial

  gSystem->Load("libRooFit");
  //signal (nominal)
  gROOT->ProcessLine(".L ./RooMyBW1.cxx++");
  gROOT->ProcessLine(".L ./RooMyBW2.cxx++");
  //3-body bkg
  gROOT->ProcessLine(".L /Users/long/workarea/Phietaetap/interference/bkg/RooMyBkg1.cxx++");
  gROOT->ProcessLine(".L /Users/long/workarea/Phietaetap/interference/bkg/RooMyBkg2.cxx++");

  double low_=1.96;//
  double up_=2.56;//

  RooRealVar m_phi_etap("m_phi_etap","m_phi_etap",low_,up_);//

  //Read data.
  //////////////////////
  //From final state 1:
  TFile *f_data_f1=new TFile ("/Users/long/workarea/G3Pi2K2/output_data_g3pi2k2.root");
  TTree *tree_f1 =(TTree*) f_data_f1->Get("TreeRes");
  RooDataSet data_f1("data_f1","data_f1",tree_f1,RooArgSet(m_phi_etap));

  //From final state 2:
  TFile *f_data_f2=new TFile ("/Users/long/workarea/G4Pi2K2/output_data_g4pi2k2.root");
  TTree *tree_f2 =(TTree*) f_data_f2->Get("TreeRes");
  RooDataSet data_f2("data_f2","data_f2",tree_f2,RooArgSet(m_phi_etap));

  //Background( f0(1500) )
  //Read the shape of Jpsi->phi f0(1500).
  TFile *f_bkg2_f1=new TFile("/Users/long/workarea/G3Pi2K2/fit_final/Pdf_g3pi2k2_6.root");
  RooKeysPdf *bkg2_f1=(RooKeysPdf*)f_bkg2_f1 -> Get("bkg2");

  //Read the shape of Jpsi->phi f0(1500)
  TFile *f_bkg2_f2=new TFile("/Users/long/workarea/G4Pi2K2/fit_final/Pdf_g4pi2k2_6.root");
  RooKeysPdf *bkg2_f2=(RooKeysPdf*)f_bkg2_f2 -> Get("bkg2");

  //Bckground sideband
  /////////////////////
  //From final state 1:
  //Sidebands of phi and etap: chebychev.
  RooRealVar c1_f1("c1_f1","",  0.15) ;
  RooRealVar c2_f1("c2_f1","", -0.548) ;
  RooRealVar c3_f1("c3_f1","", -0.029) ;
  RooChebychev bkg3_f1("bkg3_f1","",m_phi_etap,RooArgList(c1_f1,c2_f1,c3_f1));

  //From final state 2:
  //Sidebands of phi and etap: chebychev.
  RooRealVar c1_f2("c1_f2","",   1.) ;
  RooRealVar c2_f2("c2_f2","",  -1.) ;
  RooRealVar c3_f2("c3_f2","",  -0.087) ;
  RooChebychev bkg3_f2("bkg3_f2","",m_phi_etap,RooArgList(c1_f2,c2_f2,c3_f2));

  //Background: 3-body
  //First final state
  RooRealVar p1("p1","p1", 0.0815);
  RooRealVar p2("p2","p2",-0.82596);
  RooRealVar p3("p3","p3", 0.1789);
  RooRealVar p4("p4","p4",-0.21711);
  RooRealVar p5("p5","p5", 0.0663);
  RooRealVar p6("p6","p6",-0.05046);
  RooMyBkg1 bkg1_f1_1("bkg1_f1_1","",m_phi_etap,p1,p2,p3,p4,p5,p6);//without smearing

  //Second final state
  RooRealVar p1_f2("p1_f2","", 0.0453);
  RooRealVar p2_f2("p2_f2","",-0.8175);
  RooRealVar p3_f2("p3_f2","", 0.1977);
  RooRealVar p4_f2("p4_f2","",-0.23829);
  RooRealVar p5_f2("p5_f2","", 0.0644);
  RooRealVar p6_f2("p6_f2","",-0.05164);
  RooMyBkg2 bkg1_f2_1("bkg1_f2_1","",m_phi_etap,p1_f2,p2_f2,p3_f2,p4_f2,p5_f2,p6_f2);//without smearing

  ////////////////////////
  //signal 
  ////////////////////////
  RooRealVar mean("mean","mean",2.1,1.9,2.6);//
  RooRealVar Gamma("Gamma","Gamma",0.14,0.001,0.5);//

  //|BW|^2 for first final state:
  RooMyBW1 sig1_f1("sig1_f1","",m_phi_etap,mean,Gamma);

  //|BW|^2 for second final state:
  RooMyBW2 sig1_f2("sig1_f2","",m_phi_etap,mean,Gamma);

  //From final state 1:
  //Construct double gaussian: gauss_f1
  RooRealVar meang("meang","meang",0);
  RooRealVar sigma1_f1("sigma1_f1","sigma1_f1",0.01113);
  RooRealVar sigma2_f1("sigma2_f1","sigma2_f1",0.004883);
  RooRealVar frac_f1("frac_f1","frac_f1",0.272);
  RooGaussian gauss1_f1("gauss1_f1","gauss1_f1",m_phi_etap,meang,sigma1_f1);
  RooGaussian gauss2_f1("gauss2_f1","gauss2_f1",m_phi_etap,meang,sigma2_f1);
  RooAddPdf gauss_f1 ("gauss_f1","gauss_f1",RooArgList(gauss1_f1,gauss2_f1),frac_f1);

  //Construct breit-wigner (x) gauss
  RooFFTConvPdf sig_f1("sig_f1","sig1_f1 (X) gauss_f1",m_phi_etap,sig1_f1,gauss_f1) ;

  //Construct Polynomial (x) gauss
  RooFFTConvPdf bkg1_f1("bkg1_f1","bkg1_f1_1 (X) gauss_f1",m_phi_etap,bkg1_f1_1,gauss_f1) ;

  //From final state 2:
  //Construct double gaussian: gauss_f2
  RooRealVar sigma1_f2("sigma1_f2","sigma1_f2",0.00833);
  RooRealVar sigma2_f2("sigma2_f2","sigma2_f2",0.01636);
  RooRealVar frac_f2("frac_f2","frac_f2",0.604);
  RooGaussian gauss1_f2("gauss1_f2","gauss1_f2",m_phi_etap,meang,sigma1_f2);
  RooGaussian gauss2_f2("gauss2_f2","gauss2_f2",m_phi_etap,meang,sigma2_f2);
  RooAddPdf gauss_f2 ("gauss_f2","gauss_f2",RooArgList(gauss1_f2,gauss2_f2),frac_f2);

  //Construct breit-wigner (x) gauss
  RooFFTConvPdf sig_f2 ("sig_f2","sig1_f2 (X) gauss_f2",m_phi_etap,sig1_f2,gauss_f2) ;

  //Construct Polynomial (x) gauss
  RooFFTConvPdf bkg1_f2 ("bkg1_f2","bkg1_f2_1 (X) gauss_f2",m_phi_etap,bkg1_f2_1,gauss_f2) ;

  ///////////////////////////////
  //Add signal and background
  ///////////////////////////////

  RooRealVar nsig_f1("nsig_f1","nsig_f1",500,0.,1000.);
  RooRealVar nbkg1_f1("nbkg1_f1","nbkg1_f1",1100,0,3000);
  RooRealVar nbkg2_f1("nbkg2_f1","nbkg2_f1",93);//
  RooRealVar nbkg3_f1("nbkg3_f1","nbkg3_f1",275.);//
  RooAddPdf model_f1("model_f1","model_f1",RooArgList(sig_f1,bkg1_f1,*bkg2_f1,bkg3_f1),RooArgList(nsig_f1,nbkg1_f1,nbkg2_f1,nbkg3_f1));

  RooRealVar nsig_f2("nsig_f2","nsig_f2",150,0.,500.);
  RooRealVar nbkg1_f2("nbkg1_f2","nbkg1_f2",300,0,1000);
  RooRealVar nbkg2_f2("nbkg2_f2","nbkg2_f2",29.);//
  RooRealVar nbkg3_f2("nbkg3_f2","nbkg3_f2",49.);//
  RooAddPdf model_f2("model_f2","model_f2",RooArgList(sig_f2,bkg1_f2,*bkg2_f2,bkg3_f2),RooArgList(nsig_f2,nbkg1_f2,nbkg2_f2,nbkg3_f2));

  //Creat index category and join samples
  ////////////////////////////////////////
  RooCategory sample("sample","sample") ;
  sample.defineType("state1") ;
  sample.defineType("state2") ;

  //Construct combined dataset
  RooDataSet combData("combData","combined data",RooArgSet(m_phi_etap),Index(sample),Import("state1",data_f1),Import("state2",data_f2)) ;//

  //Construct a simultaneous pdf
  RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;
  simPdf.addPdf(model_f1,"state1") ;
  simPdf.addPdf(model_f2,"state2") ;

  //Perform a simultaneous fit
  RooFitResult* r = simPdf.fitTo(combData,Extended(),RooFit::Save(kTRUE));

  Double_t FCN = r->minNll();//return minimized -log(L) value
  Int_t floatParms = (r->floatParsFinal()).getSize();

  //Draw graphs
  //////////////////////////

  //Final state 1
  RooPlot *frame1=m_phi_etap.frame();
  data_f1.plotOn(frame1,Binning(30));//
  model_f1.plotOn(frame1,Components(sig_f1),LineColor(kBlue),LineStyle(kDashed));
  model_f1.plotOn(frame1,Components(bkg1_f1),LineColor(kViolet),LineStyle(kDashed));
  model_f1.plotOn(frame1,Components(*bkg2_f1),LineColor(kOrange),LineStyle(kDashed));
  model_f1.plotOn(frame1,Components(bkg3_f1),LineColor(kGreen),LineStyle(kDashed));
  model_f1.plotOn(frame1,LineColor(kRed));

  //Final state 2
  RooPlot* frame2 = m_phi_etap.frame() ;
  data_f2.plotOn(frame2,Binning(30));//
  model_f2.plotOn(frame2,Components(sig_f2),LineColor(kBlue),LineStyle(kDashed));
  model_f2.plotOn(frame2,Components(bkg1_f2),LineColor(kViolet),LineStyle(kDashed));
  model_f2.plotOn(frame2,Components(*bkg2_f2),LineColor(kOrange),LineStyle(kDashed));
  model_f2.plotOn(frame2,Components(bkg3_f2),LineColor(kGreen),LineStyle(kDashed));
  model_f2.plotOn(frame2,LineColor(kRed));

  //Pull values
  RooPlot *frame1_=m_phi_etap.frame();
  frame1_->addObject(frame1->pullHist(),"p");//residHist()
  frame1_->SetMinimum(-5);
  frame1_->SetMaximum(+5);

  RooPlot *frame2_=m_phi_etap.frame();
  frame2_->addObject(frame2->pullHist(),"p");//residHist()
  frame2_->SetMinimum(-5);
  frame2_->SetMaximum(+5);

  /////////////////////////////
  //Control the style
  /////////////////////////////

  //With number of events.
  //Final state 1
  //Control the titles of x and y.
  frame1->GetXaxis()->SetTitle("M(#phi#eta') (GeV/c^{2})");
  frame1->GetYaxis()->SetTitle("Events/(20MeV/c^{2})");//

  //Center title.
  frame1->GetXaxis()->CenterTitle(true);
  frame1->GetYaxis()->CenterTitle(true);

  //Final state 2
  //Control the titles of x and y.
  frame2->GetXaxis()->SetTitle("M(#phi#eta') (GeV/c^{2})");
  frame2->GetYaxis()->SetTitle("Events/(20MeV/c^{2})");//

  //Center title.
  frame2->GetXaxis()->CenterTitle(true);
  frame2->GetYaxis()->CenterTitle(true);

  //Range of Yaxis.
  frame1->GetYaxis()->SetRangeUser(-9.99,110);//
  frame2->GetYaxis()->SetRangeUser(-9.99,45);//

  //Label and title size
  double label_size=0.06;//
  frame1->GetXaxis()->SetLabelSize(label_size);
  frame1->GetYaxis()->SetLabelSize(label_size);
  frame2->GetXaxis()->SetLabelSize(label_size);
  frame2->GetYaxis()->SetLabelSize(label_size);

  double title_size=0.08;//
  frame1->GetXaxis()->SetTitleSize(title_size);
  frame1->GetYaxis()->SetTitleSize(title_size);
  frame2->GetXaxis()->SetTitleSize(title_size);
  frame2->GetYaxis()->SetTitleSize(title_size);

  //With pull values.
  ////////////////////////////////
  //Title with pull values.
  frame1_->GetXaxis()->SetTitle("M(#phi#eta') (GeV/c^{2})");
  frame2_->GetXaxis()->SetTitle("M(#phi#eta') (GeV/c^{2})");

  frame1_->GetYaxis()->SetTitle("Pull value ");
  frame2_->GetYaxis()->SetTitle("Pull value ");

  //Title offset
  frame1_->GetYaxis()->SetTitleOffset(0.52);//
  frame2_->GetYaxis()->SetTitleOffset(0.52);//

  //Title division
  frame1_->GetYaxis()->SetNdivisions(505);
  frame2_->GetYaxis()->SetNdivisions(505);

  //Centertitle
  frame1_->GetXaxis()->CenterTitle(true);
  frame1_->GetYaxis()->CenterTitle(true);
  frame2_->GetXaxis()->CenterTitle(true);
  frame2_->GetYaxis()->CenterTitle(true);

  //Label and title size
  double label_size=0.13;//
  frame1_->GetXaxis()->SetLabelSize(label_size);
  frame1_->GetYaxis()->SetLabelSize(label_size);
  frame2_->GetXaxis()->SetLabelSize(label_size);
  frame2_->GetYaxis()->SetLabelSize(label_size);

  double title_size=0.16;//
  frame1_->GetXaxis()->SetTitleSize(title_size);
  frame1_->GetYaxis()->SetTitleSize(title_size);
  frame2_->GetXaxis()->SetTitleSize(title_size);
  frame2_->GetYaxis()->SetTitleSize(title_size);

  ////////////////////
  //Draw the result
  ////////////////////
  TCanvas *C1=new TCanvas("c1","",0,0,700*2/1.2,500*10/8.5/1.2);//
  C1->Divide(2,1);

  C1->cd(1);
  TPad *p11=new TPad("p11","",0.0,0.0,1.0,0.3);//down
  TPad *p12=new TPad("p12","",0.0,0.3,1.0,1.0);//up
    //Margin
  p11->SetBottomMargin(0.4);//
  p12->SetBottomMargin(0.0);

  p11->Draw();
  p12->Draw();
  p11->cd();
  frame1_->Draw();

  //Line
  TLine *l1=new TLine(frame1_->GetXaxis()->GetXmin(),4,frame1_->GetXaxis()->GetXmax(),4);
  l1->SetLineStyle(2);
  l1->Draw();

  TLine *l2=new TLine(frame1_->GetXaxis()->GetXmin(),2,frame1_->GetXaxis()->GetXmax(),2);
  l2->SetLineStyle(2);
  l2->Draw();

  TLine *l3=new TLine(frame1_->GetXaxis()->GetXmin(),0,frame1_->GetXaxis()->GetXmax(),0);
  l3->SetLineStyle(2);
  l3->Draw();

  TLine *l4=new TLine(frame1_->GetXaxis()->GetXmin(),-2,frame1_->GetXaxis()->GetXmax(),-2);
  l4->SetLineStyle(2);
  l4->Draw();

  TLine *l5=new TLine(frame1_->GetXaxis()->GetXmin(),-4,frame1_->GetXaxis()->GetXmax(),-4);
  l5->SetLineStyle(2);
  l5->Draw();

  p12->cd();
  frame1->Draw();

  //Prepare for the legend.
  TH1F *Func4=new TH1F("Func4","",100,1.9,2.6);//data

  TH1F *Func5=new TH1F("Func5","",100,1.9,2.6);//model
  Func5->SetLineColor(kRed);

  TH1F *Func6=new TH1F("Func6","",100,1.9,2.6);//sideband bkg
  Func6->SetLineColor(kGreen);
  Func6->SetLineStyle(kDashed);

  TH1F *Func8=new TH1F("Func8","",100,1.9,2.6);//sig
  Func8->SetLineColor(kBlue);
  Func8->SetLineStyle(kDashed);

  TH1F *Func9=new TH1F("Func9","",100,1.9,2.6);//3-body
  Func9->SetLineColor(kViolet);
  Func9->SetLineStyle(kDashed);

  TH1F *Func10=new TH1F("Func10","",100,1.9,2.6);//3-body
  Func10->SetLineColor(kOrange);
  Func10->SetLineStyle(kDashed);

  //Legend
  TLegend *leg = new TLegend (0.5877005,0.1995047,0.7881586,0.4816452,"");
  leg->SetFillColor(0);
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->AddEntry(Func4,"Data","lpe");
  leg->AddEntry(Func5,"Model","l");
  leg->AddEntry(Func8,"Signal","l");
  leg->AddEntry(Func9,"3-body bkg","l");
  leg->AddEntry(Func6,"Sid. bkg","l");
  leg->AddEntry(Func10,"#phi f_{0}(1500)","l");
  leg->Draw("same");

  // Chisq/ndf
  TLatex *ltx1 = new TLatex();
  ltx1->SetNDC(kTRUE);
  ltx1->SetTextColor(kPink);
  ltx1->SetTextFont(22);
  ltx1->SetTextSize(0.08);
  double chisq1=frame1->chiSquare(floatParms);
  ltx1->DrawLatex(0.72,0.85,Form("#chi^{2}/ndf=%.2f",chisq1));

  p11->Modified();
  p12->Modified();
 
  C1->cd(2);
  TPad *p21=new TPad("p21","",0.0,0.0,1.0,0.3);
  TPad *p22=new TPad("p22","",0.0,0.3,1.0,1.0);
    //Margin
  p21->SetBottomMargin(0.4);//
  p22->SetBottomMargin(0.0);

  p21->Draw();
  p22->Draw();
  p21->cd();
  frame2_->Draw();

  //Line
  TLine *l6=new TLine(frame2_->GetXaxis()->GetXmin(),4,frame2_->GetXaxis()->GetXmax(),4);
  l6->SetLineStyle(2);
  l6->Draw();

  TLine *l7=new TLine(frame2_->GetXaxis()->GetXmin(),2,frame2_->GetXaxis()->GetXmax(),2);
  l7->SetLineStyle(2);
  l7->Draw();

  TLine *l8=new TLine(frame2_->GetXaxis()->GetXmin(),0,frame2_->GetXaxis()->GetXmax(),0);
  l8->SetLineStyle(2);
  l8->Draw();

  TLine *l9=new TLine(frame2_->GetXaxis()->GetXmin(),-2,frame2_->GetXaxis()->GetXmax(),-2);
  l9->SetLineStyle(2);
  l9->Draw();

  TLine *l10=new TLine(frame2_->GetXaxis()->GetXmin(),-4,frame2_->GetXaxis()->GetXmax(),-4);
  l10->SetLineStyle(2);
  l10->Draw();
 
  p22->cd();
  frame2->Draw();

  // Chisq/ndf
  TLatex *ltx2 = new TLatex();
  ltx2->SetNDC(kTRUE);
  ltx2->SetTextColor(kPink);
  ltx2->SetTextFont(22);
  ltx2->SetTextSize(0.08);
  double chisq2=frame2->chiSquare(floatParms);
  ltx2->DrawLatex(0.72,0.85,Form("#chi^{2}/ndf=%.2f",chisq2));

  p21->Modified();
  p22->Modified();

  C1->Update();

  //Output the results
  cout<<"Likelihood with signal:"<<endl;
  cout<<"-lnL=              "<<FCN<<endl;
  cout<<"floatParms=        "<<floatParms<<endl;

  cout<<"Results after the fit:"<<endl;
  cout<<"Chisq_f1=          "<<frame1->chiSquare(floatParms)<<endl;
  cout<<"Chisq_f2=          "<<frame2->chiSquare(floatParms)<<endl;

  cout<<"mean=              "<<mean<<endl;
  cout<<"Gamma=             "<<Gamma<<endl;
  cout<<"nsig_f1=           "<<nsig_f1<<endl;
  cout<<"nsig_f2=           "<<nsig_f2<<endl;
  cout<<"nbkg1_f1=          "<<nbkg1_f1<<endl;
  cout<<"nbkg1_f2=          "<<nbkg1_f2<<endl;
  cout<<"nbkg2_f1=          "<<nbkg2_f1<<endl;
  cout<<"nbkg2_f2=          "<<nbkg2_f2<<endl;
  cout<<"nbkg3_f1=          "<<nbkg3_f1<<endl;
  cout<<"nbkg3_f2=          "<<nbkg3_f2<<endl;

  C1->Print("/Users/long/Desktop/Interference_no.pdf");
  cout<<"Finished!!!"<<endl;
}
