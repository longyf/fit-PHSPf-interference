#include "components.C"
#include "components_f2.C"

using namespace RooFit;

void fit_intf_constructive_with_pull(){
// fit to M(phi etap) considering interference between resonance and 3-body bkg.

  gSystem->Load("libRooFit");
  gROOT->ProcessLine(".L RooMySig.cxx++");
  gROOT->ProcessLine(".L RooMySig_f2.cxx++");

  double low_=1.96;//
  double up_=2.56;//

  RooRealVar m_phi_etap("m_phi_etap","m_phi_etap",low_,up_);

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

  //Background( sideband )
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

  //Signal+3body+intf: f1
  ////////////////////////
  RooRealVar mean("mean","mean",2.11,1.9,2.6);//
  RooRealVar Gamma("Gamma","Gamma",0.150,0.00001,0.5);//
  RooRealVar phi("phi","phi",-1.7,-3.2,0.);//

  RooRealVar ratio_f1("ratio_f1","ratio_f1",10.,0.1,100.);//

    //polynomial (3body)
  RooRealVar p1("p1","p1", 0.0815);
  RooRealVar p2("p2","p2",-0.82596);
  RooRealVar p3("p3","p3", 0.1789);
  RooRealVar p4("p4","p4",-0.21711);
  RooRealVar p5("p5","p5", 0.0663);
  RooRealVar p6("p6","p6",-0.05046);
  RooMySig sig1_f1("sig1_f1","",m_phi_etap,mean,Gamma,phi,ratio_f1,p1,p2,p3,p4,p5,p6);

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

  //Signal+3body+intf: f2
  ////////////////////////
  RooRealVar ratio_f2("ratio_f2","ratio_f2",10.,0.1,100.);//

    //polynomial (3body)
  RooRealVar p1_f2("p1_f2","", 0.0453);
  RooRealVar p2_f2("p2_f2","",-0.8175);
  RooRealVar p3_f2("p3_f2","", 0.1977);
  RooRealVar p4_f2("p4_f2","",-0.23829);
  RooRealVar p5_f2("p5_f2","", 0.0644);
  RooRealVar p6_f2("p6_f2","",-0.05164);
  RooMySig_f2 sig1_f2("sig1_f2","",m_phi_etap,mean,Gamma,phi,ratio_f2,p1_f2,p2_f2,p3_f2,p4_f2,p5_f2,p6_f2);

  //From final state 2:
  //Construct double gaussian: gauss_f2
  RooRealVar sigma1_f2("sigma1_f2","sigma1_f2",0.00833);
  RooRealVar sigma2_f2("sigma2_f2","sigma2_f2",0.01636);
  RooRealVar frac_f2("frac_f2","frac_f2",0.604);
  RooGaussian gauss1_f2("gauss1_f2","gauss1_f2",m_phi_etap,meang,sigma1_f2);
  RooGaussian gauss2_f2("gauss2_f2","gauss2_f2",m_phi_etap,meang,sigma2_f2);
  RooAddPdf gauss_f2 ("gauss_f2","gauss_f2",RooArgList(gauss1_f2,gauss2_f2),frac_f2);

  //Construct breit-wigner (x) gauss
  RooFFTConvPdf sig_f2 ("sig_f2","sig1_f2 (X) gauss_f2",m_phi_etap,sig1_f2,gauss_f2) ;//

  //Add signal and background
  ///////////////////////////////
  RooRealVar nsig_f1("nsig_f1","nsig_f1",1700,500.,3000.);
  RooRealVar nbkg2_f1("nbkg2_f1","nbkg2_f1",93.);
  RooRealVar nbkg3_f1("nbkg3_f1","nbkg3_f1",275.);
  RooAddPdf model_f1("model_f1","model_f1",RooArgList(sig_f1,*bkg2_f1,bkg3_f1),RooArgList(nsig_f1,nbkg2_f1,nbkg3_f1));

  RooRealVar nsig_f2("nsig_f2","nsig_f2",500.,100.,1000.);
  RooRealVar nbkg2_f2("nbkg2_f2","nbkg2_f2",29.);
  RooRealVar nbkg3_f2("nbkg3_f2","nbkg3_f2",49.);
  RooAddPdf model_f2("model_f2","model_f2",RooArgList(sig_f2,*bkg2_f2,bkg3_f2),RooArgList(nsig_f2,nbkg2_f2,nbkg3_f2));

  //Creat index category and join samples
  ////////////////////////////////////////
  RooCategory sample("sample","sample") ;
  sample.defineType("state1") ;
  sample.defineType("state2") ;

  //Construct combined dataset
  RooDataSet combData("combData","combined data",RooArgSet(m_phi_etap),Index(sample),Import("state1",data_f1),Import("state2",data_f2)) ;

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
  data_f1.plotOn(frame1,Binning(30));
  model_f1.plotOn(frame1,Components(sig_f1),LineColor(kGreen));
  model_f1.plotOn(frame1,Components(*bkg2_f1),LineColor(kOrange),LineStyle(kDashed));
  model_f1.plotOn(frame1,Components(bkg3_f1),LineColor(kGreen),LineStyle(kDashed));
  model_f1.plotOn(frame1,LineColor(kRed));

  //Final state 2
  RooPlot* frame2 = m_phi_etap.frame() ;
  data_f2.plotOn(frame2,Binning(30));//
  model_f2.plotOn(frame2,Components(sig_f2),LineColor(kGreen));
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
  frame1->GetYaxis()->SetTitle("Events/(20MeV/c^{2})");

  //Center title.
  frame1->GetXaxis()->CenterTitle(true);
  frame1->GetYaxis()->CenterTitle(true);

  //Final state 2
  //Control the titles of x and y.
  frame2->GetXaxis()->SetTitle("M(#phi#eta') (GeV/c^{2})");
  frame2->GetYaxis()->SetTitle("Events/(20MeV/c^{2})");

  //Center title.
  frame2->GetXaxis()->CenterTitle(true);
  frame2->GetYaxis()->CenterTitle(true);

  //Range of Yaxis.
  frame1->GetYaxis()->SetRangeUser(-9.99,110);//
  frame2->GetYaxis()->SetRangeUser(-9.99,45);//

  //Label and title size
  double label_size=0.06;
  frame1->GetXaxis()->SetLabelSize(label_size);
  frame1->GetYaxis()->SetLabelSize(label_size);
  frame2->GetXaxis()->SetLabelSize(label_size);
  frame2->GetYaxis()->SetLabelSize(label_size);

  double title_size=0.08;
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
  frame1_->GetYaxis()->SetTitleOffset(0.52);
  frame2_->GetYaxis()->SetTitleOffset(0.52);

  //Title division
  frame1_->GetYaxis()->SetNdivisions(505);
  frame2_->GetYaxis()->SetNdivisions(505);

  //Centertitle
  frame1_->GetXaxis()->CenterTitle(true);
  frame1_->GetYaxis()->CenterTitle(true);
  frame2_->GetXaxis()->CenterTitle(true);
  frame2_->GetYaxis()->CenterTitle(true);

  //Label and title size
  double label_size=0.13;
  frame1_->GetXaxis()->SetLabelSize(label_size);
  frame1_->GetYaxis()->SetLabelSize(label_size);
  frame2_->GetXaxis()->SetLabelSize(label_size);
  frame2_->GetYaxis()->SetLabelSize(label_size);

  double title_size=0.16;
  frame1_->GetXaxis()->SetTitleSize(title_size);
  frame1_->GetYaxis()->SetTitleSize(title_size);
  frame2_->GetXaxis()->SetTitleSize(title_size);
  frame2_->GetYaxis()->SetTitleSize(title_size);

  ///////////////////////////
  //draw sig, bkg, intf (f1)
  ///////////////////////////
  double Pr0=mean.getVal();
  double Pr1=Gamma.getVal();
  double Pr2=phi.getVal();
  double Pr3=ratio_f1.getVal();
  double Pr4=p1.getVal();
  double Pr5=p2.getVal();
  double Pr6=p3.getVal();
  double Pr7=p4.getVal();
  double Pr8=p5.getVal();
  double Pr9=p6.getVal();

  double scale_f=16.5;//S1

  //all
  TF1 *Func0=new TF1("Func0", components,1.9,2.6,14);
  Func0->SetParameter(0, Pr0);
  Func0->SetParameter(1, Pr1);
  Func0->SetParameter(2, Pr2);
  Func0->SetParameter(3, Pr3);
  Func0->SetParameter(4, Pr4);
  Func0->SetParameter(5, Pr5);
  Func0->SetParameter(6, Pr6);
  Func0->SetParameter(7, Pr7);
  Func0->SetParameter(8, Pr8);
  Func0->SetParameter(9, Pr9);
  Func0->SetParameter(10, 1.0);
  Func0->SetParameter(11, 1.0);
  Func0->SetParameter(12, 1.0);
  Func0->SetParameter(13, scale_f);

  //signal
  TF1 *Func1=new TF1("Func1", components,1.9,2.6,14);
  Func1->SetParameter(0, Pr0);
  Func1->SetParameter(1, Pr1);
  Func1->SetParameter(2, Pr2);
  Func1->SetParameter(3, Pr3);
  Func1->SetParameter(4, Pr4);
  Func1->SetParameter(5, Pr5);
  Func1->SetParameter(6, Pr6);
  Func1->SetParameter(7, Pr7);
  Func1->SetParameter(8, Pr8);
  Func1->SetParameter(9, Pr9);
  Func1->SetParameter(10, 1.0);
  Func1->SetParameter(11, 0.0);
  Func1->SetParameter(12, 0.0);
  Func1->SetParameter(13, scale_f);

  Func1->SetLineColor(kBlue);
  Func1->SetLineStyle(kDashed);

  //background
  TF1 *Func2=new TF1("Func2", components,1.9,2.6,14);
  Func2->SetParameter(0, Pr0);
  Func2->SetParameter(1, Pr1);
  Func2->SetParameter(2, Pr2);
  Func2->SetParameter(3, Pr3);
  Func2->SetParameter(4, Pr4);
  Func2->SetParameter(5, Pr5);
  Func2->SetParameter(6, Pr6);
  Func2->SetParameter(7, Pr7);
  Func2->SetParameter(8, Pr8);
  Func2->SetParameter(9, Pr9);
  Func2->SetParameter(10, 0.0);
  Func2->SetParameter(11, 1.0);
  Func2->SetParameter(12, 0.0);
  Func2->SetParameter(13, scale_f);

  Func2->SetLineColor(kViolet);
  Func2->SetLineStyle(kDashed);

  //intf
  TF1 *Func3=new TF1("Func3", components,1.9,2.6,14);
  Func3->SetParameter(0, Pr0);
  Func3->SetParameter(1, Pr1);
  Func3->SetParameter(2, Pr2);
  Func3->SetParameter(3, Pr3);
  Func3->SetParameter(4, Pr4);
  Func3->SetParameter(5, Pr5);
  Func3->SetParameter(6, Pr6);
  Func3->SetParameter(7, Pr7);
  Func3->SetParameter(8, Pr8);
  Func3->SetParameter(9, Pr9);
  Func3->SetParameter(10, 0.0);
  Func3->SetParameter(11, 0.0);
  Func3->SetParameter(12, 1.0);
  Func3->SetParameter(13, scale_f);

  Func3->SetLineColor(kPink);
  Func3->SetLineStyle(4);

  ////////////////
  //f2
  ////////////////

  Pr3=ratio_f2.getVal();
  Pr4=p1_f2.getVal();
  Pr5=p2_f2.getVal();
  Pr6=p3_f2.getVal();
  Pr7=p4_f2.getVal();
  Pr8=p5_f2.getVal();
  Pr9=p6_f2.getVal();
  double scale_f2=13.5;//S2

  //all
  TF1 *Func0_2=new TF1("Func0_2", components_f2,1.9,2.6,14);
  Func0_2->SetParameter(0, Pr0);
  Func0_2->SetParameter(1, Pr1);
  Func0_2->SetParameter(2, Pr2);
  Func0_2->SetParameter(3, Pr3);
  Func0_2->SetParameter(4, Pr4);
  Func0_2->SetParameter(5, Pr5);
  Func0_2->SetParameter(6, Pr6);
  Func0_2->SetParameter(7, Pr7);
  Func0_2->SetParameter(8, Pr8);
  Func0_2->SetParameter(9, Pr9);
  Func0_2->SetParameter(10, 1.0);
  Func0_2->SetParameter(11, 1.0);
  Func0_2->SetParameter(12, 1.0);
  Func0_2->SetParameter(13, scale_f2);

  //signal
  TF1 *Func1_2=new TF1("Func1_2", components_f2,1.9,2.6,14);
  Func1_2->SetParameter(0, Pr0);
  Func1_2->SetParameter(1, Pr1);
  Func1_2->SetParameter(2, Pr2);
  Func1_2->SetParameter(3, Pr3);
  Func1_2->SetParameter(4, Pr4);
  Func1_2->SetParameter(5, Pr5);
  Func1_2->SetParameter(6, Pr6);
  Func1_2->SetParameter(7, Pr7);
  Func1_2->SetParameter(8, Pr8);
  Func1_2->SetParameter(9, Pr9);
  Func1_2->SetParameter(10, 1.0);
  Func1_2->SetParameter(11, 0.0);
  Func1_2->SetParameter(12, 0.0);
  Func1_2->SetParameter(13, scale_f2);

  Func1_2->SetLineColor(kBlue);
  Func1_2->SetLineStyle(kDashed);

  //background
  TF1 *Func2_2=new TF1("Func2_2", components_f2,1.9,2.6,14);
  Func2_2->SetParameter(0, Pr0);
  Func2_2->SetParameter(1, Pr1);
  Func2_2->SetParameter(2, Pr2);
  Func2_2->SetParameter(3, Pr3);
  Func2_2->SetParameter(4, Pr4);
  Func2_2->SetParameter(5, Pr5);
  Func2_2->SetParameter(6, Pr6);
  Func2_2->SetParameter(7, Pr7);
  Func2_2->SetParameter(8, Pr8);
  Func2_2->SetParameter(9, Pr9);
  Func2_2->SetParameter(10, 0.0);
  Func2_2->SetParameter(11, 1.0);
  Func2_2->SetParameter(12, 0.0);
  Func2_2->SetParameter(13, scale_f2);

  Func2_2->SetLineColor(kViolet);
  Func2_2->SetLineStyle(kDashed);

  //intf
  TF1 *Func3_2=new TF1("Func3_2", components_f2,1.9,2.6,14);
  Func3_2->SetParameter(0, Pr0);
  Func3_2->SetParameter(1, Pr1);
  Func3_2->SetParameter(2, Pr2);
  Func3_2->SetParameter(3, Pr3);
  Func3_2->SetParameter(4, Pr4);
  Func3_2->SetParameter(5, Pr5);
  Func3_2->SetParameter(6, Pr6);
  Func3_2->SetParameter(7, Pr7);
  Func3_2->SetParameter(8, Pr8);
  Func3_2->SetParameter(9, Pr9);
  Func3_2->SetParameter(10, 0.0);
  Func3_2->SetParameter(11, 0.0);
  Func3_2->SetParameter(12, 1.0);
  Func3_2->SetParameter(13, scale_f2);

  Func3_2->SetLineColor(kPink);
  Func3_2->SetLineStyle(4);

  ////////////////////
  //Draw the result
  ////////////////////
  TCanvas *C1=new TCanvas("c1","",0,0,700*2/1.2,500*10/8.5/1.2);
  C1->Divide(2,1);

  C1->cd(1);
  TPad *p11=new TPad("p11","",0.0,0.0,1.0,0.3);//down
  TPad *p12=new TPad("p12","",0.0,0.3,1.0,1.0);//up
    //Margin
  p11->SetBottomMargin(0.4);
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
//  Func0->Draw("same");//FUCK
  Func1->Draw("same");
  Func2->Draw("same");
  Func3->Draw("same");

  //Prepare for the legend.
  TH1F *Func4=new TH1F("Func4","",100,1.9,2.6);//data

  TH1F *Func5=new TH1F("Func5","",100,1.9,2.6);//model
  Func5->SetLineColor(kRed);

  TH1F *Func6=new TH1F("Func6","",100,1.9,2.6);//sideband bkg
  Func6->SetLineColor(kGreen);
  Func6->SetLineStyle(kDashed);

  TH1F *Func7=new TH1F("Func7","",100,1.9,2.6);//model-sideband
  Func7->SetLineColor(kGreen);

  TH1F *Func8=new TH1F("Func8","",100,1.9,2.6);//f0(1500)
  Func8->SetLineColor(kOrange);
  Func8->SetLineStyle(kDashed);

  //Legend
  TLegend *leg = new TLegend (0.6127577,0.1869651,0.8042669,0.4628358,"");
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->AddEntry(Func4,"Data","lpe");
  leg->AddEntry(Func5,"Model","l");
//  leg->AddEntry(Func7,"Model - Sid. - f_{0}(1500)","l");
  leg->AddEntry(Func1,"Signal","l");
  leg->AddEntry(Func2,"3-body bkg","l");
  leg->AddEntry(Func3,"Interference","l");
  leg->AddEntry(Func6,"Sid. bkg","l");
  leg->AddEntry(Func8,"f_{0}(1500)","l");
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
  p21->SetBottomMargin(0.4);
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
//  Func0_2->Draw("same");//FUCK
  Func1_2->Draw("same");
  Func2_2->Draw("same");
  Func3_2->Draw("same");

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
  cout<<"nbkg3_f1=          "<<nbkg3_f1<<endl;
  cout<<"nbkg3_f2=          "<<nbkg3_f2<<endl;
  cout<<"phi=            "<<phi<<endl;
  cout<<"ratio_f1=          "<<ratio_f1<<endl;
  cout<<"ratio_f2=          "<<ratio_f2<<endl;

  ///////////////////
  //scale factor of f1
  double N_tot_1=Func0->Integral(low_,up_);//
  double SCALE_F1=(1959.-275.-93)/N_tot_1;

  double N_sig_1=Func1->Integral(low_,up_);
  double N_bkg_1=Func2->Integral(low_,up_);
  double N_intf_1=Func3->Integral(low_,up_);

  cout<<"N_sig_1="<<N_sig_1*SCALE_F1<<" +- "<<TMath::Sqrt(N_sig_1*SCALE_F1)<<endl;
  cout<<"N_bkg_1="<<N_bkg_1*SCALE_F1<<" +- "<<TMath::Sqrt(N_bkg_1*SCALE_F1)<<endl;
  cout<<"N_intf_1="<<N_intf_1*SCALE_F1<<" +- "<<TMath::Sqrt(TMath::Abs(N_intf_1*SCALE_F1))<<endl;

  //scale factor of f2
  double N_tot_2=Func0_2->Integral(low_,up_);//
  double SCALE_F2=(559.-49.-29.)/N_tot_2;

  double N_sig_2=Func1_2->Integral(low_,up_);
  double N_bkg_2=Func2_2->Integral(low_,up_);
  double N_intf_2=Func3_2->Integral(low_,up_);

  cout<<"N_sig_2="<<N_sig_2*SCALE_F2<<" +- "<<TMath::Sqrt(N_sig_2*SCALE_F2)<<endl;
  cout<<"N_bkg_2="<<N_bkg_2*SCALE_F2<<" +- "<<TMath::Sqrt(N_bkg_2*SCALE_F2)<<endl;
  cout<<"N_intf_2="<<N_intf_2*SCALE_F2<<" +- "<<TMath::Sqrt(TMath::Abs(N_intf_2*SCALE_F2))<<endl;
  ///////////////////////

//  C1->Print("/Users/long/Desktop/Interference_constructive.pdf");
  cout<<"Finished!!!"<<endl;
}
