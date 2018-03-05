#include <TStyle.h>

void rootlogon(){
//Set style of graphs.
//By ShanWei.
//Modified by LongYunfei.

////////////////////////////////////
//Don't use in drawing Dalitz plot//
////////////////////////////////////

//The titles of x and y axis.
//Two plots one time (label: 0.04, title: 0.05) TitleOffset(1.5)
//Three plots one time (label: 0.05, title: 0.06) TitleOffset(1.5)
//One plot or four plots one time (label: 0.06, title: 0.08)
//Six plots one time (label: 0.08, title: 0.1)

//////////////////////////////
Double_t size_xy=0.06;
Double_t size_xy_title=0.08;
//////////////////////////////

gStyle->SetOptTitle(1);
gStyle->SetOptStat(1111111);
gStyle->SetOptFit(1111);

// BES style changed from BABAR style
TStyle *besStyle= new TStyle("BES","BES approved plots style");

// use plain black on white colors
besStyle->SetCanvasBorderMode(0);
besStyle->SetCanvasColor(0);
besStyle->SetFrameBorderMode(0);
besStyle->SetFrameBorderSize(3);
besStyle->SetFrameLineStyle(1);
besStyle->SetFrameLineWidth(2);
besStyle->SetFrameLineColor(0);
besStyle->SetPadBorderMode(0);
besStyle->SetPadColor(0);
besStyle->SetStatColor(0);

// set the paper & margin sizes
/////////////////////////////////////
besStyle->SetPaperSize(20,26);
besStyle->SetPadTopMargin(0.05);
besStyle->SetPadRightMargin(0.05);//
besStyle->SetPadBottomMargin(0.20);
besStyle->SetPadLeftMargin(0.20);//
/////////////////////////////////////

// use large Times-Roman fonts
besStyle->SetTextFont(62);
besStyle->SetTextSize(size_xy); 
besStyle->SetTitleFont(62,"x");  // set the all 3 axes title font
besStyle->SetTitleFont(62,"y");    // set the pad title font
besStyle->SetTitleFont(62,"z");    // set the pad title font
besStyle->SetLabelFont(62,"x");
besStyle->SetLabelFont(62,"y");
besStyle->SetLabelFont(62,"z");
besStyle->SetLabelSize(size_xy,"x");
besStyle->SetTitleSize(size_xy_title,"x");
besStyle->SetLabelSize(size_xy,"y");
besStyle->SetTitleSize(size_xy_title,"y");
besStyle->SetLabelSize(size_xy,"z");
besStyle->SetTitleSize(size_xy,"z");

// use bold lines and markers
besStyle->SetMarkerStyle(20);
besStyle->SetMarkerSize(1.4);//
besStyle->SetMarkerColor(1);
besStyle->SetLineWidth(2);
besStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes


// do not display any of the standard histogram decorations
besStyle->SetOptTitle(0);
besStyle->SetOptFit(0);

// put tick marks on top and RHS of plots
besStyle->SetPadTickX(1);
besStyle->SetPadTickY(1);

besStyle->SetHistLineWidth(2);

//**************************************
//some useful functions defined by Long.
//**************************************

void sethiststyle(TH1 *h){
  //Center title
  h->GetXaxis()->CenterTitle(true);
  h->GetYaxis()->CenterTitle(true);

  //TitleOffest
 // h->GetYaxis()->SetTitleOffset(iOffset);

  //Mass spectrum
  h->GetYaxis()->SetTitle(Form("Events/(%.0fMeV/c^{2})",h->GetBinWidth(1)*1000));
}

void centertitle(TH1 *h){
  //Center title
  h->GetXaxis()->CenterTitle(true);
  h->GetYaxis()->CenterTitle(true);
}

// Divisions
void setdivision(TH1 *h){
  //Ndivision
  h->GetXaxis()->SetNdivisions(505);
}

void setymass(TH1 *h){
  //Mass spectrum
  h->GetYaxis()->SetTitle(Form("Events/(%.0fMeV/c^{2})",h->GetBinWidth(1)*1000));
}

////////////////////////////////////////
//The relation between number and color
//  1  black
//  2  red
//  3  green
//  4  blue
//  5  yellow
//  6  violet
////////////////////////////////////////

// Draw arrow
void drawarrow(Double_t x,Double_t y1,Double_t y2,Double_t iColor=2){
  TArrow *ar=new TArrow(x,y1,x,y2,0.02,"<");
  ar->SetLineColor(iColor);
  ar->SetLineWidth(2);
  ar->Draw();
}

//Draw pave
//=*************************************
//The relation between number and color
//  1  black
//  2  red
//  3  green
//  4  blue
//  5  yellow
//  6  violet
//=**************************************
void drawpave(Double_t x1,Double_t y1,Double_t x2,Double_t y2,Double_t iColor){
  TPave *pave = new TPave(x1,y1,x2,y2,1,"br");
  pave->SetFillColor(0);
  pave->SetFillStyle(0);
  pave->SetLineWidth(3);
  pave->SetLineColor(iColor);
  pave->Draw();
}


// Write Entries,Mean and RMS
void box(Double_t x, Double_t y,TH1 *h,Char_t sName[]="Entries=%.0f",Char_t mName[]="Mean=%.4f",Char_t rName[]="RMS=%.4f",Double_t sizeTxt=0.08){
  h->SetStats(kFALSE);
  TLatex *ltx=new TLatex();
  ltx->SetNDC(kTRUE);
  ltx->SetTextColor(h->GetLineColor());
  ltx->SetTextFont(22);
  ltx->SetTextSize(sizeTxt);
  ltx->DrawLatex(x,y,Form(sName,h->GetEntries()));
  ltx->DrawLatex(x,y-0.1,Form(mName,h->GetMean()));
  ltx->DrawLatex(x,y-0.2,Form(rName,h->GetRMS()));
  gPad->Modified();
  gPad->Update();
}

// Write just entries.
void entry(Double_t x, Double_t y,TH1 *h,Char_t sName[]="Entries=%.0f",Double_t sizeTxt=0.08){
  h->SetStats(kFALSE);
  TLatex *ltx=new TLatex();
  ltx->SetNDC(kTRUE);
  ltx->SetTextColor(h->GetLineColor());
  ltx->SetTextFont(22);
  ltx->SetTextSize(sizeTxt);
  ltx->DrawLatex(x,y,Form(sName,h->GetEntries()));
  gPad->Modified();
  gPad->Update();
}

void txt(Double_t x, Double_t y, TString str, Int_t iColor=kBlack, Double_t fsize=0.08){
  TLatex *ltx = new TLatex();
  ltx->SetNDC(kTRUE);
  ltx->SetTextColor(iColor);
  ltx->SetTextFont(22);
  ltx->SetTextSize(fsize);
  ltx->DrawLatex(x,y,str.Data());
  gPad->Modified();
  gPad->Update();
}

//new TBrowser();
printf("Welcome to ROOT!\n");
printf("Start with new root style!\n");

gROOT->SetStyle("BES");
gROOT->ForceStyle();
}

