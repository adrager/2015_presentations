#include <iostream>
#include <vector>

#include "TArrow.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TStyle.h"



const int kNGausPar = 3;            //!< Number of parameters of Gaussian part of pdf


double mRMin;                        //!< Minimum of non-zero range of response pdf
double mRMax;                        //!< Maximum of non-zero range of response pdf
int    mNStepPar;                    //!< Number of parameters of step part of pdf
double mBinWidth;                    //!< Bin width
std::vector<double> mBinCenters;     //!< Centers of bins
std::vector<double> stepPar;         //!< Store step parameters



double gaussPart(double x, const std::vector<double>& par);
double stepPart(double x, const std::vector<double>& par, bool interpolate = true);
double binContent(double r, const std::vector<double>& par);
double interpolatedBinContent(double r, const std::vector<double>& par);
void setup(double min, double max, int n);

void plot_js_resofunc() {
  setup(0.3,1.3,7);
  
  std::vector<double> par;
  par.push_back(0.9);
  par.push_back(1.0);
  par.push_back(0.06);
  par.push_back(3E-2);
  par.push_back(5E-2);
  par.push_back(3E-1);
  par.push_back(6E-1);
  par.push_back(1);
  par.push_back(0.5);
  par.push_back(1E-1);


  std::cout << "Generating histograms... ";
  int nBins = 1000;
  double min = 0.8*mRMin;
  double max = 1.1*mRMax;
  TH1F * hStep = new TH1F("hStep",";Response R;Probability density",nBins,min,max);
  hStep->SetLineWidth(2);
  hStep->SetLineColor(4);
  TH1F * hStepInter = new TH1F("hStepInter",";Response R;Probability density",nBins,min,max);
  hStepInter->SetLineWidth(2);
  hStepInter->SetLineColor(4);
  TH1F * hStepInterFilled = new TH1F("hStepInterFilled",";Response R;Probability density",nBins,min,max);
  hStepInterFilled->SetLineWidth(2);
  hStepInterFilled->SetFillColor(4);
  hStepInterFilled->SetLineColor(4);
  TH1F * hGauss = new TH1F("hGauss",";Response R;Probability density",nBins,min,max);
  hGauss->SetLineWidth(2);
  hGauss->SetLineColor(3);
  TH1F * hResp1 = new TH1F("hResp1",";Response R;Probability density",nBins,min,max);
  hResp1->SetLineWidth(2);
  hResp1->SetLineColor(2);
  TH1F * hResp2 = new TH1F("hResp2",";Response R;Probability density",nBins,min,max);
  hResp2->SetLineWidth(2);
  hResp2->SetLineColor(2);
  std::cout << "ok\n";


  std::cout << "Filling histograms... ";
  for(int bin = 1; bin <= hStep->GetNbinsX(); bin++) {
    double x = hStep->GetBinCenter(bin);

    hStep->SetBinContent(bin,stepPart(x,par,false));
    hStepInter->SetBinContent(bin,stepPart(x,par,true));
    hGauss->SetBinContent(bin,gaussPart(x,par));
    hResp1->SetBinContent(bin,gaussPart(x,par)+stepPart(x,par,false));
    hResp2->SetBinContent(bin,gaussPart(x,par)+stepPart(x,par,true));
    if( x > mRMin+3*mBinWidth && x < mRMin+4*mBinWidth ) {
      hStepInterFilled->SetBinContent(bin,stepPart(x,par,true));
    }
    else {
      hStepInterFilled->SetBinContent(bin,0.);
    }
  }

  double norm = hResp2->Integral("width");
  hStep->Scale(1./norm);
  hStepInter->Scale(1./norm);
  hGauss->Scale(1./norm);
  hResp1->Scale(1./norm);
  hResp2->Scale(1./norm);
  
  hGauss->GetYaxis()->SetRangeUser(2E-3,10);
  hStep->GetYaxis()->SetRangeUser(2E-3,10);
  hStepInter->GetYaxis()->SetRangeUser(2E-3,10);
  hStepInterFilled->GetYaxis()->SetRangeUser(2E-3,10);
  hResp1->GetYaxis()->SetRangeUser(2E-3,10);
  hResp2->GetYaxis()->SetRangeUser(2E-3,10);
  std::cout << "ok\n";


  std::cout << "Plotting histograms... ";
  gStyle->SetOptStat(0);

  // Sum of both parts
  TCanvas * can1 = new TCanvas("can1","",1200,400);
  can1->Divide(3,1);
  can1->cd(1);
  hStep->Draw();
  hGauss->Draw("same");
  can1->cd(1)->SetLogy();

  can1->cd(2);
  TPaveText * labelSum = new TPaveText(0.08,0.55,0.88,0.65,"NDC");
  labelSum->SetFillColor(0);
  labelSum->SetTextFont(42);
  labelSum->SetTextSize(0.1);
  labelSum->SetBorderSize(0);
  labelSum->AddText("c #upoint #color[3]{G(R)} + (1-c) #upoint #color[4]{S(R)}");
  labelSum->Draw();

  TArrow * arrow = new TArrow(0.05,0.5,0.97,0.5,0.03,"|>");
  arrow->SetLineWidth(2);
  arrow->SetLineColor(1);
  arrow->SetAngle(30);
  arrow->Draw("");

  can1->cd(3);
  hResp1->Draw();
  can1->cd(3)->SetLogy();

  // Step function
  TCanvas * can2 = new TCanvas("can2","",1200,400);

  std::vector<TPad*> cPads;
  cPads.push_back(new TPad("cPad0","",0.00,0.00,0.33,1.00));
  cPads.push_back(new TPad("cPad1","",0.33,0.00,0.66,1.00));
  cPads.push_back(new TPad("cPad2","",0.66,0.00,1.00,1.00));
  cPads.push_back(new TPad("cGlobalPad","",0.,0.,1.,1.));

  for(unsigned int i = 0; i < cPads.size(); i++)
    {
      cPads.at(i)->SetFillStyle(0);
      cPads.at(i)->SetFrameFillColor(10);
      cPads.at(i)->SetFrameBorderMode(0);
      can2->cd();
      if( i < 3 ) cPads.at(i)->Draw();
    }


  can2->cd();
  cPads.at(0)->cd()->SetLogy();
  hStep->Draw();
  hGauss->Draw("same");

  can2->cd(1);
  cPads.at(1)->cd()->SetLogy();
  hStepInterFilled->Draw();
  hGauss->Draw("same");
  hStepInter->Draw("same");

  can2->cd();
  cPads.at(2)->cd()->SetLogy();
  hResp2->Draw();

  std::vector<TArrow*> connArrow;
  connArrow.push_back(new TArrow(0.3,0.6,0.45,0.65,0.03,"|>"));
  connArrow.push_back(new TArrow(0.63,0.65,0.78,0.6,0.03,"|>"));
  for(size_t i = 0; i < connArrow.size(); i++) {
    connArrow.at(i)->SetLineWidth(6);
    connArrow.at(i)->SetLineColor(2);
    connArrow.at(i)->SetFillColor(2);
    connArrow.at(i)->SetAngle(30);
  }


  can2->cd();
  cPads.at(3)->Draw("same");
  cPads.at(3)->cd();
  connArrow.at(0)->Draw();
  connArrow.at(1)->Draw();

  std::cout << "ok\n";
}


void setup(double min, double max, int n) {
  std::cout << "Setting up quantities... ";
  mRMin      = min;
  mRMax      = max;
  mNStepPar  = n;
  mBinWidth  = (mRMax - mRMin)/mNStepPar;
  for(int i = 0; i < mNStepPar+1; i++) {
    mBinCenters.push_back( mRMin + (0.5+i)*mBinWidth );
  }
  assert( 0.0 <= mRMin && mRMin < mRMax );
  std::cout << "ok\n";
}



//!  Probability density from Gaussian part
// ------------------------------------------------------------------------
double gaussPart(double x, const std::vector<double>& par) {
  double c = par.at(0);
  double u = par.at(1);
  double s = par.at(2);

  return c / sqrt(2* M_PI ) / s * exp(-pow((x - u) / s, 2) / 2);
}



//!  Probability density from step part
// ------------------------------------------------------------------------
double stepPart(double x, const std::vector<double>& par, bool interpolate) {
  std::vector<double> stepPar = std::vector<double>(mNStepPar,1.);
  double c    = par.at(0);
  double norm = 0.;
  for(int i = 0; i < mNStepPar; i++) {
    stepPar.at(i) = par.at(kNGausPar+i);
    if( stepPar.at(i) < 0. ) stepPar.at(i) = 0.;
    norm += stepPar.at(i);
  }
  norm = (1.-c)/norm/mBinWidth;
    
  double p = 0.;
  if( interpolate ) p = norm*interpolatedBinContent(x,stepPar);
  else              p = norm*binContent(x,stepPar);

  return p;
}



// ------------------------------------------------------------------------
double binContent(double r, const std::vector<double>& par) {
  double p = 0.;

  // Check that r is in range of binning +/- 0.5*mBinWidth
  if( r > mRMin  &&  r < mRMax ) {
    // Find index i of the bin
    int i = 0;                              
    while( r > mRMin + (1+i)*mBinWidth ) i++;
    assert( i < mNStepPar );
    p = par.at(i);
  }

  return p;
}



//!  \brief Linear interpolation between two bins
//!
//!  Returns the linearly weighted mean value of two adjacent bin
//!  contents, where the two bins are the ones whose bin centers
//!  are closest to r. The contents are weighted with the distance
//!  of r from the corresponding bin center. The interpolation
//!  assumes mNStepPar bins. At the left and right edges of the
//!  binned range, the interpolation is done assuming 0 bin content
//!  outside the binned range.
//!
//!  \param r   Response
//!  \return    Interpolated bin content
// ------------------------------------------------------------------------
double interpolatedBinContent(double r, const std::vector<double>& par) {
  double p = 0.;     // The interpolated parameter

  // Check that r is in range of binning +/- 0.5*mBinWidth
  if( r > mRMin - 0.5*mBinWidth  &&  r < mRMax + 0.5*mBinWidth ) {
    // Find index i of the next larger bin center
    unsigned int i = 0;                              
    while( r > mBinCenters.at(i) ) i++;
    assert( i < mBinCenters.size() );

    double dx  = mBinCenters.at(i) - r;  // Distance to next larger bin center
    dx        /= mBinWidth;              // dx relative to bin width
    double a   = 0.;
    double b   = 0.;

    // Find bin contents to be interpolated
    if( i == 0 ) { // Left edge
      a = 0.;
      b = par.front();
    }
    else if( i == mBinCenters.size()-1 ) { // Right edge
      a = par.back();
      b = 0.;
    }
    else { // Central part
      a = par.at(i-1);
      b = par.at(i);
    }

    // Interpolate contents
    p = dx * a + (1-dx) * b;
  }

  return p;
}
