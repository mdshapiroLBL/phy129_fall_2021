// Filename:  siliconResolution.C
// Author:    Marjorie Shapiro
// Date       September 2016
//
// This macro is an example showing how to write a simple
// Monte Carlo using ROOT.
//
// To run this macro startup ROOT and then do:
//
// .L siliconResolution.C++
//  siliconResoltion(10000,"siliconResolution.png");
//
//
// To exit ROOT type
// .q
//
//
// Note:  Comments  in ++ start with two backslashes
// You need to include at the top of your file the header files
// defining any classes you use.
# 
#include <iostream>      // This is the C++ class for print things
#include <string>        // This is the C++ class for handling strings
                         // All Root classes have names starting with T
#include <TROOT.h>
#include "TRandom.h"     // This is the random number class in Root
#include "TH1F.h"        // A profile plot gives the mean and variance per bin
#include "TCanvas.h"
#include "TMath.h"

using std::cout;
using std::endl;
using std::string;

// This declares our class
// The input is the number of events, the min and max momentum 
// You need to declare for each variable if it is an integer, float or double
// the function type void means that EMShower doesnt' return a value 
//
// General rules:
//   All statements end with a ;
//   All variables and functions must be declared as integer, float or double
//   You can initialize variables when you declare them. 
//
// We've set this up so that the number of events, Pmin and Pmax are passed
// to the code but the lifetime is set in the code
//
void siliconResolution(int numEvt, string fileName) {

  // Create histograms 
  // Arguments are: Name of histogram, Title, number of bins,min, max
  TH1F hista("Parta","Part a ",100,-1.0,1.0);
  hista.GetXaxis()->SetTitle("(x_{meas}-x_{true})");
  hista.GetYaxis()->SetTitle("Number of Entries");

  TH1F histb("Partb","Part b ",100,-1.0,1.0);
  histb.GetXaxis()->SetTitle("(x_{meas}-x_{true})");
  histb.GetYaxis()->SetTitle("Number of Entries");

  TH1F histc("Partc.i","Part c.i ",100,-1.0,1.0);
  histc.GetXaxis()->SetTitle("(x_{meas}-x_{true})");
  histc.GetYaxis()->SetTitle("Number of Entries");

  TH1F histd("Partc.ii","Part c.iii ",100,-1.0,1.0);
  histd.GetXaxis()->SetTitle("(x_{meas}-x_{true})");
  histd.GetYaxis()->SetTitle("Number of Entries");

  // This next command makes Root calcuate the mean and sigma using the
  // data rather than the center of the bins.
  hista.Sumw2();
  histb.Sumw2();
  histc.Sumw2();
  histd.Sumw2();

  // We need to create an object that is a random number generator
  TRandom* random = new TRandom();
  random->SetSeed(12345);

  cout << "About to generate " << numEvt << " events "  << endl;

  double sigmaM = 1.0;
  double thresh1 = 0.2;
  double thresh2 = 0.1;
  double sigmaN1 = 0.05;
  double sigmaN2 = 0.025;

  // We'll make our histograms in units of strip pitch.  So, our
  // beam is made of particles with -0.5<x<0.5
 for (int i=0; i<numEvt; i++) {
      // TRandom::Uniform(x) throws a flat distribution from -0.5 to 0.5
    double xtrue = random->Uniform(1.0)-0.5; 
    double xmeasa = 0;  // Measurement at center of strip for binary readout
    hista.Fill(xmeasa-xtrue);

    double xmeasb;
    if(fabs(xtrue)<1.0/3.0)  xmeasb=0;
    else if(xtrue>0) xmeasb = 0.5;
    else xmeasb=-0.5;
    histb.Fill(xmeasb-xtrue);

    double xTot1=0;
    double chTot1=0;
    double xTot2=0;
    double chTot2=0;
    for(int strip=0; strip<7; strip++) {
      // The strip edges need to be calculated relative to xtrue
      double lowEdge=-3.5+(double) strip-xtrue;
      double highEdge=-3.5+(double) strip+1-xtrue;
      double int0, int1;

      // TMath::Erf(x) is Computation of the error function erf(x).
      //   Erf(x) = (2/sqrt(pi)) Integral(exp(-t^2))dt between 0 and x
      // TMath::Erfc(x) is Computation the complementary error function erfc(x)
      //   Erfc(x) = (2/sqrt(pi)) Integral(exp(-t^2))dt between x and infinity
      if(lowEdge<0) int0 = 0.5*TMath::Erfc(fabs(lowEdge));
      else int0 = 0.5+0.5*TMath::Erf(lowEdge);
      if(highEdge>0) int1 = 0.5*TMath::Erfc(fabs(highEdge));
      else int1 = 0.5+0.5*TMath::Erf(fabs(highEdge));

      // ch is the charge on the strip
      double noise1 = random->Gaus(0.0,sigmaN1);  
      double ch1 = 1.0-int0-int1+noise1;
      double noise2 = random->Gaus(0.0,sigmaN2);  
      double ch2 = 1.0-int0-int1+noise2;

      if(ch1<thresh1) ch1=0;
      chTot1+=ch1;
      xTot1+=ch1*(-3.0+strip);
      if(ch2<thresh2) ch2=0;
      chTot2+=ch2;
      xTot2+=ch2*(-3.0+strip);
    }
    double xmeasc1=xTot1/chTot1;
    histc.Fill(xmeasc1-xtrue);
    double xmeasc2=xTot2/chTot2;
    histd.Fill(xmeasc2-xtrue);

  }
  // Make a canvase and put the plot on it.

  TCanvas c1;
  c1.Divide(2,2);  // divide canvas into 4 squares
  c1.cd(1);
  hista.Draw();
  c1.cd(2);
  histb.Draw();
  c1.cd(3);
  histc.Draw();
  c1.cd(4);
  histd.Draw();

  c1.Print(fileName.c_str());

}           //End of macro
