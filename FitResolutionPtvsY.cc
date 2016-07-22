#include "iostream"
#include "iomanip"
#include "fstream"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TColor.h"
#include "TApplication.h"
#include "stdio.h"
#include "cstring"
#include "TMath.h"
#include <math.h>
#include <TLine.h>
#include <sys/stat.h>

using namespace std; 

#define pi 3.14159
#define set_title 0  			       

const int Ny = 7;

int black = 1;			   
int red = 2;				   
int blue = 4;				   
int fuchsia = 6;			   
int light_blue = 8;			   
int white = 10;			   
int color[4] = {blue,red,fuchsia,black}; 
int marker_style[4] = {24,25,26,20};     

double round(double x) {
  return floor(x * 10 + 0.5) / 10;
}

void SetStyle() {
  gStyle->SetErrorX(0.0);
  gStyle->SetOptStat(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetStatColor(kWhite);
  gStyle->SetPalette(1,0);  
  gStyle->SetFillColor(0);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleBorderSize(0);
  gROOT->ForceStyle();
}

void CreateDirectory(struct stat dir, TString dirname, bool wait) {
  using namespace std;
  if(stat(dirname, &dir) == 0 && S_ISDIR(dir.st_mode)) {
    cout<<endl;
    cout<<"directory "<<dirname<<" is already present!"<<endl;
  }
  else {
    if(mkdir(dirname,0755)==0) {
      cout<<endl;
      cout<<"directory "<<dirname<<" created successfully!"<<endl;
    }
  }
  
  cout<<"press enter to continue"<<endl;
  if(wait) getchar();
}

double Sigmoid(Double_t *x, Double_t *par){  
  double val = 1/(1 + par[0]*TMath::Exp(-par[1]*x[0]));
  return(val);
}

TF1* FitSigmoid(TH1F* histo, double min, double max) {

  TF1* f = new TF1("fSigmoid",Sigmoid,min,max,2);
  f->SetParNames("normalisation");

  f->SetLineColor(4);
  f->SetLineWidth(2);

  f->SetRange(min,max);

  for(int i = 0; i < 10; i++)
    histo->Fit("fSigmoid","RQNO");

  return f;
}

double Gaussian(Double_t *x, Double_t *par){

  //-- par[0] = normalisation
  //-- par[1] = mean
  //-- par[2] = sigma

  double value = par[0]*TMath::Exp(-TMath::Power(x[0]-par[1],2)/(2*par[2]*par[2]));
  return(value);
}


TF1* FitGaussian(TH1F* histo) {

  double parGa[3];
  parGa[0] = histo->GetBinContent(histo->GetMaximumBin());
  parGa[1] = histo->GetMean();
  parGa[2] = histo->GetRMS();
  
  double min = parGa[1] - 5 * parGa[2]; 
  double max = parGa[1] + 5 * parGa[2]; 
  
  TF1* g = new TF1("fGauss",Gaussian,min,max,3);
  g->SetParameters(parGa);
  g->SetParNames("normalisation","mu","sigma");
  
  for(int i = 0; i < 5; i++) {
    g->SetRange(min,max);
    histo->Fit("fGauss","RQNO");
    min = g->GetParameter(1) - 5 * g->GetParameter(2);
    max = g->GetParameter(1) + 5 * g->GetParameter(2);
  }

  g->SetLineColor(red);
  g->SetLineWidth(3);

  g->SetRange(2*min,2*max);
  
  cout<<endl<<g->GetName()<<": xi2/ndf = "<<g->GetChisquare()/g->GetNDF()<<", range = "<<min<<" - "<<max<<endl<<endl;
  
  return(g);
}

  
double Novosibirsk(Double_t *x, Double_t *par){  

  //-- par[0] = normalisation
  //-- par[1] = peak
  //-- par[2] = width
  //-- par[3] = tail

  if (TMath::Abs(par[3]) < 1.e-7) 
    return par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2],2));
  
  Double_t arg = 1.0 - (x[0]-par[1]) * par[3] / par[2];

  if (arg < 1.e-7) 
    return 0.0;
   
  Double_t log = TMath::Log(arg);    
  static const Double_t xi = 2.3548200450309494; //-- 2 Sqrt(Ln(4))
   
  Double_t width_zero = (2.0/xi)*TMath::ASinH(par[3]*xi*0.5);
  
  Double_t width_zero2 = width_zero * width_zero;
  Double_t exponent = -0.5*log*log/width_zero2 -0.5*width_zero2;
  Double_t value = par[0]*TMath::Exp(exponent);
  
  return value;
}

TF1* FitNovosibirsk(TH1F* histo, double min, double max, double init[4]) {
  
  TF1* f = new TF1("fNovosibirsk",Novosibirsk,min,max,4);
  f->SetParNames("normalisation","peak","width","tail");

  f->SetParameters(init);
  f->SetParLimits(2,0,100);

  f->SetLineColor(4);
  f->SetLineWidth(2);

  f->SetRange(min,max);

  for(int i = 0; i < 10; i++)
    histo->Fit("fNovosibirsk","RQNO");

  return f;
}

double EMG(Double_t *x, Double_t *par) {

  //-- par[0] = normalisation
  //-- par[1] = mean
  //-- par[2] = sigma
  //-- par[3] = lambda

  double exponent = (par[3]/2) * (2*par[1] + par[3]*par[2]*par[2] - 2*x[0]);
  double A = (par[3]/2)*TMath::Exp(exponent);

  double arg = (par[1] + par[3]*par[2]*par[2] - x[0]) / (TMath::Sqrt(2)*par[2]);
  double B = 1 - TMath::Erf(arg);
  
  double value = par[0]*A*B;

  return(value);
}

double DoubleEMG(Double_t *x, Double_t *par) {
  return(EMG(x,par)+EMG(x,&par[4]));
}
                  
TF1* FitEMG(TH1F* histo, double nmin, double nmax) {

  int Niter = 5;
  int Nrange = 20;

  //-- gaussian fit (nmin = 4 or min = -0.4, nmax = 5 or max = 1.2) 

  double parGa[3];
  parGa[0] = histo->GetBinContent(histo->GetMaximumBin());
  parGa[1] = histo->GetMean();
  parGa[2] = histo->GetRMS();
  
  double minGa = parGa[1] - 4 * parGa[2]; 
  double maxGa = parGa[1] + 4 * parGa[2]; 
  
  TF1* g = new TF1("fGauss",Gaussian,minGa,maxGa,3);
  g->SetParameters(parGa);
  
  for(int iter = 0; iter < Niter; iter++) {
    g->SetRange(minGa,maxGa);
    histo->Fit("fGauss","RQNO");
    minGa = g->GetParameter(1) - 4 * g->GetParameter(2);
    maxGa = g->GetParameter(1) + 4 * g->GetParameter(2);
  }
  
  double chi2Ga = g->GetChisquare()/g->GetNDF();

  //-- EMG fit

  double min0 = g->GetParameter(1) - nmin * g->GetParameter(2);
  double max0 = g->GetParameter(1) + nmax * g->GetParameter(2);
  
  TF1* f = new TF1("fEMG",EMG,min0,max0,4);
  f->SetParNames("normalisation","mu","sigma","lambda");
  
  double parEMG[4];
  parEMG[0] = g->GetParameter(0);
  parEMG[1] = g->GetParameter(1);
  parEMG[2] = g->GetParameter(2);
  parEMG[3] = 10;

  f->SetParLimits(0,0,100000000); //-- normalisation > 0
  f->SetParLimits(3,0,45);        //-- 0 < lambda < 45

  double chi2min = 100; 
  int lowbest = 0;             
  int upbest = 0;
  
  for(int low = 0; low < Nrange; low++) {

    double bmin = min0 - 0.25 * low * f->GetParameter(2);

    for(int up = 0; up < Nrange; up++) {
      
      double bmax = max0 + 0.25 * up * f->GetParameter(2);

      f->SetParameters(parEMG);
      f->SetRange(bmin,bmax);
    
      for(int iter = 0; iter < Niter; iter++) 
	histo->Fit("fEMG","RQNO");
    
    
      if(f->GetChisquare()/f->GetNDF() < chi2min) {
	chi2min = f->GetChisquare()/f->GetNDF();
	f->GetParameters(parEMG);
	lowbest = low;
	upbest = up;
      }

    }

  }
  
  //-- retrieve the best function

  double min = min0 - 0.25 * lowbest * f->GetParameter(2);
  double max = max0 + 0.25 * upbest * f->GetParameter(2);

  f->SetRange(min,max);
  f->SetParameters(parEMG);	

  for(int iter = 0; iter < Niter; iter++) 
    histo->Fit("fEMG","RQNO");

  f->SetLineColor(blue); 
  f->SetLineWidth(3);

  f->SetRange(2*min,2*max);

  cout<<endl<<f->GetName()<<": xi2/ndf = "<<f->GetChisquare()/f->GetNDF()<<", range = "<<min<<" - "<<max<<endl<<endl;

  return(f);
}

TF1* FitDoubleEMG(TH1F* histo, double nmin, double nmax) {

  int Niter = 5;
  int Nrange = 10;
  
  //-- gaussian fit (nmin = 4 or min = -0.4, nmax = 5 or max = 1.2) 

  double parGa[3];
  parGa[0] = histo->GetBinContent(histo->GetMaximumBin());
  parGa[1] = histo->GetMean();
  parGa[2] = histo->GetRMS();
  
  double minGa = parGa[1] - 4 * parGa[2]; 
  double maxGa = parGa[1] + 4 * parGa[2]; 
  
  TF1* g = new TF1("fGauss",Gaussian,minGa,maxGa,3);
  g->SetParameters(parGa);
  
  for(int iter = 0; iter < Niter; iter++) {
    g->SetRange(minGa,maxGa);
    histo->Fit("fGauss","RQNO");
    minGa = g->GetParameter(1) - 4 * g->GetParameter(2);
    maxGa = g->GetParameter(1) + 4 * g->GetParameter(2);
  }
  
  double chi2Ga = g->GetChisquare()/g->GetNDF();

  //-- EMG fit

  double min0 = g->GetParameter(1) - nmin * g->GetParameter(2);
  double max0 = g->GetParameter(1) + nmax * g->GetParameter(2);
  
  TF1* f1 = new TF1("fEMG1",EMG,min0,max0,4);
  TF1* f2 = new TF1("fEMG2",EMG,min0,max0,4);
  TF1* ftot = new TF1("fEMGtot",DoubleEMG,min0,max0,8);

  f1->SetParNames("normalisation","mu","sigma","lambda");
  f2->SetParNames("normalisation","mu","sigma","lambda");
  ftot->SetParNames("normalisation 1","mu 1","sigma 1","lambda 1","normalisation 2","mu 2","sigma 2","lambda 2");

  double parEMGtot[8];
  
  double parEMG1[4];
  parEMG1[0] = g->GetParameter(0);
  parEMG1[1] = g->GetParameter(1);
  parEMG1[2] = g->GetParameter(2);
  parEMG1[3] = 10;

  double parEMG2[4];			     
  parEMG2[0] = g->GetParameter(0);	     
  parEMG2[1] = g->GetParameter(1);	       
  parEMG2[2] = g->GetParameter(2);	     
  parEMG2[3] = 10;                           

  //-- normalisation > 0
  f1->SetParLimits(0,0,1000000);  
  f2->SetParLimits(0,0,1000000);  

  ftot->SetParLimits(0,0,1000000); 
  ftot->SetParLimits(4,0,1000000); 

  //-- sigma > 0
  f1->SetParLimits(2,0,1);   
  f2->SetParLimits(2,0,1);   

  ftot->SetParLimits(2,0,1); 
  ftot->SetParLimits(6,0,1); 

  //-- 0 < lambda < 30
  f1->SetParLimits(3,0,30); 
  f2->SetParLimits(3,0,30); 

  ftot->SetParLimits(3,0,30);
  ftot->SetParLimits(7,0,30);
  
  ftot->SetParameter(1,f1->GetParameter(1));
  ftot->SetParameter(2,f1->GetParameter(2));
  ftot->SetParameter(3,f1->GetParameter(3));
	
  ftot->SetParameter(5,f2->GetParameter(1));
  ftot->SetParameter(6,f2->GetParameter(2));
  ftot->SetParameter(7,f2->GetParameter(3));

  
  double chi2min = 100; 
  int lowbest = 0;             
  int upbest = 0;
  
  for(int low = 0; low < Nrange; low++) {

    double bmin = min0 - 0.25 * low * g->GetParameter(2);

    for(int up = 0; up < Nrange; up++) {

      double bmax = max0 + 0.25 * up * g->GetParameter(2);
    
      f1->SetParameters(parEMG1);
      f1->SetRange(bmin,g->GetParameter(1));
    
      f2->SetParameters(parEMG2);
      f2->SetRange(g->GetParameter(1),bmax);
      
      ftot->SetRange(bmin,bmax);
      
      for(int iter = 0; iter < Niter; iter++) {
	histo->Fit("fEMG1","RQNO");
	histo->Fit("fEMG2","RQNO+");
	
	ftot->SetParameter(1,f1->GetParameter(1));
	ftot->SetParameter(2,f1->GetParameter(2));
	ftot->SetParameter(3,f1->GetParameter(3));
	
	ftot->SetParameter(5,f2->GetParameter(1));
	ftot->SetParameter(6,f2->GetParameter(2));
	ftot->SetParameter(7,f2->GetParameter(3));
	
	histo->Fit("fEMGtot","RQNO+");
	
	if(ftot->GetChisquare()/ftot->GetNDF() < chi2min) {
	  chi2min = ftot->GetChisquare()/ftot->GetNDF();
	  f1->GetParameters(parEMG1);
	  f2->GetParameters(parEMG2);
	  ftot->GetParameters(parEMGtot);
	  lowbest = low;
	  upbest = up;
	}
      }
    }
  }
 
  //-- retrieve the best function
    
  double min = min0 - 0.25 * lowbest * g->GetParameter(2);
  double max = max0 + 0.25 * upbest  * g->GetParameter(2);  

  ftot->SetRange(min,max);
  ftot->SetParameters(parEMGtot);	
  
  for(int iter = 0; iter < Niter; iter++) 
    histo->Fit("fEMGtot","RQNO");
  
  ftot->SetLineColor(fuchsia);
  ftot->SetLineWidth(3);
  
  ftot->SetRange(2*min,2*max);

  cout<<endl<<ftot->GetName()<<": xi2/ndf = "<<ftot->GetChisquare()/ftot->GetNDF()<<", range = "<<min<<" - "<<max<<endl<<endl;

  return(ftot);
}
  
void SetGraphStyle(TGraphErrors* g,int marker_style,double marker_size,int marker_color,int line_color) {
  g->SetMarkerStyle(marker_style);
  g->SetMarkerSize(marker_size);
  g->SetMarkerColor(marker_color);
  g->SetLineColor(line_color);
}

void SetGraphStyle(TGraphAsymmErrors* g,int marker_style,double marker_size,int marker_color,int line_color) {
  g->SetMarkerStyle(marker_style);
  g->SetMarkerSize(marker_size);
  g->SetMarkerColor(marker_color);
  g->SetLineColor(line_color);
}

void SetLegendStyle(TLegend* legend, int color) {
  legend->SetTextSize(0.03);
  legend->SetFillColor(kWhite);
  legend->SetLineColor(kWhite);
  legend->SetTextColor(color);
}


void SetHistoStyleData(TH1F* histo,int color,int marker_style,double marker_size) {
  histo->SetLineColor(color);

  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(marker_style);
  histo->SetMarkerSize(marker_size);
}

void SetHistoStyleData(TH1F* histo,int marker_color,int line_color,int fill_color) {
  histo->SetMarkerColor(marker_color);
  histo->SetLineColor(line_color);
  histo->SetFillColor(fill_color);
}

void SetHistoStyleData(TH2D* histo,int color,int marker_style,double marker_size) {
  histo->SetLineColor(color);

  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(marker_style);
  histo->SetMarkerSize(marker_size);
}

void SetHistoStyleMC(TH1F* histo,int line_style,int line_color,double line_width) {
  histo->SetLineStyle(line_style);
  histo->SetLineColor(line_color);
  histo->SetLineWidth(line_width);
}

void SetHistoStyleMC(TH2F* histo,int line_style,int line_color,double line_width) {
  histo->SetLineStyle(line_style);
  histo->SetLineColor(line_color);
  histo->SetLineWidth(line_width);
}

void SetNullHisto(TH2D* histo, double x_offset, double y_offset, TString x_title, TString y_title) {

  histo->GetXaxis()->SetTitleOffset(x_offset);
  histo->GetYaxis()->SetTitleOffset(y_offset);
  histo->GetXaxis()->SetTitle(x_title);
  histo->GetYaxis()->SetTitle(y_title);
}

void SetHistoStyleMC(TH2D* histo, TString x_title, TString y_title, TString z_title, double title_offset, double title_size, double label_size) {

  histo->GetXaxis()->SetTitle(x_title);
  histo->GetXaxis()->SetTitleSize(title_size);
  histo->GetXaxis()->SetTitleOffset(title_offset);
  histo->GetXaxis()->SetLabelSize(label_size);
  
  histo->GetYaxis()->SetTitle(y_title);
  histo->GetYaxis()->SetTitleSize(title_size);
  histo->GetYaxis()->SetTitleOffset(1.25*title_offset);              
  histo->GetYaxis()->SetLabelSize(label_size);

  histo->GetZaxis()->SetTitle(z_title);  
  histo->GetZaxis()->SetTitleSize(title_size);
  histo->GetZaxis()->SetTitleOffset(title_offset);
  histo->GetZaxis()->SetLabelSize(label_size);
}

void ComputeProduct(double a, double b, double error_a, double error_b, double &product, double &error) {

  product = 0;
  error = 0;

  product = a*b;
  
  if(a*b != 0) {
    double error2 = TMath::Power(product,2)*(TMath::Power(error_a/a,2)+TMath::Power(error_b/b,2));
    error = TMath::Sqrt(error2);
  }

}

void ComputeHistoProduct(TH1D* ha, TH1D* hb, TH1D* hproduct) {

  for(int ibin = 1; ibin <= ha->GetNbinsX(); ibin++) {

    double a = ha->GetBinContent(ibin);
    double b = hb->GetBinContent(ibin);

    double error_a = ha->GetBinError(ibin);
    double error_b = hb->GetBinError(ibin);

    double product = 0;
    double error = 0;
    ComputeProduct(a,b,error_a,error_b,product,error);

    hproduct->SetBinContent(ibin,product);
    hproduct->SetBinError(ibin,error);
  }

}

void GetMinMaxHisto(TH1F* histo_down, TH1F* histo_up, TH1F* histo_min,TH1F* histo_max) {

  int Nbin = histo_down->GetNbinsX();

  double value[4];
  double value_min;
  double value_max;
  
  for(int ibin = 1; ibin <= Nbin; ibin++) {

    value[0] = histo_down->GetBinContent(ibin) + histo_down->GetBinError(ibin);
    value[1] = histo_down->GetBinContent(ibin) - histo_down->GetBinError(ibin);

    value[2] = histo_up->GetBinContent(ibin) + histo_up->GetBinError(ibin);
    value[3] = histo_up->GetBinContent(ibin) - histo_up->GetBinError(ibin);

    value_min = value[0];
    value_max = value[0];
    
    for(int i = 0; i < 4; i++) {
      if(value[i] > value_max) value_max = value[i];
      if(value[i] < value_min) value_min = value[i];
    }

    histo_min->SetBinContent(ibin,value_min);
    histo_max->SetBinContent(ibin,value_max);      
  }

} 

int main(int argc, char *argv[]) {

  struct stat dir;
  
  printf("argc = %d\n", argc);
  for (int i = 0; i<argc; i++){
    static int k = 0;
    printf("argv[%d] = %s\n", k, argv[k]);
    k++;
  }

  if(argc < 2) {
    cout<<endl;
    cout<<"usage: "<<"./Correlation-Plot \"wait\" "<<endl;
    cout<<"wait: true or false "<<endl;
    return(0);
  }

  TString decision = TString(argv[1]);
  bool wait = false;
  if (strcmp(decision,"true")==0) wait = true;
    
  int rebin = 1;
  double margin = 0.12;
  double canvas_size = 700;

  TH1F* histoResoPtvsY[Ny];
  TH1F* histoYmean[Ny];

  TH2D* hnull;
  
  TString title = "";
  TString canvas_name, canvas_title,canvas_label;
  TString legend_title;

  double ymax = 0;
  
  double marker_size = 0.6; //-- 0.8 
  double marker_size_medium = 0.6;
  int line_width = 2;

  int circle = 20;
  int opencircle = 24;
  
  TString histo_name;

  SetStyle();

  TString inputdir = "/Users/roland/DESY-DAS-School/Inclusive-Jet/InputFile/";

  TString inputfile = inputdir;
  inputfile+="OutputDAS-MC.root";

  TFile file(inputfile);
  
  file.cd("efficiency");
    
  TString outputdir = "/Users/roland/DESY-DAS-School/Inclusive-Jet/Result/ResolutionPtvsY/"; 
  TString save_label = "";

  CreateDirectory(dir,outputdir,true);
  TApplication* rootapp = new TApplication("app",&argc, argv);

  cout<<endl;

  //-- y bin

  TString ybin[Ny] = {" | y | < 0.5"," 0.5 < | y | < 1"," 1 < | y | < 1.5"," 1.5 < | y | < 2"," 2 < | y | < 2.5"," 2.5 < | y | < 3"," 3.2 < | y | < 4.7"};

  //-- fit

  TF1* fit_y_GA[Ny];
  TF1* fit_y_EMG[Ny];
  TF1* fit_y_DoubleEMG[Ny];

  double chi2_y_GA[Ny];
  double chi2_y_EMG[Ny];
  double chi2_y_DoubleEMG[Ny];
  
  //-- legend

  TLegend* legend_GA;
  TLegend* legend_EMG;
  TLegend* legend_DoubleEMG;
  
  TLegend* legend_mc;
  TLegend* legend_y;

  //----------------------------//
  //-- resolution pT versus y --//
  //----------------------------//

  TCanvas* canvas_sigpt_y[Ny];
  
  double sigmapt_y[Ny];
  double sigmapt_y_error[Ny];

  double ymean[Ny];
  double ymean_error[Ny];

  //-- loop over y bins
  for(int iy = 0; iy < Ny; iy++) {

    canvas_label = "canvas_sigpt_y_" + TString::Format("%d",iy+1);
    canvas_title = "pT resolution versus y - " + TString::Format("%d",iy+1);
    
    canvas_sigpt_y[iy] = new TCanvas(canvas_label,canvas_title,canvas_size,canvas_size);	
    canvas_sigpt_y[iy]->Divide(1,1);
    canvas_sigpt_y[iy]->SetLogy();
    ((TPad*) canvas_sigpt_y[iy])->cd();
    ((TPad*) canvas_sigpt_y[iy])->SetTicky();
    canvas_sigpt_y[iy]->SetLeftMargin(margin);  
    canvas_sigpt_y[iy]->SetRightMargin(margin); 
    
    histo_name = "efficiency/ResolutionPtvsY_" + TString::Format("%d",iy+1) + "bin";

    ymax = 0;

    histoResoPtvsY[iy] = (TH1F*) file.Get(histo_name);
    histoResoPtvsY[iy]->Rebin(2);
    SetHistoStyleData(histoResoPtvsY[iy],black,circle,marker_size_medium);
    
    if(histoResoPtvsY[iy]->GetBinContent(histoResoPtvsY[iy]->GetMaximumBin()) > ymax) ymax = histoResoPtvsY[iy]->GetBinContent(histoResoPtvsY[iy]->GetMaximumBin());

    histo_name = "efficiency/y_GENJetReso_" + TString::Format("%d",iy+1) + "bin";
    histoYmean[iy] = (TH1F*) file.Get(histo_name);
    ymean[iy] = histoYmean[iy]->GetMean();
    ymean_error[iy] = histoYmean[iy]->GetMeanError();;
	
    ymax*=100;
    
    hnull = new TH2D("hnull",title,100,-1,1,100,ymax/2000000,ymax);   
    SetNullHisto(hnull,1.3,1.4,"(p_{T} reco - p_{T} gen) /p_{T} gen","number of jets");
    hnull->DrawCopy();

    //-- Gaussian
    fit_y_GA[iy] = FitGaussian(histoResoPtvsY[iy]);
    fit_y_GA[iy]->Draw("same");
    chi2_y_GA[iy] = fit_y_GA[iy]->GetChisquare()/fit_y_GA[iy]->GetNDF();

    sigmapt_y[iy] = 100*fit_y_GA[iy]->GetParameter(2);
    sigmapt_y_error[iy] = 100*fit_y_GA[iy]->GetParError(2);
    
    for(int ipar = 0; ipar < 3; ipar++) 															     
      cout<<fit_y_GA[iy]->GetParName(ipar)<<" = "<<fit_y_GA[iy]->GetParameter(ipar)<<" +/- "<<fit_y_GA[iy]->GetParError(ipar)<<endl;
      
    //-- Exponentially modified gaussian     
    fit_y_EMG[iy] = FitEMG(histoResoPtvsY[iy],5,5);											  
    fit_y_EMG[iy]->Draw("same");															  
    chi2_y_EMG[iy] = fit_y_EMG[iy]->GetChisquare()/fit_y_EMG[iy]->GetNDF();
    
    for(int ipar = 0; ipar < 4; ipar++) 														  
      cout<<fit_y_EMG[iy]->GetParName(ipar)<<" = "<<fit_y_EMG[iy]->GetParameter(ipar)<<" +/- "<<fit_y_EMG[iy]->GetParError(ipar)<<endl;
      
    //-- Double exponentially modified gaussian     																  
    fit_y_DoubleEMG[iy] = FitDoubleEMG(histoResoPtvsY[iy],5,5);											  	  
    fit_y_DoubleEMG[iy]->Draw("same");
    chi2_y_DoubleEMG[iy] = fit_y_DoubleEMG[iy]->GetChisquare()/fit_y_DoubleEMG[iy]->GetNDF();

    for(int ipar = 0; ipar < 8; ipar++) {													   
      cout<<fit_y_DoubleEMG[iy]->GetParName(ipar)<<" = "<<fit_y_DoubleEMG[iy]->GetParameter(ipar)<<" +/- "<<fit_y_DoubleEMG[iy]->GetParError(ipar)<<endl;
      if(ipar == 3) cout<<endl;
    }
    
    histoResoPtvsY[iy]->Draw("E1same");
    
    legend_mc = new TLegend(0.17,0.82,0.22,0.87);	 
    SetLegendStyle(legend_mc,black);			 
    legend_mc->AddEntry("","PYTHIA8 CUETM1","");             
    legend_mc->Draw("same");                             
    
    legend_y = new TLegend(0.17,0.77,0.22,0.82);     
    SetLegendStyle(legend_y,black);		       
    legend_y->AddEntry("",ybin[iy],"");          
    legend_y->Draw("same");

    legend_GA = new TLegend(0.48,0.82,0.68,0.87);
    SetLegendStyle(legend_GA,red);
    legend_GA->AddEntry(fit_y_GA[iy]," Gaussian -  #chi^{2} / ndf = " + TString::Format("%2.1f",chi2_y_GA[iy]),"l");
    legend_GA->Draw("same");

    legend_EMG = new TLegend(0.48,0.77,0.68,0.82);		     
    SetLegendStyle(legend_EMG,blue);
    legend_EMG->AddEntry(fit_y_EMG[iy]," EMG -  #chi^{2} / ndf = " + TString::Format("%2.1f",chi2_y_EMG[iy]),"l");
    legend_EMG->Draw("same");

    legend_DoubleEMG = new TLegend(0.48,0.72,0.68,0.77);		    
    SetLegendStyle(legend_DoubleEMG,fuchsia);				    
    legend_DoubleEMG->AddEntry(fit_y_DoubleEMG[iy]," DEMG -  #chi^{2} / ndf = " + TString::Format("%2.1f",chi2_y_DoubleEMG[iy]),"l");
    legend_DoubleEMG->Draw("same");                                       

    ((TPad*) canvas_sigpt_y[iy])->RedrawAxis();

    cout<<endl;
    canvas_sigpt_y[iy]->Print(outputdir + "resolution-pt-vs-y-" + TString::Format("%d",iy+1) + ".pdf");
    if(iy == 0) canvas_sigpt_y[iy]->Print(outputdir + "resolution-pt-vs-y.pdf(","pdf");
    else canvas_sigpt_y[iy]->Print(outputdir + "resolution-pt-vs-y.pdf","pdf");

    if(wait) canvas_sigpt_y[iy]->WaitPrimitive();
    
    delete hnull;

    delete legend_mc;    
    delete legend_y;

    delete legend_GA;
    delete legend_EMG;
    delete legend_DoubleEMG;
  }

  //------------------------------------//
  //-- resolution pT versus y - curve --//  
  //------------------------------------//  

  TCanvas* canvas_sigpt_y_curve = new TCanvas("canvas_sigpt_y_curve","pT resolution versus y",canvas_size,canvas_size);	
  canvas_sigpt_y_curve->Divide(1,1);
  ((TPad*) canvas_sigpt_y_curve)->cd();
  ((TPad*) canvas_sigpt_y_curve)->SetTicky();
  canvas_sigpt_y_curve->SetLeftMargin(margin);  
  canvas_sigpt_y_curve->SetRightMargin(margin); 

  hnull = new TH2D("hnull",title,100,0,4,100,0,30);   
  SetNullHisto(hnull,1.3,1.4,"y","#sigma p_{T} / p_{T} [%]");
  hnull->DrawCopy();
  
  TGraphErrors* gsigpt_y = new TGraphErrors(Ny,ymean,sigmapt_y,ymean_error,sigmapt_y_error);
  SetGraphStyle(gsigpt_y,circle,marker_size,black,black);
  gsigpt_y->Draw("Psame");

  legend_mc = new TLegend(0.17,0.82,0.22,0.87);	 
  SetLegendStyle(legend_mc,black);			 
  legend_mc->AddEntry("","PYTHIA8 CUETM1","");           
  legend_mc->Draw("same");                             
  
  ((TPad*) canvas_sigpt_y_curve)->RedrawAxis();
  
  canvas_sigpt_y_curve->Print(outputdir + "resolution-pt-vs-y-curve.pdf");
  canvas_sigpt_y_curve->Print(outputdir + "resolution-pt-vs-y.pdf)","pdf");
  if(wait) canvas_sigpt_y_curve->WaitPrimitive();

  delete hnull;
  
  delete legend_mc;

  return(0);
}
