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

const int CutFit = 100;

const int Npt = 80;

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

  TH1F* histoResoPtvsPt[Npt];
  TH1F* histoPtmean[Npt];

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
  file.ls();
    
  TString outputdir = "/Users/roland/DESY-DAS-School/Inclusive-Jet/Result/ResolutionPtvsPt/"; 
  TString save_label = "";

  CreateDirectory(dir,outputdir,true);
  TApplication* rootapp = new TApplication("app",&argc, argv);

  cout<<endl;

  //-- pt bin

  TString ptbin[Npt+1] = {"0", "1", "5", "6", "8", "10", "12", "15", "18", "21", "24", "28", "32", "37", "43", "49", "56", "64", "74", "84","97", "114", "133", "153", "174", "196", "220", "245",
			  "272", "300", "330", "362", "395", "430", "468","507", "548", "592", "638", "686", "737", "790", "846", "905", "967","1032", "1101", "1172", "1248", "1327", "1410", "1497"
			  , "1588","1684", "1784", "1890", "2000","2116", "2238", "2366", "2500", "2640", "2787", "2941", "3103", "3273", "3450", "3637", "3832","4037", "4252", "4477", "4713", "4961",   
			  "5220", "5492", "5777", "6076", "6389", "6717", "7000"};                                                                                    

  //-- fit

  TF1* fit_pt_GA[Npt];
  TF1* fit_pt_EMG[Npt];
  TF1* fit_pt_DoubleEMG[Npt];

  double chi2_pt_GA[Npt];
  double chi2_pt_EMG[Npt];
  double chi2_pt_DoubleEMG[Npt];
  
  //-- legend

  TLegend* legend_GA;
  TLegend* legend_EMG;
  TLegend* legend_DoubleEMG;
  
  TLegend* legend_mc;
  TLegend* legend_pt;

  //-----------------------------//
  //-- resolution pT versus pT --//
  //-----------------------------//

  int ifilled = 0;
  
  TCanvas* canvas_sigpt_pt[Npt];
  
  double sigmapt_pt[Npt];
  double sigmapt_pt_error[Npt];

  double ptmean[Npt];
  double ptmean_error[Npt];       


  //-- loop over pt bins
  for(int ipt = 0; ipt < Npt; ipt++) {

    sigmapt_pt[ipt] = 0;	  
    sigmapt_pt_error[ipt] = 0;	  
    
    ptmean[ipt] = 0;		  
    ptmean_error[ipt] = 0;           
        
    canvas_label = "canvas_sigpt_pt_" + TString::Format("%d",ipt+1);			       
    canvas_title = "pT resolution versus pT - " + TString::Format("%d",ipt+1);		       
    											       
    canvas_sigpt_pt[ipt] = new TCanvas(canvas_label,canvas_title,canvas_size,canvas_size);     
    canvas_sigpt_pt[ipt]->Divide(1,1);							       
    canvas_sigpt_pt[ipt]->SetLogy();							       
    ((TPad*) canvas_sigpt_pt[ipt])->cd();						       
    ((TPad*) canvas_sigpt_pt[ipt])->SetTicky();						       
    canvas_sigpt_pt[ipt]->SetLeftMargin(margin);  					       
    canvas_sigpt_pt[ipt]->SetRightMargin(margin);                                              
    
    histo_name = "efficiency/ResolutionPtvsPt_" + TString::Format("%d",ipt+1) + "bin";
    
    ymax = 0;

    histoResoPtvsPt[ipt] = (TH1F*) file.Get(histo_name);    
    histoResoPtvsPt[ipt]->Rebin(2);
    SetHistoStyleData(histoResoPtvsPt[ipt],black,circle,marker_size_medium);
    
    if(histoResoPtvsPt[ipt]->GetBinContent(histoResoPtvsPt[ipt]->GetMaximumBin()) > ymax) ymax = histoResoPtvsPt[ipt]->GetBinContent(histoResoPtvsPt[ipt]->GetMaximumBin());

    bool IsFilled = histoResoPtvsPt[ipt]->Integral() > 0.01;
    
    histo_name = "efficiency/pt_GENJetReso_" + TString::Format("%d",ipt+1) + "bin";
    histoPtmean[ipt] = (TH1F*) file.Get(histo_name);
    if(histoPtmean[ipt]->Integral() > 0) ptmean[ipt] = histoPtmean[ipt]->GetMean();
    ptmean_error[ipt] = histoPtmean[ipt]->GetMeanError();;
    
    ymax*=100;

    hnull = new TH2D("hnull",title,100,-1,1,100,ymax/2000000,ymax);   
    SetNullHisto(hnull,1.3,1.4,"(p_{T} reco - p_{T} gen) /p_{T} gen","number of jets");
    hnull->DrawCopy();

    if(!IsFilled) continue;
      
    //-- Gaussian
    fit_pt_GA[ipt] = FitGaussian(histoResoPtvsPt[ipt]);
    chi2_pt_GA[ipt] = fit_pt_GA[ipt]->GetChisquare()/fit_pt_GA[ipt]->GetNDF();
    
    if(chi2_pt_GA[ipt] < CutFit) fit_pt_GA[ipt]->Draw("same");

    sigmapt_pt[ipt] = 100*fit_pt_GA[ipt]->GetParameter(2);
    sigmapt_pt_error[ipt] = 100*fit_pt_GA[ipt]->GetParError(2);
    
    for(int ipar = 0; ipar < 3; ipar++) 															     
      cout<<fit_pt_GA[ipt]->GetParName(ipar)<<" = "<<fit_pt_GA[ipt]->GetParameter(ipar)<<" +/- "<<fit_pt_GA[ipt]->GetParError(ipar)<<endl;
    
    //-- Exponentially modified gaussian     
    fit_pt_EMG[ipt] = FitEMG(histoResoPtvsPt[ipt],5,5);
    chi2_pt_EMG[ipt] = fit_pt_EMG[ipt]->GetChisquare()/fit_pt_EMG[ipt]->GetNDF();

    if(chi2_pt_EMG[ipt] > CutFit) fit_pt_EMG[ipt] = FitEMG(histoResoPtvsPt[ipt],4.5,4.5);
    chi2_pt_EMG[ipt] = fit_pt_EMG[ipt]->GetChisquare()/fit_pt_EMG[ipt]->GetNDF();

    if(chi2_pt_EMG[ipt] > CutFit) fit_pt_EMG[ipt] = FitEMG(histoResoPtvsPt[ipt],4.0,4.0);
    chi2_pt_EMG[ipt] = fit_pt_EMG[ipt]->GetChisquare()/fit_pt_EMG[ipt]->GetNDF();
      
    if(chi2_pt_EMG[ipt] < CutFit) fit_pt_EMG[ipt]->Draw("same");
      
    for(int ipar = 0; ipar < 4; ipar++) 														  
      cout<<fit_pt_EMG[ipt]->GetParName(ipar)<<" = "<<fit_pt_EMG[ipt]->GetParameter(ipar)<<" +/- "<<fit_pt_EMG[ipt]->GetParError(ipar)<<endl;
    
    //-- Double exponentially modified gaussian     																  
    fit_pt_DoubleEMG[ipt] = FitDoubleEMG(histoResoPtvsPt[ipt],5,5);											  	  
    chi2_pt_DoubleEMG[ipt] = fit_pt_DoubleEMG[ipt]->GetChisquare()/fit_pt_DoubleEMG[ipt]->GetNDF();

    if(ipt == 39) fit_pt_DoubleEMG[ipt] = FitDoubleEMG(histoResoPtvsPt[ipt],4.5,4.5);				   
    chi2_pt_DoubleEMG[ipt] = fit_pt_DoubleEMG[ipt]->GetChisquare()/fit_pt_DoubleEMG[ipt]->GetNDF();
      
    if(chi2_pt_DoubleEMG[ipt] < CutFit) fit_pt_DoubleEMG[ipt]->Draw("same");
    
    for(int ipar = 0; ipar < 8; ipar++) {													   
      cout<<fit_pt_DoubleEMG[ipt]->GetParName(ipar)<<" = "<<fit_pt_DoubleEMG[ipt]->GetParameter(ipar)<<" +/- "<<fit_pt_DoubleEMG[ipt]->GetParError(ipar)<<endl;
      if(ipar == 3) cout<<endl;
    }

    histoResoPtvsPt[ipt]->Draw("E1same");
    
    legend_mc = new TLegend(0.17,0.82,0.22,0.87);	 
    SetLegendStyle(legend_mc,black);			 
    legend_mc->AddEntry("","PYTHIA8 CUETM1","");             
    legend_mc->Draw("same");                             
    
    legend_pt = new TLegend(0.17,0.77,0.22,0.82);     
    SetLegendStyle(legend_pt,black);		       
    legend_pt->AddEntry("",ptbin[ipt] + " <  p_{ T} < " + ptbin[ipt+1] + " GeV","");          
    legend_pt->Draw("same");
    
    TString chi2_label_GA = "";
    if(chi2_pt_GA[ipt] < CutFit) chi2_label_GA = " Gaussian -  #chi^{2} / ndf = " + TString::Format("%2.1f",chi2_pt_GA[ipt]);
    else chi2_label_GA = " Gaussian -  #chi^{2} / ndf > " + TString::Format("%d",CutFit);

    TString chi2_label_EMG = "";
    if(chi2_pt_EMG[ipt] < CutFit) chi2_label_EMG = " EMG -  #chi^{2} / ndf = " + TString::Format("%2.1f",chi2_pt_EMG[ipt]);
    else chi2_label_EMG = " EMG -  #chi^{2} / ndf > " + TString::Format("%d",CutFit);

    TString chi2_label_DoubleEMG = "";
    if(chi2_pt_DoubleEMG[ipt] < CutFit) chi2_label_DoubleEMG = " DEMG -  #chi^{2} / ndf = " + TString::Format("%2.1f",chi2_pt_DoubleEMG[ipt]);
    else chi2_label_DoubleEMG = " DEMG -  #chi^{2} / ndf > " + TString::Format("%d",CutFit);
      
    legend_GA = new TLegend(0.48,0.82,0.68,0.87);
    SetLegendStyle(legend_GA,red);
    legend_GA->AddEntry(fit_pt_GA[ipt],chi2_label_GA,"l");
    legend_GA->Draw("same");
    
    legend_EMG = new TLegend(0.48,0.77,0.68,0.82);		     
    SetLegendStyle(legend_EMG,blue);
    legend_EMG->AddEntry(fit_pt_EMG[ipt],chi2_label_EMG,"l");
    legend_EMG->Draw("same");
    
    legend_DoubleEMG = new TLegend(0.48,0.72,0.68,0.77);		    
    SetLegendStyle(legend_DoubleEMG,fuchsia);				    
    legend_DoubleEMG->AddEntry(fit_pt_DoubleEMG[ipt],chi2_label_DoubleEMG,"l");
    legend_DoubleEMG->Draw("same");                                       
  
  
    ((TPad*) canvas_sigpt_pt[ipt])->RedrawAxis();
    
    cout<<endl;

    canvas_sigpt_pt[ipt]->Print(outputdir + "resolution-pt-vs-pt-" + TString::Format("%d",ifilled+1) + ".pdf");
    if(ifilled == 0) canvas_sigpt_pt[ipt]->Print(outputdir + "resolution-pt-vs-pt.pdf(","pdf");
    else canvas_sigpt_pt[ipt]->Print(outputdir + "resolution-pt-vs-pt.pdf","pdf");
    
    if(wait) canvas_sigpt_pt[ipt]->WaitPrimitive();

    delete hnull;
    
    delete legend_mc;    
    delete legend_pt;
    
    delete legend_GA;
    delete legend_EMG;
    delete legend_DoubleEMG;
    ifilled++;    
  }

  //-------------------------------------//
  //-- resolution pT versus pT - curve --//  
  //-------------------------------------//  

  for(int ipt = 0; ipt < Npt; ipt++) {
    cout<<ptbin[ipt]<<" <  pT < "<<ptbin[ipt+1]<<" GeV: pT mean = "<<ptmean[ipt]<<" GeV - pT resolution "<<sigmapt_pt[ipt]<<" +/- "<<sigmapt_pt_error[ipt]<<" %"<<endl;
  }
  
  TCanvas* canvas_sigpt_pt_curve = new TCanvas("canvas_sigpt_pt_curve","pT resolution versus y",canvas_size,canvas_size);	
  canvas_sigpt_pt_curve->Divide(1,1);
  ((TPad*) canvas_sigpt_pt_curve)->cd();
  ((TPad*) canvas_sigpt_pt_curve)->SetTicky();
  canvas_sigpt_pt_curve->SetLeftMargin(margin);  
  canvas_sigpt_pt_curve->SetRightMargin(margin); 

  hnull = new TH2D("hnull",title,100,80,1500,100,0,16);   
  SetNullHisto(hnull,1.3,1.4,"p_{T} [GeV]","#sigma p_{T} / p_{T} [%]");
  hnull->DrawCopy();
  
  TGraphErrors* gsigpt_pt = new TGraphErrors(Npt,ptmean,sigmapt_pt,ptmean_error,sigmapt_pt_error);
  SetGraphStyle(gsigpt_pt,circle,marker_size,black,black);
  gsigpt_pt->Draw("Psame");

  legend_mc = new TLegend(0.17,0.82,0.22,0.87);	 
  SetLegendStyle(legend_mc,black);			 
  legend_mc->AddEntry("","PYTHIA8 CUETM1","");           
  legend_mc->Draw("same");                             
  
  ((TPad*) canvas_sigpt_pt_curve)->RedrawAxis();
  
  canvas_sigpt_pt_curve->Print(outputdir + "resolution-pt-vs-pt-curve.pdf");
  canvas_sigpt_pt_curve->Print(outputdir + "resolution-pt-vs-pt.pdf)","pdf");
  if(wait) canvas_sigpt_pt_curve->WaitPrimitive();

  delete hnull;
  
  delete legend_mc;

  return(0);
}
