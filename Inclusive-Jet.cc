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

const int Nfile = 2;	       
const int Ntrigger = 9;

const int Ny = 7;
const int Npt = 80;

const int CutFit = 100;

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

void NormalizeTH2Column(TH2F* histo) {

  int Nbinx = histo->GetNbinsX();
  int Nbiny = histo->GetNbinsY();

  double column[Nbinx];
  
  for(int ibinx = 1; ibinx <= Nbinx; ibinx++) {
    column[ibinx-1] = 0;
    for(int ibiny = 1; ibiny <= Nbiny; ibiny++) {
      column[ibinx-1]+=histo->GetBinContent(ibinx,ibiny);
    }
  }

  for(int ibinx = 1; ibinx <= Nbinx; ibinx++) {
    vector<double> result_temp;
    
    for(int ibiny = 1; ibiny <= Nbiny; ibiny++) {
      double content = 0;
      if(column[ibinx-1] > 0) content = 100*histo->GetBinContent(ibinx,ibiny)/column[ibinx-1];

      histo->SetBinContent(ibinx,ibiny,content);
    }
  }

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

  TLine* line;
  
  TH1F* hR[2];
  
  TH1F* heffi[Ntrigger];
  
  TH1F* histo[Nfile];
  double scaling[Nfile];
  
  TH2F* histoPUVX;
  TH2F* histoResponse[Ny];
  
  TH1F* histoPt[Nfile][Ny];
  double scalingPt[Nfile][Ny];  

  TH1F* histoPtDOWN[Nfile][Ny];
  double scalingPtDOWN[Nfile][Ny];
  
  TH1F* histoPtUP[Nfile][Ny];  
  double scalingPtUP[Nfile][Ny];
  
  TH1F* histoPtMIN[Nfile][Ny];
  TH1F* histoPtMAX[Nfile][Ny];

  TH1F* histoPtUncor[Nfile][Ny];
  double scalingPtUncor[Nfile][Ny];
  
  TH1F* histoPtHadron[Nfile][Ny];

  TH1F* histoResoPtvsY[Ny];	   
  TH1F* histoYmean[Ny];			   
					   
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
  
  int black = 1;
  int blue = 4;
  int red = 2;
  int fuchsia = 6;
  int light_blue = 8;
  int white = 10;
  int color[4] = {blue,red,fuchsia,black};
  int marker_style[4] = {24,25,26,20};
  
  TString histo_name;

  SetStyle();

  TString inputdir = "/Users/roland/DESY-DAS-School/Inclusive-Jet/InputFile/";

  TString inputfile_mc = inputdir;
  TString inputfile_data =  inputdir;

  TString inputfile[Nfile];
  for(int i = 0; i < Nfile; i++) inputfile[i] = inputdir;

  inputfile_mc+="OutputDAS-MC.root";
  inputfile_data+="OutputDAS-Data.root";
 
  inputfile[0]+="OutputDAS-MC.root";
  inputfile[1]+="OutputDAS-Data.root";
   
  TFile file_mc(inputfile_mc);
  TFile file_data(inputfile_data);
  
  file_mc.cd("efficiency");
  file_data.cd("efficiency");
    
  TFile* file[Nfile];
  for(int ifile = 0; ifile < Nfile; ifile++) {
    file[ifile] = TFile::Open(inputfile[ifile]);
    file[ifile]->cd("efficiency");
    file[ifile]->ls();
  }
  
  TString outputdir = "/Users/roland/DESY-DAS-School/Inclusive-Jet/Result/"; 
  TString save_label = "";

  CreateDirectory(dir,outputdir,true);
  TApplication* rootapp = new TApplication("app",&argc, argv);

  cout<<endl;

  //-- y bin

  double ptmin = 50;
  
  TString ybin[Ny] = {" | y | < 0.5"," 0.5 < | y | < 1"," 1 < | y | < 1.5"," 1.5 < | y | < 2"," 2 < | y | < 2.5"," 2.5 < | y | < 3"," 3.2 < | y | < 4.7"};
  TString sample[Nfile] = {"pythia8","data"};

  //-- pt bin

  TString ptbin[Npt+1] = {"0", "1", "5", "6", "8", "10", "12", "15", "18", "21", "24", "28", "32", "37", "43", "49", "56", "64", "74", "84","97", "114", "133", "153", "174", "196", "220", "245",
			  "272", "300", "330", "362", "395", "430", "468","507", "548", "592", "638", "686", "737", "790", "846", "905", "967","1032", "1101", "1172", "1248", "1327", "1410", "1497"
			  , "1588","1684", "1784", "1890", "2000","2116", "2238", "2366", "2500", "2640", "2787", "2941", "3103", "3273", "3450", "3637", "3832","4037", "4252", "4477", "4713", "4961",   
			  "5220", "5492", "5777", "6076", "6389", "6717", "7000"};                                                                                    

  //-- fit

  TF1* fit_y_GA[Ny];
  double chi2_y_GA[Ny];
  
  TF1* fit_pt_GA[Npt];  
  double chi2_pt_GA[Npt];         
  
  //-- legend

  TLegend* legend_fit1;
  TLegend* legend_fit2;
  
  TLegend* legend_level;	
  TLegend* legend_JES;
  TLegend* legend_y;
  TLegend* legend_pt;
  TLegend* legend_trigger[Ntrigger];
  
  TLegend* legend_mc;
  TLegend* legend_mc_uncor;		
  TLegend* legend_mc_hadron;
  
  TLegend* legend_data;
  TLegend* legend_data_uncor;

  //---------------------------------------//
  //-- detector level - (un)corrected pT --//
  //---------------------------------------//

  TCanvas* canvas_pt_det[Ny];
   
  //-- loop over y bins
  for(int iy = 0; iy < Ny; iy++) {

    canvas_label = "canvas_pt_det_" + TString::Format("%d",iy+1);
    canvas_title = "pT detector level - " + TString::Format("%d",iy+1);
    
    canvas_pt_det[iy] = new TCanvas(canvas_label,canvas_title,canvas_size,canvas_size);	
    canvas_pt_det[iy]->Divide(1,1);
    canvas_pt_det[iy]->SetLogx();
    canvas_pt_det[iy]->SetLogy();
    ((TPad*) canvas_pt_det[iy])->cd();
    ((TPad*) canvas_pt_det[iy])->SetTicky();
    canvas_pt_det[iy]->SetLeftMargin(margin);                                             
    canvas_pt_det[iy]->SetRightMargin(margin); 
    
    //-- JES nominal
    
    histo_name = "efficiency/pt_DETJet_" + TString::Format("%d",iy+1) + "bin";

    ymax = 0;

    //-- scaling
    for(int ifile = 0; ifile < Nfile; ifile++) 
      histoPt[ifile][iy] = (TH1F*) file[ifile]->Get(histo_name);
    
    for(int ifile = 0; ifile < Nfile; ifile++) 
      scalingPt[ifile][iy] = histoPt[1][iy]->Integral()/histoPt[ifile][iy]->Integral();
    
    //-- loop on data and monte carlo
    for(int ifile = 0; ifile < Nfile; ifile++) {
      histoPt[ifile][iy] = (TH1F*) file[ifile]->Get(histo_name);
      histoPt[ifile][iy]->Scale(scalingPt[ifile][iy]);
      if(ifile == 1) SetHistoStyleData(histoPt[ifile][iy],black,circle,marker_size_medium);
      else SetHistoStyleMC(histoPt[ifile][iy],1,blue,line_width);

      if(histoPt[ifile][iy]->GetBinContent(histoPt[ifile][iy]->GetMaximumBin()) > ymax) ymax = histoPt[ifile][iy]->GetBinContent(histoPt[ifile][iy]->GetMaximumBin());
    }

    //-- JES UP
    
    histo_name = "efficiency/pt_DETJetUP_" + TString::Format("%d",iy+1) + "bin";


    //-- scaling
    for(int ifile = 0; ifile < Nfile; ifile++) 
      histoPtUP[ifile][iy] = (TH1F*) file[ifile]->Get(histo_name);
    
    for(int ifile = 0; ifile < Nfile; ifile++) 
      scalingPtUP[ifile][iy] = histoPtUP[1][iy]->Integral()/histoPtUP[ifile][iy]->Integral();
    
    //-- loop on data and monte carlo
    for(int ifile = 0; ifile < Nfile; ifile++) {
      histoPtUP[ifile][iy] = (TH1F*) file[ifile]->Get(histo_name);
      histoPtUP[ifile][iy]->Scale(scalingPtUP[ifile][iy]);
      if(ifile == 1) SetHistoStyleData(histoPtUP[ifile][iy],black,circle,marker_size_medium);
      else SetHistoStyleMC(histoPtUP[ifile][iy],1,blue,line_width);

      if(histoPtUP[ifile][iy]->GetBinContent(histoPtUP[ifile][iy]->GetMaximumBin()) > ymax) ymax = histoPtUP[ifile][iy]->GetBinContent(histoPtUP[ifile][iy]->GetMaximumBin());
    }

    //-- JES DOWN
    
    histo_name = "efficiency/pt_DETJetDOWN_" + TString::Format("%d",iy+1) + "bin";


    //-- scaling
    for(int ifile = 0; ifile < Nfile; ifile++) 
      histoPtDOWN[ifile][iy] = (TH1F*) file[ifile]->Get(histo_name);
    
    for(int ifile = 0; ifile < Nfile; ifile++) 
      scalingPtDOWN[ifile][iy] = histoPtDOWN[1][iy]->Integral()/histoPtDOWN[ifile][iy]->Integral();
    
    //-- loop on data and monte carlo
    for(int ifile = 0; ifile < Nfile; ifile++) {
      histoPtDOWN[ifile][iy] = (TH1F*) file[ifile]->Get(histo_name);
      histoPtDOWN[ifile][iy]->Scale(scalingPtDOWN[ifile][iy]);
      if(ifile == 1) SetHistoStyleData(histoPtDOWN[ifile][iy],black,circle,marker_size_medium);
      else SetHistoStyleMC(histoPtDOWN[ifile][iy],1,blue,line_width);

      if(histoPtDOWN[ifile][iy]->GetBinContent(histoPtDOWN[ifile][iy]->GetMaximumBin()) > ymax) ymax = histoPtDOWN[ifile][iy]->GetBinContent(histoPtDOWN[ifile][iy]->GetMaximumBin());
    }

    //-- UNCORRECTED
    
    histo_name = "efficiency/pt_DETJetUncor_" + TString::Format("%d",iy+1) + "bin";

    ymax = 0;

    //-- scaling
    for(int ifile = 0; ifile < Nfile; ifile++) 
      histoPtUncor[ifile][iy] = (TH1F*) file[ifile]->Get(histo_name);
    
    for(int ifile = 0; ifile < Nfile; ifile++) 
      scalingPtUncor[ifile][iy] = histoPtUncor[1][iy]->Integral()/histoPtUncor[ifile][iy]->Integral();
    
    //-- loop on data and monte carlo
    for(int ifile = 0; ifile < Nfile; ifile++) {
      histoPtUncor[ifile][iy] = (TH1F*) file[ifile]->Get(histo_name);
      histoPtUncor[ifile][iy]->Scale(scalingPtUncor[ifile][iy]);
      if(ifile == 1) SetHistoStyleData(histoPtUncor[ifile][iy],black,opencircle,marker_size_medium);
      else SetHistoStyleMC(histoPtUncor[ifile][iy],1,red,line_width);

      if(histoPtUncor[ifile][iy]->GetBinContent(histoPtUncor[ifile][iy]->GetMaximumBin()) > ymax) ymax = histoPtUncor[ifile][iy]->GetBinContent(histoPtUncor[ifile][iy]->GetMaximumBin());
    }

    ymax*=100000;
    double ptmax = 0;
    if(iy == 0 || iy == 1 || iy == 2) ptmax = 2000;
    if(iy == 3) ptmax = 1500;
    if(iy == 4) ptmax = 1000;
    if(iy == 5) ptmax = 800;
    if(iy == 6) ptmax = 500;
    
    hnull = new TH2D("hnull",title,100,ptmin,ptmax,100,1,ymax);   
    SetNullHisto(hnull,1.3,1.4,"p_{T} [GeV]","number of jets");
    hnull->DrawCopy();

    for(int ifile = 0; ifile < Nfile; ifile++) {
      histo_name = "histoPtMIN_" + TString::Format("%d",ifile) + "_" + TString::Format("%d",iy);      
      histoPtMIN[ifile][iy] = (TH1F*) histoPt[ifile][iy]->Clone(histo_name);

      histo_name = "histoPtMAX_" + TString::Format("%d",ifile) + "_" + TString::Format("%d",iy);  
      histoPtMAX[ifile][iy] = (TH1F*) histoPt[ifile][iy]->Clone(histo_name);                      
    
      GetMinMaxHisto(histoPtDOWN[ifile][iy],histoPtUP[ifile][iy],histoPtMIN[ifile][iy],histoPtMAX[ifile][iy]);
    }

    histoPtMAX[1][iy]->SetFillColor(7);
    histoPtMAX[1][iy]->SetLineColor(7);
    histoPtMAX[1][iy]->Draw("][histosame");

    histoPtMIN[1][iy]->SetFillColor(10);    
    histoPtMIN[1][iy]->SetLineColor(10);        
    histoPtMIN[1][iy]->Draw("][histosame");

    for(int ifile = 0; ifile < Nfile; ifile++) {
      if(ifile == 1) histoPt[ifile][iy]->Draw("E1same");
      else histoPt[ifile][iy]->Draw("histosame][");
      
      if(ifile == 1) histoPtUncor[ifile][iy]->Draw("E1same");
      else histoPtUncor[ifile][iy]->Draw("histosame][");           
    }

    legend_level = new TLegend(0.17,0.82,0.22,0.87);	 
    SetLegendStyle(legend_level,black);			 
    legend_level->AddEntry("","detector level","");             
    legend_level->Draw("same");                              
    
    legend_data = new TLegend(0.38,0.82,0.63,0.87);
    SetLegendStyle(legend_data,black);
    legend_data->AddEntry(histoPt[1][iy]," data p_{ T} corrected","lp");     
    legend_data->Draw("same");

    legend_data_uncor = new TLegend(0.38,0.77,0.63,0.82);
    SetLegendStyle(legend_data_uncor,black);
    legend_data_uncor->AddEntry(histoPtUncor[1][iy]," data p_{ T} raw","lp");     
    legend_data_uncor->Draw("same");

    legend_mc = new TLegend(0.38,0.72,0.63,0.77);					
    SetLegendStyle(legend_mc,blue);							
    legend_mc->AddEntry(histoPt[0][iy]," PYTHIA8 CUETM1 p_{   T} corrected","lp");        
    legend_mc->Draw("same");								
                                                                                    
    legend_mc_uncor = new TLegend(0.38,0.67,0.63,0.72);					
    SetLegendStyle(legend_mc_uncor,red);							
    legend_mc_uncor->AddEntry(histoPtUncor[0][iy]," PYTHIA8 CUETM1 p_{   T} raw","lp");    
    legend_mc_uncor->Draw("same");
		       
    legend_JES = new TLegend(0.38,0.62,0.63,0.67);		       
    SetLegendStyle(legend_JES,black);
    legend_JES->AddEntry(histoPtMAX[1][iy]," JES uncertainty","f");
    legend_JES->Draw("same");
    
    legend_y = new TLegend(0.38,0.57,0.63,0.62);     
    SetLegendStyle(legend_y,black);		       
    legend_y->AddEntry("",ybin[iy],"");          
    legend_y->Draw("same");                          

    ((TPad*) canvas_pt_det[iy])->RedrawAxis();
    
    canvas_pt_det[iy]->Print(outputdir + "pt-det-" + TString::Format("%d",iy+1) + ".pdf");
    if(iy == 0) canvas_pt_det[iy]->Print(outputdir + "inclusive-jet.pdf(","pdf");
    else canvas_pt_det[iy]->Print(outputdir + "inclusive-jet.pdf","pdf");
    if(wait) canvas_pt_det[iy]->WaitPrimitive();

    delete hnull;

    delete legend_level;
    
    delete legend_mc;
    delete legend_data;

    delete legend_JES;
    delete legend_y;
  }

  //-----------------------//
  //-- hadron level - pT --//
  //-----------------------//

  TCanvas* canvas_pt_had[Ny];
   
  //-- loop over y bins
  for(int iy = 0; iy < Ny; iy++) {

    canvas_label = "canvas_pt_had_" + TString::Format("%d",iy+1);
    canvas_title = "pT hadron level - " + TString::Format("%d",iy+1);
    
    canvas_pt_had[iy] = new TCanvas(canvas_label,canvas_title,canvas_size,canvas_size);	
    canvas_pt_had[iy]->Divide(1,1);
    canvas_pt_had[iy]->SetLogx();
    canvas_pt_had[iy]->SetLogy();
    ((TPad*) canvas_pt_had[iy])->cd();
    ((TPad*) canvas_pt_had[iy])->SetTicky();
    canvas_pt_had[iy]->SetLeftMargin(margin);  
    canvas_pt_had[iy]->SetRightMargin(margin); 
    
    //-- GENERATED
    
    histo_name = "efficiency/pt_GENJet_" + TString::Format("%d",iy+1) + "bin";

    ymax = 0;

    //-- loop on monte carlo
    for(int ifile = 0; ifile < Nfile; ifile++) {
      if(ifile == 1) continue;
      histoPtHadron[ifile][iy] = (TH1F*) file[ifile]->Get(histo_name);
      SetHistoStyleData(histoPtHadron[ifile][iy],black,circle,marker_size_medium);

      if(histoPtHadron[ifile][iy]->GetBinContent(histoPtHadron[ifile][iy]->GetMaximumBin()) > ymax) ymax = histoPtHadron[ifile][iy]->GetBinContent(histoPtHadron[ifile][iy]->GetMaximumBin());
    }

    ymax*=100;
    double ptmax = 0;
    if(iy == 0 || iy == 1 || iy == 2) ptmax = 2000;
    if(iy == 3) ptmax = 1500;
    if(iy == 4) ptmax = 1000;
    if(iy == 5) ptmax = 800;
    if(iy == 6) ptmax = 500;
    
    hnull = new TH2D("hnull",title,100,ptmin,ptmax,100,0.01,ymax);   
    SetNullHisto(hnull,1.3,1.4,"p_{T} [GeV]","number of jets");
    hnull->DrawCopy();

    for(int ifile = 0; ifile < Nfile; ifile++) {
      if(ifile == 1) continue;
      histoPtHadron[ifile][iy]->Draw("E1same");

      histoPt[ifile][iy]->Scale(1/scalingPt[ifile][iy]);
      histoPt[ifile][iy]->Draw("histosame][");

      histoPtUncor[ifile][iy]->Scale(1/scalingPtUncor[ifile][iy]);
      histoPtUncor[ifile][iy]->Draw("histosame][");
    }

    legend_level = new TLegend(0.17,0.82,0.22,0.87);	 
    SetLegendStyle(legend_level,black);			 
    legend_level->AddEntry("","PYTHIA8 CUETM1","");             
    legend_level->Draw("same");                             
    
    legend_mc_hadron = new TLegend(0.48,0.82,0.68,0.87);
    SetLegendStyle(legend_mc_hadron,black);
    legend_mc_hadron->AddEntry(histoPtHadron[0][iy]," p_{T} hadron level","lp");     
    legend_mc_hadron->Draw("same");

    legend_mc = new TLegend(0.48,0.77,0.68,0.82);			       
    SetLegendStyle(legend_mc,blue);					       
    legend_mc->AddEntry(histoPt[0][iy]," p_{T} detector level corrected","lp");   
    legend_mc->Draw("same");

    legend_mc_uncor = new TLegend(0.48,0.72,0.68,0.77);			       	     
    SetLegendStyle(legend_mc_uncor,red);					       	     
    legend_mc_uncor->AddEntry(histoPtUncor[0][iy]," p_{T} detector level raw","lp"); 
    legend_mc_uncor->Draw("same");                                                           
    
    legend_y = new TLegend(0.48,0.67,0.68,0.72);     
    SetLegendStyle(legend_y,black);		       
    legend_y->AddEntry("",ybin[iy],"");          
    legend_y->Draw("same");                          

    ((TPad*) canvas_pt_had[iy])->RedrawAxis();
    
    canvas_pt_had[iy]->Print(outputdir + "pt-hadron-" + TString::Format("%d",iy+1) + ".pdf");
    canvas_pt_had[iy]->Print(outputdir + "inclusive-jet.pdf","pdf");
    if(wait) canvas_pt_had[iy]->WaitPrimitive();

    delete hnull;

    delete legend_level;
    
    delete legend_mc;
    delete legend_mc_hadron;
    delete legend_mc_uncor;

    delete legend_y;
  }

  //----------------------//
  //-- Jet multiplicity --//
  //----------------------//

  TCanvas* canvas_num_jet = new TCanvas("canvas_num_jet","jet multiplicity",canvas_size,canvas_size);	
  canvas_num_jet->Divide(1,1);
  canvas_num_jet->SetLogy();
  ((TPad*) canvas_num_jet)->cd();
  ((TPad*) canvas_num_jet)->SetTicky();
  canvas_num_jet->SetLeftMargin(margin);  
  canvas_num_jet->SetRightMargin(margin); 

  histo_name = "efficiency/Multiplicity_DETJet";
  
  ymax = 0;
  
  //-- scaling
  for(int ifile = 0; ifile < Nfile; ifile++) 
    histo[ifile] = (TH1F*) file[ifile]->Get(histo_name);
  
  for(int ifile = 0; ifile < Nfile; ifile++) 
    scaling[ifile] = histo[1]->Integral()/histo[ifile]->Integral();
  
  //-- loop on data and monte carlo
  for(int ifile = 0; ifile < Nfile; ifile++) {
    histo[ifile] = (TH1F*) file[ifile]->Get(histo_name);
    histo[ifile]->Scale(scaling[ifile]);
    if(ifile == 1) SetHistoStyleData(histo[ifile],black,circle,marker_size_medium);
    else SetHistoStyleMC(histo[ifile],1,blue,line_width);
    
    if(histo[ifile]->GetBinContent(histo[ifile]->GetMaximumBin()) > ymax) ymax = histo[ifile]->GetBinContent(histo[ifile]->GetMaximumBin());
  }

  ymax*=10;
  hnull = new TH2D("hnull",title,100,0.8,10.2,100,10,ymax);   
  SetNullHisto(hnull,1.3,1.4,"jet multiplicity","number of events");
  hnull->DrawCopy();
  
  for(int ifile = 0; ifile < Nfile; ifile++) {
    if(ifile == 1) histo[ifile]->Draw("E1same");
    else histo[ifile]->Draw("histosame][");
  }

  legend_level = new TLegend(0.20,0.82,0.25,0.87);	     
  SetLegendStyle(legend_level,black);			       
  legend_level->AddEntry("","detector level","");          
  legend_level->Draw("same");                                
  
  legend_data = new TLegend(0.50,0.82,0.70,0.87);     
  SetLegendStyle(legend_data,black);		      
  legend_data->AddEntry(histo[1]," data","lp");       
  legend_data->Draw("same");                          
  
  legend_mc = new TLegend(0.50,0.77,0.70,0.82);
  SetLegendStyle(legend_mc,blue);
  legend_mc->AddEntry(histo[0]," PYTHIA8 CUETM1","lp");     
  legend_mc->Draw("same");                                           
  
  ((TPad*) canvas_num_jet)->RedrawAxis();
    
  canvas_num_jet->Print(outputdir + "jet-multiplicity.pdf");
  canvas_num_jet->Print(outputdir + "inclusive-jet.pdf","pdf");
  if(wait) canvas_num_jet->WaitPrimitive();

  delete hnull;

  delete legend_level;

  delete legend_mc;
  
  delete legend_data;

  //----------------//
  //-- Missing ET --//
  //----------------//
  
  TCanvas* canvas_MET = new TCanvas("canvas_MET","Missing ET",canvas_size,canvas_size);	
  canvas_MET->Divide(1,1);
  canvas_MET->SetLogy();
  ((TPad*) canvas_MET)->cd();
  ((TPad*) canvas_MET)->SetTicky();
  canvas_MET->SetLeftMargin(margin);  
  canvas_MET->SetRightMargin(margin); 
  
  histo_name = "efficiency/MET_DET";
  
  ymax = 0;
  
  //-- scaling
  for(int ifile = 0; ifile < Nfile; ifile++) 
    histo[ifile] = (TH1F*) file[ifile]->Get(histo_name);
  
  for(int ifile = 0; ifile < Nfile; ifile++) 
    scaling[ifile] = histo[1]->Integral()/histo[ifile]->Integral();
  
  //-- loop on data and monte carlo
  for(int ifile = 0; ifile < Nfile; ifile++) {
    histo[ifile] = (TH1F*) file[ifile]->Get(histo_name);
    histo[ifile]->Scale(scaling[ifile]);
    if(ifile == 1) SetHistoStyleData(histo[ifile],black,circle,marker_size_medium);
    else SetHistoStyleMC(histo[ifile],1,blue,line_width);
    
    if(histo[ifile]->GetBinContent(histo[ifile]->GetMaximumBin()) > ymax) ymax = histo[ifile]->GetBinContent(histo[ifile]->GetMaximumBin());
  }

  ymax*=10;
  hnull = new TH2D("hnull",title,100,0,220,100,200,ymax);   
  SetNullHisto(hnull,1.3,1.4,"missing E_{ T} [GeV]","number of events");
  hnull->DrawCopy();
  
  for(int ifile = 0; ifile < Nfile; ifile++) {
    if(ifile == 1) histo[ifile]->Draw("E1same");
    else histo[ifile]->Draw("histosame][");
  }

  legend_level = new TLegend(0.20,0.82,0.25,0.87);	     
  SetLegendStyle(legend_level,black);			       
  legend_level->AddEntry("","detector level","");          
  legend_level->Draw("same");                                
  
  legend_data = new TLegend(0.50,0.82,0.70,0.87);     
  SetLegendStyle(legend_data,black);		      
  legend_data->AddEntry(histo[1]," data","lp");       
  legend_data->Draw("same");                          
  
  legend_mc = new TLegend(0.50,0.77,0.70,0.82);
  SetLegendStyle(legend_mc,blue);
  legend_mc->AddEntry(histo[0]," PYTHIA8 CUETM1","lp");     
  legend_mc->Draw("same");                                           
  
  ((TPad*) canvas_MET)->RedrawAxis();
    
  canvas_MET->Print(outputdir + "MET.pdf");
  canvas_MET->Print(outputdir + "inclusive-jet.pdf","pdf");
  if(wait) canvas_MET->WaitPrimitive();

  delete hnull;

  delete legend_level;

  delete legend_mc;
  
  delete legend_data;

  //-------------------------//
  //-- Missing ET fraction --//
  //-------------------------//
  
  TCanvas* canvas_MET_fraction = new TCanvas("canvas_MET_fraction","Missing ET fraction",canvas_size,canvas_size);	
  canvas_MET_fraction->Divide(1,1);
  canvas_MET_fraction->SetLogy();
  ((TPad*) canvas_MET_fraction)->cd();
  ((TPad*) canvas_MET_fraction)->SetTicky();
  canvas_MET_fraction->SetLeftMargin(margin);  
  canvas_MET_fraction->SetRightMargin(margin); 
  
  histo_name = "efficiency/FractionMET_DET";
  
  ymax = 0;
  
  //-- scaling
  for(int ifile = 0; ifile < Nfile; ifile++) 
    histo[ifile] = (TH1F*) file[ifile]->Get(histo_name);
  
  for(int ifile = 0; ifile < Nfile; ifile++) 
    scaling[ifile] = histo[1]->Integral()/histo[ifile]->Integral();
  
  //-- loop on data and monte carlo
  for(int ifile = 0; ifile < Nfile; ifile++) {
    histo[ifile] = (TH1F*) file[ifile]->Get(histo_name);
    histo[ifile]->Scale(scaling[ifile]);
    if(ifile == 1) SetHistoStyleData(histo[ifile],black,circle,marker_size_medium);
    else SetHistoStyleMC(histo[ifile],1,blue,line_width);
    
    if(histo[ifile]->GetBinContent(histo[ifile]->GetMaximumBin()) > ymax) ymax = histo[ifile]->GetBinContent(histo[ifile]->GetMaximumBin());
  }

  ymax*=10;
  hnull = new TH2D("hnull",title,100,0,0.2,100,100,ymax);   
  SetNullHisto(hnull,1.3,1.4,"missing E_{ T} fraction","number of events");
  hnull->DrawCopy();
  
  for(int ifile = 0; ifile < Nfile; ifile++) {
    if(ifile == 1) histo[ifile]->Draw("E1same");
    else histo[ifile]->Draw("histosame][");
  }

  legend_level = new TLegend(0.20,0.82,0.25,0.87);	     
  SetLegendStyle(legend_level,black);			       
  legend_level->AddEntry("","detector level","");          
  legend_level->Draw("same");                                

  legend_data = new TLegend(0.50,0.82,0.70,0.87);     
  SetLegendStyle(legend_data,black);		      
  legend_data->AddEntry(histo[1]," data","lp");       
  legend_data->Draw("same");                          
  
  legend_mc = new TLegend(0.50,0.77,0.70,0.82);
  SetLegendStyle(legend_mc,blue);
  legend_mc->AddEntry(histo[0]," PYTHIA8 CUETM1","lp");     
  legend_mc->Draw("same");                                           
  
  ((TPad*) canvas_MET_fraction)->RedrawAxis();
    
  canvas_MET_fraction->Print(outputdir + "MET-fraction.pdf");
  canvas_MET_fraction->Print(outputdir + "inclusive-jet.pdf","pdf");
  if(wait) canvas_MET_fraction->WaitPrimitive();

  delete hnull;

  delete legend_level;

  delete legend_mc;
  
  delete legend_data;

  //------------------------------//
  //-- good vertex multiplicity --//
  //------------------------------//
  
  TCanvas* canvas_num_vertex = new TCanvas("canvas_num_vertex","vertex multiplicity",canvas_size,canvas_size);	
  canvas_num_vertex->Divide(1,1);
  canvas_num_vertex->SetLogy();
  ((TPad*) canvas_num_vertex)->cd();
  ((TPad*) canvas_num_vertex)->SetTicky();
  canvas_num_vertex->SetLeftMargin(margin);  
  canvas_num_vertex->SetRightMargin(margin); 
  
  histo_name = "efficiency/num_of_VtxGood";
  
  ymax = 0;
  
  //-- scaling
  for(int ifile = 0; ifile < Nfile; ifile++) 
    histo[ifile] = (TH1F*) file[ifile]->Get(histo_name);
  
  for(int ifile = 0; ifile < Nfile; ifile++) 
    scaling[ifile] = histo[1]->Integral()/histo[ifile]->Integral();
  
  //-- loop on data and monte carlo
  for(int ifile = 0; ifile < Nfile; ifile++) {
    histo[ifile] = (TH1F*) file[ifile]->Get(histo_name);
    histo[ifile]->Scale(scaling[ifile]);
    if(ifile == 1) SetHistoStyleData(histo[ifile],black,circle,marker_size_medium);
    else SetHistoStyleMC(histo[ifile],1,blue,line_width);
    
    if(histo[ifile]->GetBinContent(histo[ifile]->GetMaximumBin()) > ymax) ymax = histo[ifile]->GetBinContent(histo[ifile]->GetMaximumBin());
  }

  ymax*=10;
  hnull = new TH2D("hnull",title,100,1,25,100,100,ymax);   
  SetNullHisto(hnull,1.3,1.4,"vertex multiplicity","number of events");
  hnull->DrawCopy();
  
  for(int ifile = 0; ifile < Nfile; ifile++) {
    if(ifile == 1) histo[ifile]->Draw("E1same");
    else histo[ifile]->Draw("histosame][");
  }

  legend_level = new TLegend(0.20,0.82,0.25,0.87);	     
  SetLegendStyle(legend_level,black);			       
  legend_level->AddEntry("","detector level","");          
  legend_level->Draw("same");                              
  
  legend_data = new TLegend(0.50,0.82,0.70,0.87);     
  SetLegendStyle(legend_data,black);		      
  legend_data->AddEntry(histo[1]," data","lp");       
  legend_data->Draw("same");                          
  
  legend_mc = new TLegend(0.50,0.77,0.70,0.82);
  SetLegendStyle(legend_mc,blue);
  legend_mc->AddEntry(histo[0]," PYTHIA8 CUETM1","lp");     
  legend_mc->Draw("same");                                           
  
  ((TPad*) canvas_num_vertex)->RedrawAxis();
    
  canvas_num_vertex->Print(outputdir + "vertex-multiplicity.pdf");
  canvas_num_vertex->Print(outputdir + "inclusive-jet.pdf","pdf");
  if(wait) canvas_num_vertex->WaitPrimitive();

  delete hnull;

  delete legend_level;
  
  delete legend_mc;
  
  delete legend_data;

  //----------------------------------------//
  //-- pile up versus vertex multiplicity --//
  //----------------------------------------//
  
  TCanvas* canvas_pu_vertex = new TCanvas("canvas_pu_vertex","pile up versus vertex multiplicity",canvas_size,canvas_size);	
  canvas_pu_vertex->Divide(1,1);
  canvas_pu_vertex->cd();
  canvas_pu_vertex->SetLogz();
  canvas_pu_vertex->SetLeftMargin(margin);
  canvas_pu_vertex->SetRightMargin(margin);
  
  histo_name = "efficiency/PileUpVSVertex";
  
  histoPUVX = (TH2F*) file[0]->Get(histo_name);
  SetHistoStyleMC(histoPUVX,1,blue,line_width);
    
  hnull = new TH2D("hnull",title,100,-0.5,30.5,100,-0.5,30.5);   
  SetNullHisto(hnull,1.3,1.4,"pile up","vertex multiplicity");
  hnull->SetTitle("PYTHIA8 CUETM1");
  hnull->DrawCopy();
  
  histoPUVX->Draw("colz same");
    
  ((TPad*) canvas_pu_vertex)->RedrawAxis();
    
  canvas_pu_vertex->Print(outputdir + "vertex-pu.pdf");
  canvas_pu_vertex->Print(outputdir + "inclusive-jet.pdf","pdf");
  if(wait) canvas_pu_vertex->WaitPrimitive();

  delete hnull;

  //------------------------//
  //-- trigger efficiency --//
  //------------------------//

  TString HLT[Ntrigger] = {"60","80","140","200","260","320","400","450","500"};
  TF1* fit_effi[Ntrigger];
  
  TCanvas* canvas_trigger[Ntrigger];
   
  //-- loop over trigger
  for(int itrig = 0; itrig < Ntrigger; itrig++) {

    canvas_label = "canvas_trigger_" + TString::Format("%d",itrig+1);
    canvas_title = "trigger efficiency - " + TString::Format("%d",itrig+1);
    
    canvas_trigger[itrig] = new TCanvas(canvas_label,canvas_title,canvas_size,canvas_size);	
    canvas_trigger[itrig]->Divide(1,1);
    canvas_trigger[itrig]->SetLogx();
    ((TPad*) canvas_trigger[itrig])->cd();
    ((TPad*) canvas_trigger[itrig])->SetTicky();
    canvas_trigger[itrig]->SetLeftMargin(margin);  
    canvas_trigger[itrig]->SetRightMargin(margin); 
    
    //-- data only
    
    histo_name = "efficiency/leading_pt_HLTeffi_Jet" + HLT[itrig];
    heffi[itrig] = (TH1F*) file[1]->Get(histo_name);
    SetHistoStyleData(heffi[itrig],black,circle,marker_size_medium);

    if(itrig == 0) hnull = new TH2D("hnull",title,100,40,400,100,0,1.4);
    else if(itrig == 1 || itrig == 2) hnull = new TH2D("hnull",title,100,60,1000,100,0,1.4);
    else if(itrig == 3) hnull = new TH2D("hnull",title,100,150,1000,100,0,1.4);
    else if(itrig == 4) hnull = new TH2D("hnull",title,100,200,1000,100,0,1.4);
    else if(itrig == 5) hnull = new TH2D("hnull",title,100,280,1000,100,0,1.4);
    else if(itrig == 6) hnull = new TH2D("hnull",title,100,350,1000,100,0,1.4);
    else if(itrig == 7) hnull = new TH2D("hnull",title,100,400,1000,100,0,1.4);
    else if(itrig == 8) hnull = new TH2D("hnull",title,100,460,1000,100,0,1.4);
    
    else hnull = new TH2D("hnull",title,100,0,3000,100,0,1.4);   
    SetNullHisto(hnull,1.3,1.4,"p_{T} corrected [GeV]","trigger efficiency");
    hnull->DrawCopy();

    heffi[itrig]->Draw("E1same");

    double xmin = 0, xmax = 3000;

    if(itrig == 0) { xmin = 50; xmax = 200;}
    else if(itrig == 1) {xmin = 60; xmax = 500;}
    else if(itrig == 2) {xmin = 110; xmax = 500;}
    else if(itrig == 3) {xmin = 100; xmax = 600;}
    else if(itrig == 4) {xmin = 220; xmax = 600;}
    else if(itrig == 5) {xmin = 100; xmax = 800;}
    else if(itrig == 6) {xmin = 100; xmax = 800;}
    else if(itrig == 7) {xmin = 100; xmax = 800;}
    else if(itrig == 8) {xmin = 300; xmax = 800;}
    
    fit_effi[itrig] = FitSigmoid(heffi[itrig],xmin,xmax);
    fit_effi[itrig]->Draw("same");                         
    
    legend_level = new TLegend(0.20,0.82,0.25,0.87);	 
    SetLegendStyle(legend_level,black);			 
    legend_level->AddEntry("","data detector level","");             
    legend_level->Draw("same");                              
    
    legend_trigger[itrig] = new TLegend(0.50,0.82,0.70,0.87);	    
    SetLegendStyle(legend_trigger[itrig],black);
    TString trigger_name = " HLT" + HLT[itrig] + "U";
    legend_trigger[itrig]->AddEntry(heffi[itrig],trigger_name,"lp");     
    legend_trigger[itrig]->Draw("same");                              
     
    ((TPad*) canvas_trigger[itrig])->RedrawAxis();
    
    canvas_trigger[itrig]->Print(outputdir + "trigger-efficiency-" + TString::Format("%d",itrig+1) + ".pdf");
    canvas_trigger[itrig]->Print(outputdir + "inclusive-jet.pdf","pdf");
    if(wait) canvas_trigger[itrig]->WaitPrimitive();

    delete hnull;

    delete legend_level;

    delete legend_trigger[itrig];   
  }

  //-------------------//
  //-- trigger curve --//
  //-------------------//
  
  TCanvas* canvas_trigger_curve = new TCanvas("canvas_trigger_curve","trigger efficiency curve",canvas_size,canvas_size);	
  canvas_trigger_curve->Divide(1,1);
  ((TPad*) canvas_trigger_curve)->cd();
  ((TPad*) canvas_trigger_curve)->SetTicky();
  canvas_trigger_curve->SetLeftMargin(margin);  
  canvas_trigger_curve->SetRightMargin(margin); 
    
  hnull = new TH2D("hnull",title,100,40,700,100,0,1.8);
  SetNullHisto(hnull,1.3,1.4,"p_{T} corrected [GeV]","trigger efficiency");
  hnull->DrawCopy();

  for(int itrig = 0; itrig < Ntrigger; itrig++) {
    fit_effi[itrig]->SetLineColor(itrig+1);
    fit_effi[itrig]->SetLineStyle(itrig+1);
    fit_effi[itrig]->Draw("same");

    if(itrig < 3) {
      legend_trigger[itrig] = new TLegend(0.18,0.72-itrig*0.05,0.38,0.77-itrig*0.05);	    
      SetLegendStyle(legend_trigger[itrig],black);
      TString trigger_name = " HLT" + HLT[itrig] + "U";
      legend_trigger[itrig]->AddEntry(fit_effi[itrig],trigger_name,"l");     
      legend_trigger[itrig]->Draw("same");
    }

    else if(itrig < 6) {
      legend_trigger[itrig] = new TLegend(0.40,0.72-(itrig-3)*0.05,0.60,0.77-(itrig-3)*0.05);	    
      SetLegendStyle(legend_trigger[itrig],black);
      TString trigger_name = " HLT" + HLT[itrig] + "U";
      legend_trigger[itrig]->AddEntry(fit_effi[itrig],trigger_name,"l");     
      legend_trigger[itrig]->Draw("same");
    }

    else if(itrig < 9) {
      legend_trigger[itrig] = new TLegend(0.62,0.72-(itrig-6)*0.05,0.82,0.77-(itrig-6)*0.05);	    
      SetLegendStyle(legend_trigger[itrig],black);
      TString trigger_name = " HLT" + HLT[itrig] + "U";
      legend_trigger[itrig]->AddEntry(fit_effi[itrig],trigger_name,"l");     
      legend_trigger[itrig]->Draw("same");
    }

      
  }

  legend_level = new TLegend(0.20,0.82,0.25,0.87);	 
  SetLegendStyle(legend_level,black);			 
  legend_level->AddEntry("","data detector level","");             
  legend_level->Draw("same");                              
  
  ((TPad*) canvas_trigger_curve)->RedrawAxis();
    
  canvas_trigger_curve->Print(outputdir + "trigger-efficiency-curve.pdf");
  canvas_trigger_curve->Print(outputdir + "inclusive-jet.pdf","pdf");
  if(wait) canvas_trigger_curve->WaitPrimitive();
  
  delete hnull;
  
  delete legend_level;

  for(int itrig = 0; itrig < Ntrigger; itrig++) delete legend_trigger[itrig];   

  //------------------------------//
  //-- R - delta y - delta phi  --//
  //------------------------------//
  
  TCanvas* canvas_Rpt = new TCanvas("canvas_Rpt","R - delta y - delta phi",canvas_size,canvas_size);	
  canvas_Rpt->Divide(1,1);
  canvas_Rpt->SetLogy();
  ((TPad*) canvas_Rpt)->cd();
  ((TPad*) canvas_Rpt)->SetTicky();
  canvas_Rpt->SetLeftMargin(margin);  
  canvas_Rpt->SetRightMargin(margin); 
  
  histo_name = "efficiency/DeltaRpt";
  
  hR[0] = (TH1F*) file[0]->Get(histo_name);
  hR[0]->Rebin(100);
  SetHistoStyleMC(hR[0],1,blue,line_width);

  histo_name = "efficiency/DeltaRptmin";	   
  					   
  hR[1] = (TH1F*) file[0]->Get(histo_name);
  hR[1]->Rebin(100);			     
  SetHistoStyleData(hR[1],black,circle,marker_size_medium);

  ymax = 0;
  if(hR[0]->GetBinContent(hR[0]->GetMaximumBin()) > ymax) ymax = hR[0]->GetBinContent(hR[0]->GetMaximumBin());
  if(hR[1]->GetBinContent(hR[1]->GetMaximumBin()) > ymax) ymax = hR[1]->GetBinContent(hR[1]->GetMaximumBin());
  
  ymax*=10;
  hnull = new TH2D("hnull",title,100,0,6,100,20,ymax);   
  SetNullHisto(hnull,1.3,1.4,"#Delta R(#Delta y, #Delta #varphi)","number of jet pairs");
  hnull->DrawCopy();
  
  hR[0]->Draw("histosame][");
  hR[1]->Draw("E1same");
  
  line = new TLine(0.2,0,0.2,ymax);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw("same");                
  
  legend_level = new TLegend(0.25,0.82,0.30,0.87);	     
  SetLegendStyle(legend_level,black);			       
  legend_level->AddEntry("","PYTHIA8 CUETM1","");          
  legend_level->Draw("same");                              
  
  legend_mc = new TLegend(0.55,0.82,0.75,0.87);
  SetLegendStyle(legend_mc,blue);
  legend_mc->AddEntry(hR[0]," #Delta R(#Delta y, #Delta #varphi)","lp");     
  legend_mc->Draw("same");

  legend_data = new TLegend(0.55,0.77,0.75,0.82);				  
  SetLegendStyle(legend_data,black);					  
  legend_data->AddEntry(hR[1]," #Delta R(#Delta y, #Delta #varphi) min","lp");  
  legend_data->Draw("same");                                                  
  
  ((TPad*) canvas_Rpt)->RedrawAxis();
    
  canvas_Rpt->Print(outputdir + "delta-R-pt.pdf");
  canvas_Rpt->Print(outputdir + "inclusive-jet.pdf","pdf");
  if(wait) canvas_Rpt->WaitPrimitive();

  delete hnull;

  delete legend_level;
  
  delete legend_mc;

  delete legend_data;

  delete line;

  //-------------------------------//
  //-- R - delta pt - delta phi  --//
  //-------------------------------//
  
  TCanvas* canvas_Ry = new TCanvas("canvas_Ry","R - delta pT - delta phi",canvas_size,canvas_size);	
  canvas_Ry->Divide(1,1);
  canvas_Ry->SetLogy();
  ((TPad*) canvas_Ry)->cd();
  ((TPad*) canvas_Ry)->SetTicky();
  canvas_Ry->SetLeftMargin(margin);  
  canvas_Ry->SetRightMargin(margin); 
  
  histo_name = "efficiency/DeltaRy";
  
  hR[0] = (TH1F*) file[0]->Get(histo_name);
  hR[0]->Rebin(100);
  SetHistoStyleMC(hR[0],1,blue,line_width);

  histo_name = "efficiency/DeltaRymin";	   
  					   
  hR[1] = (TH1F*) file[0]->Get(histo_name);
  hR[1]->Rebin(100);			     
  SetHistoStyleData(hR[1],black,circle,marker_size_medium);

  ymax = 0;
  if(hR[0]->GetBinContent(hR[0]->GetMaximumBin()) > ymax) ymax = hR[0]->GetBinContent(hR[0]->GetMaximumBin());
  if(hR[1]->GetBinContent(hR[1]->GetMaximumBin()) > ymax) ymax = hR[1]->GetBinContent(hR[1]->GetMaximumBin());
  
  ymax*=10;
  hnull = new TH2D("hnull",title,100,0,4,100,20,ymax);   
  SetNullHisto(hnull,1.3,1.4,"#Delta R(#Delta p_{T}, #Delta #varphi)","number of jet pairs");
  hnull->DrawCopy();
  
  hR[0]->Draw("histosame][");
  hR[1]->Draw("E1same");

  line = new TLine(0.4,0,0.4,ymax);  
  line->SetLineStyle(2);	       
  line->SetLineWidth(2);	     
  line->Draw("same");                  
  
  legend_level = new TLegend(0.25,0.82,0.30,0.87);	     
  SetLegendStyle(legend_level,black);			       
  legend_level->AddEntry("","PYTHIA8 CUETM1","");          
  legend_level->Draw("same");                              
  
  legend_mc = new TLegend(0.55,0.82,0.75,0.87);
  SetLegendStyle(legend_mc,blue);
  legend_mc->AddEntry(hR[0]," #Delta R(#Delta p_{T}, #Delta #varphi)","lp");     
  legend_mc->Draw("same");

  legend_data = new TLegend(0.55,0.77,0.75,0.82);				  
  SetLegendStyle(legend_data,black);					  
  legend_data->AddEntry(hR[1]," #Delta R(#Delta p_{T}, #Delta #varphi) min","lp");  
  legend_data->Draw("same");                                                  
  
  ((TPad*) canvas_Ry)->RedrawAxis();
    
  canvas_Ry->Print(outputdir + "delta-R-y.pdf");
  canvas_Ry->Print(outputdir + "inclusive-jet.pdf","pdf");
  if(wait) canvas_Ry->WaitPrimitive();

  delete hnull;

  delete legend_level;
  
  delete legend_mc;

  delete legend_data;

  delete line;

  //------------------------------------//
  //-- resolution pT versus y - curve --//
  //------------------------------------//

  double sigmapt_y[Ny];
  double sigmapt_y_error[Ny];

  double ymean[Ny];
  double ymean_error[Ny];
  
  //-- loop over y bins
  for(int iy = 0; iy < Ny; iy++) {

    histo_name = "efficiency/ResolutionPtvsY_" + TString::Format("%d",iy+1) + "bin";

    histoResoPtvsY[iy] = (TH1F*) file[0]->Get(histo_name);
    histoResoPtvsY[iy]->Rebin(2);

    histo_name = "efficiency/y_GENJetReso_" + TString::Format("%d",iy+1) + "bin";
    histoYmean[iy] = (TH1F*) file[0]->Get(histo_name);
    ymean[iy] = histoYmean[iy]->GetMean();
    ymean_error[iy] = 0;
    
    //-- Gaussian
    fit_y_GA[iy] = FitGaussian(histoResoPtvsY[iy]);
    fit_y_GA[iy]->Draw("same");
    chi2_y_GA[iy] = fit_y_GA[iy]->GetChisquare()/fit_y_GA[iy]->GetNDF();
    
    sigmapt_y[iy] = 100*fit_y_GA[iy]->GetParameter(2);
    sigmapt_y_error[iy] = 100*fit_y_GA[iy]->GetParError(2);
    
    for(int ipar = 0; ipar < 3; ipar++)                                                                                                                              
      cout<<fit_y_GA[iy]->GetParName(ipar)<<" = "<<fit_y_GA[iy]->GetParameter(ipar)<<" +/- "<<fit_y_GA[iy]->GetParError(ipar)<<endl;

  }
  
  TCanvas* canvas_sigpt_y_curve = new TCanvas("canvas_sigpt_y_curve","pT resolution verus y",canvas_size,canvas_size);	
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
  
  legend_level = new TLegend(0.17,0.82,0.22,0.87);	 
  SetLegendStyle(legend_level,black);			 
  legend_level->AddEntry("","PYTHIA8 CUETM1","");             
  legend_level->Draw("same");

  ((TPad*) canvas_sigpt_y_curve)->RedrawAxis();
  
  canvas_sigpt_y_curve->Print(outputdir + "resolution-pt-vs-y-curve.pdf");
  canvas_sigpt_y_curve->Print(outputdir + "inclusive-jet.pdf","pdf");
  if(wait) canvas_sigpt_y_curve->WaitPrimitive();

  delete hnull;
  
  delete legend_level;

  //-------------------------------------//
  //-- resolution pT versus pT - curve --//
  //-------------------------------------//

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
        
    histo_name = "efficiency/ResolutionPtvsPt_" + TString::Format("%d",ipt+1) + "bin";

    histoResoPtvsPt[ipt] = (TH1F*) file[0]->Get(histo_name);    
    histoResoPtvsPt[ipt]->Rebin(2);
    
    bool filled = histoResoPtvsPt[ipt]->Integral() > 0.01;
    
    histo_name = "efficiency/pt_GENJetReso_" + TString::Format("%d",ipt+1) + "bin";
    histoPtmean[ipt] = (TH1F*) file[0]->Get(histo_name);
    if(histoPtmean[ipt]->Integral() > 0) ptmean[ipt] = histoPtmean[ipt]->GetMean();
    ptmean_error[ipt] = histoPtmean[ipt]->GetMeanError();;
	
    if(filled) {

      //-- Gaussian
      fit_pt_GA[ipt] = FitGaussian(histoResoPtvsPt[ipt]);
      chi2_pt_GA[ipt] = fit_pt_GA[ipt]->GetChisquare()/fit_pt_GA[ipt]->GetNDF();
      
      if(chi2_pt_GA[ipt] < CutFit) fit_pt_GA[ipt]->Draw("same");
      
      sigmapt_pt[ipt] = 100*fit_pt_GA[ipt]->GetParameter(2);
      sigmapt_pt_error[ipt] = 100*fit_pt_GA[ipt]->GetParError(2);
      
      for(int ipar = 0; ipar < 3; ipar++) 															     
	cout<<fit_pt_GA[ipt]->GetParName(ipar)<<" = "<<fit_pt_GA[ipt]->GetParameter(ipar)<<" +/- "<<fit_pt_GA[ipt]->GetParError(ipar)<<endl;    
    }
    
  }

  //-- for(int ipt = 0; ipt < Npt; ipt++) 
  //--  cout<<ptbin[ipt]<<" <  pT < "<<ptbin[ipt+1]<<" GeV: pT mean = "<<ptmean[ipt]<<" GeV - pT resolution "<<sigmapt_pt[ipt]<<" +/- "<<sigmapt_pt_error[ipt]<<" %"<<endl;

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
  canvas_sigpt_pt_curve->Print(outputdir + "inclusive-jet.pdf","pdf");
  if(wait) canvas_sigpt_pt_curve->WaitPrimitive();

  delete hnull;
  
  delete legend_mc;

  //--------------------------//
  //-- pT response versus y --//
  //--------------------------//

  TCanvas* canvas_pt_response[Ny];

  //-- loop over y bins
  for(int iy = 0; iy < Ny; iy++) {

    canvas_label = "canvas_pt_response_" + TString::Format("%d",iy+1);
    canvas_title = "pT response versus y - " + TString::Format("%d",iy+1);

    canvas_pt_response[iy] = new TCanvas(canvas_label,canvas_title,canvas_size,canvas_size);	
    canvas_pt_response[iy]->Divide(1,1);
    canvas_pt_response[iy]->cd();
    canvas_pt_response[iy]->SetLogx();
    canvas_pt_response[iy]->SetLogy();
    canvas_pt_response[iy]->SetLogz();
    canvas_pt_response[iy]->SetLeftMargin(margin*1.2);
    canvas_pt_response[iy]->SetRightMargin(margin*1.2);
  
    histo_name = "efficiency/pt_JetMatched_" + TString::Format("%d",iy+1) + "bin";
  
    histoResponse[iy] = (TH2F*) file[0]->Get(histo_name);
    SetHistoStyleMC(histoResponse[iy],1,blue,line_width);
    NormalizeTH2Column(histoResponse[iy]);
    histoResponse[iy]->GetZaxis()->SetRangeUser(0.01,100);

    if(iy == 2) hnull = new TH2D("hnull",title,100,50,2500,100,50,2500);
    else if(iy == 3) hnull = new TH2D("hnull",title,100,50,2000,100,50,2000);
    else if(iy == 4) hnull = new TH2D("hnull",title,100,50,1400,100,50,1400);
    else if(iy == 5) hnull = new TH2D("hnull",title,100,50,800,100,50,800);
    else if(iy == 6) hnull = new TH2D("hnull",title,100,50,400,100,50,400);
    else hnull = new TH2D("hnull",title,100,50,3000,100,50,3000);
    SetNullHisto(hnull,1.3,1.4,"p_{T} reco [GeV]","p_{T} gen [GeV]");
    hnull->SetTitle("PYTHIA8 CUETM1 - " + ybin[iy]);
    hnull->DrawCopy();
  
    histoResponse[iy]->Draw("colz same");
    
    ((TPad*) canvas_pt_response[iy])->RedrawAxis();
    
    canvas_pt_response[iy]->Print(outputdir + "pt-response-" + TString::Format("%d",iy+1) + ".pdf");
    canvas_pt_response[iy]->Print(outputdir + "inclusive-jet.pdf","pdf");
    if(wait) canvas_pt_response[iy]->WaitPrimitive();

    delete hnull;

  }  

  //--------------------//
  //-- Stability - pT --//
  //--------------------//
  
  TCanvas* canvas_pt_stability[Ny];
  
  //-- loop over y bins
  for(int iy = 0; iy < Ny; iy++) {

    canvas_label = "canvas_pt_stability_" + TString::Format("%d",iy+1);
    canvas_title = "pT stability - " + TString::Format("%d",iy+1);
    
    canvas_pt_stability[iy] = new TCanvas(canvas_label,canvas_title,canvas_size,canvas_size);	
    canvas_pt_stability[iy]->Divide(1,1);
    canvas_pt_stability[iy]->SetLogx();
    ((TPad*) canvas_pt_stability[iy])->cd();
    ((TPad*) canvas_pt_stability[iy])->SetTicky();
    canvas_pt_stability[iy]->SetLeftMargin(margin);  
    canvas_pt_stability[iy]->SetRightMargin(margin); 
        
    histo_name = "efficiency/pt_GENJetStability_" + TString::Format("%d",iy+1) + "bin";

    ymax = 0;

    //-- loop on monte carlo
    for(int ifile = 0; ifile < Nfile; ifile++) {
      if(ifile == 1) continue;
      histoPt[ifile][iy] = (TH1F*) file[ifile]->Get(histo_name);
      histoPt[ifile][iy]->Scale(100);
      SetHistoStyleMC(histoPt[ifile][iy],1,blue,line_width);

      if(histoPt[ifile][iy]->GetBinContent(histoPt[ifile][iy]->GetMaximumBin()) > ymax) ymax = histoPt[ifile][iy]->GetBinContent(histoPt[ifile][iy]->GetMaximumBin());
    }

    
    ymax*=1.2;
    double ptmax = 0;
    
    if(iy == 0 || iy == 1 || iy == 2) ptmax = 1950;
    if(iy == 3) ptmax = 1450;
    if(iy == 4) ptmax = 1000;
    if(iy == 5) ptmax = 780;
    if(iy == 6) ptmax = 340;
    
    hnull = new TH2D("hnull",title,100,ptmin,ptmax,100,0,100);   
    SetNullHisto(hnull,1.3,1.4,"p_{T} gen [GeV]","Stability [%]");
    hnull->DrawCopy();

    for(int ifile = 0; ifile < Nfile; ifile++) {
      if(ifile == 1) continue;
      histoPt[ifile][iy]->Draw("histosame][");
    }

    legend_level = new TLegend(0.17,0.82,0.22,0.87);	 
    SetLegendStyle(legend_level,black);			 
    legend_level->AddEntry("","PYTHIA8 CUETM1","");             
    legend_level->Draw("same");                             
        
    legend_y = new TLegend(0.48,0.82,0.68,0.87);     
    SetLegendStyle(legend_y,black);		       
    legend_y->AddEntry("",ybin[iy],"");          
    legend_y->Draw("same");                          

    ((TPad*) canvas_pt_stability[iy])->RedrawAxis();
    
    canvas_pt_stability[iy]->Print(outputdir + "pt-stability-" + TString::Format("%d",iy+1) + ".pdf");
    canvas_pt_stability[iy]->Print(outputdir + "inclusive-jet.pdf","pdf");
    if(wait) canvas_pt_stability[iy]->WaitPrimitive();

    delete hnull;

    delete legend_level;
    
    delete legend_y;
  }

  //---------------//
  //-- Miss - pT --//
  //---------------//

  TCanvas* canvas_pt_miss[Ny];
   
  //-- loop over y bins
  for(int iy = 0; iy < Ny; iy++) {

    canvas_label = "canvas_pt_miss_" + TString::Format("%d",iy+1);
    canvas_title = "pT miss - " + TString::Format("%d",iy+1);
    
    canvas_pt_miss[iy] = new TCanvas(canvas_label,canvas_title,canvas_size,canvas_size);	
    canvas_pt_miss[iy]->Divide(1,1);
    canvas_pt_miss[iy]->SetLogx();
    canvas_pt_miss[iy]->SetLogy();
    ((TPad*) canvas_pt_miss[iy])->cd();
    ((TPad*) canvas_pt_miss[iy])->SetTicky();
    canvas_pt_miss[iy]->SetLeftMargin(margin);  
    canvas_pt_miss[iy]->SetRightMargin(margin); 
        
    histo_name = "efficiency/pt_GENJetMiss_" + TString::Format("%d",iy+1) + "bin";

    ymax = 0;

    //-- loop on monte carlo
    for(int ifile = 0; ifile < Nfile; ifile++) {
      if(ifile == 1) continue;
      histoPt[ifile][iy] = (TH1F*) file[ifile]->Get(histo_name);
      histoPt[ifile][iy]->Scale(100);
      SetHistoStyleMC(histoPt[ifile][iy],1,red,line_width);

      if(histoPt[ifile][iy]->GetBinContent(histoPt[ifile][iy]->GetMaximumBin()) > ymax) ymax = histoPt[ifile][iy]->GetBinContent(histoPt[ifile][iy]->GetMaximumBin());
    }

    
    ymax*=1.2;
    double ptmax = 0;
    
    if(iy == 0 || iy == 1 || iy == 2) ptmax = 1780;
    if(iy == 3) ptmax = 1100;
    if(iy == 4) ptmax = 850;
    if(iy == 5) ptmax = 630;
    if(iy == 6) ptmax = 220;

    if(iy == 4) hnull = new TH2D("hnull",title,100,ptmin,ptmax,100,0.1,100);   
    else if(iy == 5) hnull = new TH2D("hnull",title,100,ptmin,ptmax,100,1,100);   
    else if(iy == 6) hnull = new TH2D("hnull",title,100,ptmin,ptmax,100,0.05,100);   
    else hnull = new TH2D("hnull",title,100,ptmin,ptmax,100,0.01,100);   
    SetNullHisto(hnull,1.3,1.4,"p_{T} gen [GeV]","Miss [%]");
    hnull->DrawCopy();

    for(int ifile = 0; ifile < Nfile; ifile++) {
      if(ifile == 1) continue;
      histoPt[ifile][iy]->Draw("histosame][");
    }

    legend_level = new TLegend(0.17,0.82,0.22,0.87);	 
    SetLegendStyle(legend_level,black);			 
    legend_level->AddEntry("","PYTHIA8 CUETM1","");             
    legend_level->Draw("same");                             
        
    legend_y = new TLegend(0.48,0.82,0.68,0.87);     
    SetLegendStyle(legend_y,black);		       
    legend_y->AddEntry("",ybin[iy],"");          
    legend_y->Draw("same");                          

    ((TPad*) canvas_pt_miss[iy])->RedrawAxis();
    
    canvas_pt_miss[iy]->Print(outputdir + "pt-miss-" + TString::Format("%d",iy+1) + ".pdf");
    canvas_pt_miss[iy]->Print(outputdir + "inclusive-jet.pdf","pdf");
    if(wait) canvas_pt_miss[iy]->WaitPrimitive();

    delete hnull;

    delete legend_level;
    
    delete legend_y;
  }

  //---------------------//
  //-- Acceptance - pT --//
  //---------------------//

  TCanvas* canvas_pt_acceptance[Ny];
   
  //-- loop over y bins
  for(int iy = 0; iy < Ny; iy++) {

    canvas_label = "canvas_pt_acceptance_" + TString::Format("%d",iy+1);
    canvas_title = "pT acceptance - " + TString::Format("%d",iy+1);
    
    canvas_pt_acceptance[iy] = new TCanvas(canvas_label,canvas_title,canvas_size,canvas_size);	
    canvas_pt_acceptance[iy]->Divide(1,1);
    canvas_pt_acceptance[iy]->SetLogx();
    //canvas_pt_acceptance[iy]->SetLogy();
    ((TPad*) canvas_pt_acceptance[iy])->cd();
    ((TPad*) canvas_pt_acceptance[iy])->SetTicky();
    canvas_pt_acceptance[iy]->SetLeftMargin(margin);  
    canvas_pt_acceptance[iy]->SetRightMargin(margin); 
        
    histo_name = "efficiency/pt_GENJetAcceptance_" + TString::Format("%d",iy+1) + "bin";

    ymax = 0;

    //-- loop on monte carlo
    for(int ifile = 0; ifile < Nfile; ifile++) {
      if(ifile == 1) continue;
      histoPt[ifile][iy] = (TH1F*) file[ifile]->Get(histo_name);
      histoPt[ifile][iy]->Scale(100);
      SetHistoStyleMC(histoPt[ifile][iy],1,blue,line_width);

      if(histoPt[ifile][iy]->GetBinContent(histoPt[ifile][iy]->GetMaximumBin()) > ymax) ymax = histoPt[ifile][iy]->GetBinContent(histoPt[ifile][iy]->GetMaximumBin());
    }

    
    ymax*=1.2;
    double ptmax = 0;
    
    if(iy == 0 || iy == 1 || iy == 2) ptmax = 1780;
    if(iy == 3) ptmax = 1100;
    if(iy == 4) ptmax = 850;
    if(iy == 5) ptmax = 630;
    if(iy == 6) ptmax = 220;
    
    hnull = new TH2D("hnull",title,100,ptmin,ptmax,100,50,120);   
    SetNullHisto(hnull,1.3,1.4,"p_{T} gen [GeV]","Acceptance [%]");
    hnull->DrawCopy();

    for(int ifile = 0; ifile < Nfile; ifile++) {
      if(ifile == 1) continue;
      histoPt[ifile][iy]->Draw("histosame][");
    }

    legend_level = new TLegend(0.17,0.82,0.22,0.87);	 
    SetLegendStyle(legend_level,black);			 
    legend_level->AddEntry("","PYTHIA8 CUETM1","");             
    legend_level->Draw("same");                             
        
    legend_y = new TLegend(0.48,0.82,0.68,0.87);     
    SetLegendStyle(legend_y,black);		       
    legend_y->AddEntry("",ybin[iy],"");          
    legend_y->Draw("same");                          

    ((TPad*) canvas_pt_acceptance[iy])->RedrawAxis();
    
    canvas_pt_acceptance[iy]->Print(outputdir + "pt-acceptance-" + TString::Format("%d",iy+1) + ".pdf");
    canvas_pt_acceptance[iy]->Print(outputdir + "inclusive-jet.pdf","pdf");
    if(wait) canvas_pt_acceptance[iy]->WaitPrimitive();

    delete hnull;

    delete legend_level;
    
    delete legend_y;
  }


  //-----------------//
  //-- Purity - pT --//
  //-----------------//

  TCanvas* canvas_pt_purity[Ny];
   
  //-- loop over y bins
  for(int iy = 0; iy < Ny; iy++) {

    canvas_label = "canvas_pt_purity_" + TString::Format("%d",iy+1);
    canvas_title = "pT purity - " + TString::Format("%d",iy+1);
    
    canvas_pt_purity[iy] = new TCanvas(canvas_label,canvas_title,canvas_size,canvas_size);	
    canvas_pt_purity[iy]->Divide(1,1);
    canvas_pt_purity[iy]->SetLogx();
    ((TPad*) canvas_pt_purity[iy])->cd();
    ((TPad*) canvas_pt_purity[iy])->SetTicky();
    canvas_pt_purity[iy]->SetLeftMargin(margin);  
    canvas_pt_purity[iy]->SetRightMargin(margin); 
        
    histo_name = "efficiency/pt_DETJetPurity_" + TString::Format("%d",iy+1) + "bin";

    ymax = 0;

    //-- loop on monte carlo
    for(int ifile = 0; ifile < Nfile; ifile++) {
      if(ifile == 1) continue;
      histoPt[ifile][iy] = (TH1F*) file[ifile]->Get(histo_name);
      histoPt[ifile][iy]->Scale(100);
      SetHistoStyleMC(histoPt[ifile][iy],1,red,line_width);

      if(histoPt[ifile][iy]->GetBinContent(histoPt[ifile][iy]->GetMaximumBin()) > ymax) ymax = histoPt[ifile][iy]->GetBinContent(histoPt[ifile][iy]->GetMaximumBin());
    }

    
    ymax*=1.2;
    double ptmax = 0;
    
    if(iy == 0 || iy == 1 || iy == 2) ptmax = 1950;
    if(iy == 3) ptmax = 1450;
    if(iy == 4) ptmax = 1000;
    if(iy == 5) ptmax = 780;
    if(iy == 6) ptmax = 340;
    
    hnull = new TH2D("hnull",title,100,ptmin,ptmax,100,0,100);   
    SetNullHisto(hnull,1.3,1.4,"p_{T} reco [GeV]","Purity [%]");
    hnull->DrawCopy();

    for(int ifile = 0; ifile < Nfile; ifile++) {
      if(ifile == 1) continue;
      histoPt[ifile][iy]->Draw("histosame][");
    }

    legend_level = new TLegend(0.17,0.82,0.22,0.87);	 
    SetLegendStyle(legend_level,black);			 
    legend_level->AddEntry("","PYTHIA8 CUETM1","");             
    legend_level->Draw("same");                             
        
    legend_y = new TLegend(0.48,0.82,0.68,0.87);     
    SetLegendStyle(legend_y,black);		       
    legend_y->AddEntry("",ybin[iy],"");          
    legend_y->Draw("same");                          

    ((TPad*) canvas_pt_purity[iy])->RedrawAxis();
    
    canvas_pt_purity[iy]->Print(outputdir + "pt-purity-" + TString::Format("%d",iy+1) + ".pdf");
    canvas_pt_purity[iy]->Print(outputdir + "inclusive-jet.pdf","pdf");
    if(wait) canvas_pt_purity[iy]->WaitPrimitive();

    delete hnull;

    delete legend_level;
    
    delete legend_y;
  }

  //-----------------//
  //-- Fake - pT --//
  //-----------------//

  TCanvas* canvas_pt_fake[Ny];
   
  //-- loop over y bins
  for(int iy = 0; iy < Ny; iy++) {

    canvas_label = "canvas_pt_fake_" + TString::Format("%d",iy+1);
    canvas_title = "pT fake - " + TString::Format("%d",iy+1);
    
    canvas_pt_fake[iy] = new TCanvas(canvas_label,canvas_title,canvas_size,canvas_size);	
    canvas_pt_fake[iy]->Divide(1,1);
    canvas_pt_fake[iy]->SetLogx();
    //    canvas_pt_fake[iy]->SetLogy();
    ((TPad*) canvas_pt_fake[iy])->cd();
    ((TPad*) canvas_pt_fake[iy])->SetTicky();
    canvas_pt_fake[iy]->SetLeftMargin(margin);  
    canvas_pt_fake[iy]->SetRightMargin(margin); 
        
    histo_name = "efficiency/pt_DETJetFake_" + TString::Format("%d",iy+1) + "bin";

    ymax = 0;

    //-- loop on monte carlo
    for(int ifile = 0; ifile < Nfile; ifile++) {
      if(ifile == 1) continue;
      histoPt[ifile][iy] = (TH1F*) file[ifile]->Get(histo_name);
      histoPt[ifile][iy]->Scale(100);
      SetHistoStyleMC(histoPt[ifile][iy],1,blue,line_width);

      if(histoPt[ifile][iy]->GetBinContent(histoPt[ifile][iy]->GetMaximumBin()) > ymax) ymax = histoPt[ifile][iy]->GetBinContent(histoPt[ifile][iy]->GetMaximumBin());
    }

    
    ymax*=1.2;
    double ptmax = 0;
    
    if(iy == 0 || iy == 1 || iy == 2) ptmax = 150;
    if(iy == 3) ptmax = 150;
    if(iy == 4) ptmax = 200;
    if(iy == 5) ptmax = 150;
    if(iy == 6) ptmax = 150;
    
    hnull = new TH2D("hnull",title,100,ptmin,ptmax,100,0,120);   
    SetNullHisto(hnull,1.3,1.4,"p_{T} reco [GeV]","Fake [%]");
    hnull->DrawCopy();

    for(int ifile = 0; ifile < Nfile; ifile++) {
      if(ifile == 1) continue;
      histoPt[ifile][iy]->Draw("histosame][");
    }

    legend_level = new TLegend(0.17,0.82,0.22,0.87);	 
    SetLegendStyle(legend_level,black);			 
    legend_level->AddEntry("","PYTHIA8 CUETM1","");             
    legend_level->Draw("same");                             
        
    legend_y = new TLegend(0.48,0.82,0.68,0.87);     
    SetLegendStyle(legend_y,black);		       
    legend_y->AddEntry("",ybin[iy],"");          
    legend_y->Draw("same");                          

    ((TPad*) canvas_pt_fake[iy])->RedrawAxis();
    
    canvas_pt_fake[iy]->Print(outputdir + "pt-fake-" + TString::Format("%d",iy+1) + ".pdf");
    if(iy < Ny-1) canvas_pt_fake[iy]->Print(outputdir + "inclusive-jet.pdf","pdf");
    else canvas_pt_fake[iy]->Print(outputdir + "inclusive-jet.pdf)","pdf");
    if(wait) canvas_pt_fake[iy]->WaitPrimitive();

    delete hnull;

    delete legend_level;
    
    delete legend_y;
  }

  return(0);
}
