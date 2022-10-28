//root
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TMath.h"
#include "TFree.h"
#include "TBranch.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"

//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>

#include <time.h>

using namespace std;

//https://en.wikipedia.org/wiki/Barometric_formula
//
const int n = 7;
//const double h_b[7] = {      0.00, 11000.00, 20000.00, 32000.00, 47000.00, 51000.00, 71000.00};
//const double p_b[7] = { 101325.00, 22632.10,  5474.89,   868.02,   110.91,    66.94,     3.96};
//const double T_b[7] = {    288.15,   216.65,   216.65,   228.65,   270.65,   270.65,   214.65};
//const double L_b[7] = {   -0.0065,      0.0,    0.001,   0.0028,      0.0,  -0.0028,   -0.002};

const double p_max = 1.5e5;
//const double p_min = 1.0e2;
const double p_min = 1.0e-5;

const double rho_max = 2.0;
//const double rho_min = 1.0e-2;
const double rho_min = 1.0e-20;

const double h_min = 0.0;
const double h_max = 175000.0; //m
const int n_sp = 1000000;
//const double h_max = 30000.0; //m
const double h_b[7] =   {      0.00, 11000.00, 20000.00, 32000.00, 47000.00, 51000.00, 71000.00};
const double p_b[7] =   { 101325.00, 22632.10,  5474.89,   868.02,   110.91,    66.94,     3.96};
const double rho_b[7] = {    1.2250,  0.36391,  0.08803,  0.01322,  0.00143,  0.00086, 0.000064};
const double T_b[7] =   {    288.15,   216.65,   216.65,   228.65,   270.65,   270.65,   214.65};
const double Tc_b[7] =   {288.15-273.15,   216.65-273.15,   216.65-273.15,   228.65-273.15,   270.65-273.15,   270.65-273.15,   214.65-273.15};
const double L_b[7] =   {   -0.0065,    0.001,    0.001,   0.0028,  -0.0028,  -0.0028,   -0.002};

//
const Double_t Re = 6371300;       //m
const Double_t Me = 5.9742*1.0e24; //kg

const double p0 = 101325;      // Sea level standard atmospheric pressure 101325 Pa
const double L = 0.00976;      // Temperature lapse rate, = g/cp for dry air ~ 0.00976 K/m
const double cp = 1004.68506;  // Constant-pressure specific heat 1004.68506 J/(kg·K)
const double T0 = 288.16;      // Sea level standard temperature 288.16 K
const double g0 = 9.80665;     // Earth-surface gravitational acceleration 9.80665 m/s2
const double M = 0.02896968;   // Molar mass of dry air 0.02896968 kg/mol
const double R0 = 8.314462618; // Universal gas constant 8.314462618 J/(mol·K)

double g_gravity(Double_t R);

//h Height above Earth-surface m
double get_atm_preasure(double h);
double get_atm_rho(double h);
double get_atm_approx01(double h);
double get_atm_approx02(double h);
void get_atm_model_from_corsika(TString filen, TGraph *gr_corsika_rho, TGraph *gr_corsika_thick, TGraph *gr_corsika_n, TGraph *gr_corsika_T, TGraph *gr_corsika_p, TGraph *gr_corsika_pw);

struct sphere_layer_str {
  //
  double h;
  double d_r;
  double impact_parameter; 
  double rho_atm; //kg/m/m/m
  //
  double r;
  double r_min;
  double r_max;
  double x_int_min;
  double y_int_min;
  double x_int_max;
  double y_int_max;
  double d_l;
  double g_cm2; //g/cm/cm
  //
  sphere_layer_str(){
    h = 0.0;
    d_r = 0.0;
    impact_parameter = 0.0;
    rho_atm = 0.0;
    //
    r = 0.0;
    r_min = 0.0;
    r_max = 0.0;
    x_int_min = 0.0;
    y_int_min = 0.0;
    y_int_min = 0.0;
    y_int_max = 0.0;
    d_l = 0.0;
    g_cm2 = 0.0;
  }
  //
  sphere_layer_str( double alt, double d_alt, double impact_par, double rho_atmosphere){
    h = alt;
    d_r = d_alt;
    impact_parameter = impact_par;
    rho_atm = rho_atmosphere;
    //
    r = Re + h;
    r_min = r - d_r/2.0;
    r_max = r + d_r/2.0;
    if(!get_intersection( r_min, impact_parameter, x_int_min, y_int_min)){
      x_int_min = -999.0;
      y_int_min = -999.0;
    }
    if(!get_intersection( r_max, impact_parameter, x_int_max, y_int_max)){
      x_int_max = -999.0;
      y_int_max = -999.0;
    }
    d_l = get_dl(x_int_min, y_int_min, x_int_max, y_int_max);
    g_cm2 = get_g_per_cm2(d_l,rho_atm);
  }
  //
  bool get_intersection(double r,double impactpar,double &x_int,double &y_int){
    if(impactpar>=0.0 && impactpar<=r){
      x_int = impactpar;
      y_int = TMath::Sqrt(r*r - x_int*x_int);
      return true;
    }
    else{
      x_int = -999.0;
      y_int = -999.0;
      return false;
    }
  }
  //
  double get_dl(double x_int_s,double y_int_s,double x_int_e,double y_int_e){
    if(x_int_s < 0.0 || x_int_e < 0.0 || y_int_s < 0.0 || y_int_e < 0.0)
      return 0.0;
    if(x_int_s<=Re)
      return TMath::Sqrt((x_int_s - x_int_e)*(x_int_s - x_int_e) + (y_int_s - y_int_e)*(y_int_s - y_int_e));
    else
      return 2*TMath::Sqrt((x_int_s - x_int_e)*(x_int_s - x_int_e) + (y_int_s - y_int_e)*(y_int_s - y_int_e));
  }
  //
  double get_g_per_cm2(double d_l,double rho_atm){
    return d_l*rho_atm/10.0; 
  }
  double get_r(double x,double y){
    return TMath::Sqrt(x*x + y*y);
  }
  //
  void printInfo(){
    cout<<"h                "<<h<<endl
      	<<"d_r              "<<d_r<<endl
	<<"impact_parameter "<<impact_parameter<<endl
	<<"rho_atm          "<<rho_atm<<endl;
    cout<<"r         "<<r<<endl
	<<"r_min     "<<r_min<<endl
	<<"r_max     "<<r_max<<endl
	<<"x_int_min "<<x_int_min<<endl
	<<"y_int_min "<<y_int_min<<endl
	<<"y_int_min "<<y_int_min<<endl
	<<"y_int_max "<<y_int_max<<endl
	<<"d_l       "<<d_l<<endl
      	<<"g_cm2     "<<g_cm2<<endl;
  }
};

int main(){
  //
  //
  //
  TGraph *gr_corsika_rho = new TGraph();
  TGraph *gr_corsika_thick = new TGraph();
  TGraph *gr_corsika_n = new TGraph();
  TGraph *gr_corsika_T = new TGraph();
  TGraph *gr_corsika_p = new TGraph();
  TGraph *gr_corsika_pw = new TGraph();
  gr_corsika_rho->SetNameTitle("gr_corsika_rho","gr_corsika_rho");
  gr_corsika_thick->SetNameTitle("gr_corsika_thick","gr_corsika_thick");
  gr_corsika_n->SetNameTitle("gr_corsika_n","gr_corsika_n");
  gr_corsika_T->SetNameTitle("gr_corsika_T","gr_corsika_T");
  gr_corsika_p->SetNameTitle("gr_corsika_p","gr_corsika_p");
  gr_corsika_pw->SetNameTitle("gr_corsika_pw","gr_corsika_pw");
  get_atm_model_from_corsika("./atmospheric_model/atmprof_ecmwf_north_winter_nexp.dat", gr_corsika_rho, gr_corsika_thick, gr_corsika_n, gr_corsika_T, gr_corsika_p, gr_corsika_pw);
  //
  //
  //
  int nn = 10000;
  double h = 0.0;
  //
  TGraph *gr_pa = new TGraph();
  gr_pa->SetNameTitle("gr_pa","gr_pa");
  TGraph *gr_pa_approx01 = new TGraph();
  gr_pa_approx01->SetNameTitle("gr_pa_approx01","gr_pa_approx01");
  TGraph *gr_pa_approx02 = new TGraph();
  gr_pa_approx02->SetNameTitle("gr_pa_approx02","gr_pa_approx02");
  TGraph *gr_rho = new TGraph();
  gr_rho->SetNameTitle("gr_rho","gr_rho");
  TGraph *gr_t = new TGraph(n,h_b,Tc_b);
  gr_t->SetNameTitle("gr_t","gr_t");
  //
  for(Int_t i = 0;i<nn;i++){
    h = h_min + (h_max-h_min)/(nn-1)*i;
    gr_pa->SetPoint(gr_pa->GetN(),h,get_atm_preasure(h));
    gr_pa_approx01->SetPoint(gr_pa_approx01->GetN(),h,get_atm_approx01(h));
    gr_pa_approx02->SetPoint(gr_pa_approx02->GetN(),h,get_atm_approx02(h));
    gr_rho->SetPoint(gr_rho->GetN(),h,get_atm_rho(h));
  }
  //
  gr_pa->SetLineColor(kBlack);
  gr_pa_approx01->SetLineColor(kRed);
  gr_pa_approx02->SetLineColor(kBlue);
  gr_corsika_p->SetLineColor(kGreen+2);
  //
  gr_pa->SetMarkerColor(kBlack);
  gr_pa_approx01->SetMarkerColor(kRed);
  gr_pa_approx02->SetMarkerColor(kBlue);
  gr_corsika_p->SetMarkerColor(kGreen+2);
  //
  gr_pa->SetLineWidth(2.0);
  gr_pa_approx01->SetLineWidth(2.0);
  gr_pa_approx02->SetLineWidth(2.0);
  gr_corsika_p->SetLineWidth(2.0);
  //
  //gr_pa->SetMinimum(0.0);
  //gr_pa->Draw("APL");
  TCanvas *c1 = new TCanvas("c1","c1",10,10,800,800);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetOptStat(kFALSE);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  c1->SetRightMargin(0.07);
  c1->SetLeftMargin(0.12);
  c1->SetTopMargin(0.02);
  c1->SetBottomMargin(0.08);
  //
  TMultiGraph *mg01 = new TMultiGraph();
  mg01->Add(gr_pa);
  mg01->Add(gr_pa_approx01);
  mg01->Add(gr_pa_approx02);
  mg01->Add(gr_corsika_p);
  mg01->SetMaximum(p_max);
  mg01->SetMinimum(p_min);
  mg01->Draw("APL");
  mg01->GetYaxis()->SetTitle("Atmospheric pressure, Pa");
  mg01->GetXaxis()->SetTitle("Altitude, m");
  //
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  //
  TLegend *leg01 = new TLegend(0.50,0.85,0.92,0.97,"","brNDC");
  leg01->AddEntry(gr_pa,          "Barometric formula", "pl");
  leg01->AddEntry(gr_pa_approx01, "Approximation 1 for small altitudes (<8 km)", "pl");
  leg01->AddEntry(gr_pa_approx02, "Approximation 2 for small altitudes (<8 km), exp", "pl");
  leg01->AddEntry(gr_corsika_p, "Corsika atmospheric model (ECMWF)", "pl");
  leg01->Draw();
  //
  TCanvas *c2 = new TCanvas("c2","c2",10,10,800,800);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetOptStat(kFALSE);
  //
  c2->SetRightMargin(0.07);
  c2->SetLeftMargin(0.13);
  c2->SetTopMargin(0.02);
  c2->SetBottomMargin(0.09);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  //
  gr_rho->SetMaximum(rho_max);
  gr_rho->SetMinimum(rho_min);
  //
  gr_rho->SetTitle("");
  gr_rho->SetLineColor(kBlack);
  gr_rho->SetMarkerColor(kBlack);
  gr_rho->SetLineWidth(2.0);
  //
  gr_corsika_rho->SetLineColor(kRed);
  gr_corsika_rho->SetMarkerColor(kRed);
  gr_corsika_rho->SetLineWidth(2.0);
  //gr_rho->Draw("APL");
  //gr_rho->GetYaxis()->SetTitle("Atmosphere density, kg/m^3");
  //gr_rho->GetXaxis()->SetTitle("Altitude, m");
  //
  //
  TMultiGraph *mg02 = new TMultiGraph();
  mg02->Add(gr_rho);
  mg02->Add(gr_corsika_rho);
  mg02->SetMaximum(rho_max);
  mg02->SetMinimum(rho_min);
  mg02->Draw("APL");
  mg02->GetYaxis()->SetTitle("Atmosphere density, kg/m^3");
  mg02->GetXaxis()->SetTitle("Altitude, m");
  //
  TLegend *leg02 = new TLegend(0.50,0.85,0.92,0.97,"","brNDC");
  leg02->AddEntry(gr_t, "Used in barometric formula", "pl");
  leg02->AddEntry(gr_corsika_T, "Corsika atmospheric model (ECMWF)", "pl");
  leg02->Draw();
  //
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  //
  TCanvas *c3 = new TCanvas("c3","c3",10,10,800,800);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetOptStat(kFALSE);
  //
  c3->SetRightMargin(0.03);
  c3->SetLeftMargin(0.13);
  c3->SetTopMargin(0.02);
  c3->SetBottomMargin(0.09);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  //
  gr_t->SetTitle("");
  gr_t->SetMarkerColor(kBlack);
  gr_t->SetLineColor(kBlack);
  gr_t->SetLineWidth(2.0);
  gr_corsika_T->SetMarkerColor(kRed);
  gr_corsika_T->SetLineColor(kRed);
  gr_corsika_T->SetLineWidth(2.0);
  //gr_t->SetMaximum(20.0);
  //gr_t->SetMinimum(-100);
  //gr_t->Draw("APL");
  //gr_t->GetYaxis()->SetTitle("Temperature, o^C");
  //gr_t->GetXaxis()->SetTitle("Altitude, m");
  //
  TMultiGraph *mg03 = new TMultiGraph();
  mg03->Add(gr_corsika_T);
  mg03->Add(gr_t);
  mg03->SetMaximum(65.0);
  mg03->SetMinimum(-100.0);
  mg03->Draw("APL");
  mg03->GetYaxis()->SetTitle("Temperature, o^C");
  mg03->GetXaxis()->SetTitle("Altitude, m");
  //
  TLegend *leg03 = new TLegend(0.50,0.85,0.92,0.97,"","brNDC");
  leg03->AddEntry(gr_t, "Used in barometric formula", "pl");
  leg03->AddEntry(gr_corsika_T, "Corsika atmospheric model (ECMWF)", "pl");
  leg03->Draw();
  //
  //
  //
  double alt;
  double d_alt;
  double impact_par;
  double rho_atmosphere;
  //
  int n_theta = 10;
  double theta;
  double theta_deg_min = 0.0;
  double theta_deg_max = 90.0;
  double theta_deg;
  //
  TGraph *gr_g_cm2_tot = new TGraph();
  gr_g_cm2_tot->SetNameTitle("gr_g_cm2_tot","gr_g_cm2_tot");
  TCanvas *c4 = new TCanvas("c4","c4",10,10,800,800);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetOptStat(kFALSE);
  //
  c4->SetRightMargin(0.03);
  c4->SetLeftMargin(0.13);
  c4->SetTopMargin(0.02);
  c4->SetBottomMargin(0.09);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  //
  for(int j = 0;j<n_theta;j++){
    if(j%10 == 0)
      cout<<j<<endl;
    theta_deg = theta_deg_min + (theta_deg_max - theta_deg_min)/(n_theta-1)*j;
    theta = theta_deg/180.0*TMath::Pi();
    impact_par=TMath::Sin(theta)*Re;
    d_alt = (h_max - h_min)/n_sp;
    vector<sphere_layer_str> v_sphere_layer;
    for(int i = 0;i<n_sp;i++){
      alt = d_alt*(1.0/2.0 + i);
      rho_atmosphere = get_atm_rho(alt);
      v_sphere_layer.push_back(sphere_layer_str(alt,d_alt,impact_par,rho_atmosphere));
    }
    double g_cm2_tot = 0.0;
    for(int i = 0;i<n_sp;i++){
      g_cm2_tot += v_sphere_layer.at(i).g_cm2;
      //v_sphere_layer.at(i).printInfo();
    }
    gr_g_cm2_tot->SetPoint(j,theta_deg,g_cm2_tot);
    //cout<<"g_cm2_tot "<<g_cm2_tot<<endl;
  }
  //
  gr_g_cm2_tot->SetTitle("");
  gr_g_cm2_tot->SetMarkerColor(kBlack);
  gr_g_cm2_tot->SetMarkerColor(kBlack);
  gr_g_cm2_tot->SetLineWidth(2.0);
  gr_g_cm2_tot->SetMaximum(40.0*1.0e+3);
  //gr_g_cm2_tot->SetMinimum(900.0);
  gr_g_cm2_tot->SetMinimum(0.0);
  gr_g_cm2_tot->Draw("APL");
  gr_g_cm2_tot->GetXaxis()->SetRangeUser(0.0,93.0);
  gr_g_cm2_tot->GetYaxis()->SetTitle("Air mass, g/cm/cm");
  gr_g_cm2_tot->GetXaxis()->SetTitle("Zenith, deg");
  //
  //
  //
  TGraph *gr_g_cm2_tot_vs_impact = new TGraph();
  gr_g_cm2_tot_vs_impact->SetNameTitle("gr_g_cm2_tot_vs_impact","gr_g_cm2_tot_vs_impact");
  //
  int n_impact = 2;
  double impact_par_min = 0.0;
  double impact_par_max = 6350000;
  for(int j = 0;j<n_impact;j++){
    if(j%10 == 0)
      cout<<j<<endl;
    impact_par=impact_par_min + (impact_par_max - impact_par_min)/(n_impact-1)*j;
    d_alt = (h_max - h_min)/n_sp;
    vector<sphere_layer_str> v_sphere_layer;
    for(int i = 0;i<n_sp;i++){
      alt = d_alt*(1.0/2.0 + i);
      rho_atmosphere = get_atm_rho(alt);
      v_sphere_layer.push_back(sphere_layer_str(alt,d_alt,impact_par,rho_atmosphere));
    }
    double g_cm2_tot = 0.0;
    for(int i = 0;i<n_sp;i++)
      g_cm2_tot += v_sphere_layer.at(i).g_cm2;
    gr_g_cm2_tot_vs_impact->SetPoint(gr_g_cm2_tot_vs_impact->GetN(),impact_par-Re,g_cm2_tot);
  }
  n_impact = 2;
  impact_par_min = impact_par_max;
  impact_par_max = 6450000;
  for(int j = 0;j<n_impact;j++){
    if(j%10 == 0)
      cout<<j<<endl;
    impact_par=impact_par_min + (impact_par_max - impact_par_min)/(n_impact-1)*j;
    d_alt = (h_max - h_min)/n_sp;
    vector<sphere_layer_str> v_sphere_layer;
    for(int i = 0;i<n_sp;i++){
      alt = d_alt*(1.0/2.0 + i);
      rho_atmosphere = get_atm_rho(alt);
      v_sphere_layer.push_back(sphere_layer_str(alt,d_alt,impact_par,rho_atmosphere));
    }
    double g_cm2_tot = 0.0;
    for(int i = 0;i<n_sp;i++)
      g_cm2_tot += v_sphere_layer.at(i).g_cm2;
    gr_g_cm2_tot_vs_impact->SetPoint(gr_g_cm2_tot_vs_impact->GetN(),impact_par-Re,g_cm2_tot);
  }
  //
  TCanvas *c5 = new TCanvas("c5","c5",10,10,800,800);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetOptStat(kFALSE);
  //
  c5->SetRightMargin(0.03);
  c5->SetLeftMargin(0.13);
  c5->SetTopMargin(0.02);
  c5->SetBottomMargin(0.09);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  //
  gr_g_cm2_tot_vs_impact->SetTitle("");
  gr_g_cm2_tot_vs_impact->SetMarkerColor(kBlack);
  gr_g_cm2_tot_vs_impact->SetMarkerColor(kBlack);
  gr_g_cm2_tot_vs_impact->SetLineWidth(2.0);
  //gr_g_cm2_tot_vs_impact->SetMaximum(40.0*1.0e+3);
  //gr_g_cm2_tot->SetMinimum(900.0);
  //gr_g_cm2_tot_vs_impact->SetMinimum(0.0);
  gr_g_cm2_tot_vs_impact->Draw("APL");
  //gr_g_cm2_tot_vs_impact->GetXaxis()->SetRangeUser(0.0,93.0);
  gr_g_cm2_tot_vs_impact->GetYaxis()->SetTitle("Air mass, g/cm/cm");
  gr_g_cm2_tot_vs_impact->GetXaxis()->SetTitle("Impact parameter, m");
  //
  //
  //
  TString rootFileName = "eas.root";
  TFile* rootFile = new TFile(rootFileName.Data(), "RECREATE", " Histograms", 1);
  rootFile->cd();
  if (rootFile->IsZombie()){
    std::cout<<"  ERROR ---> file "<<rootFileName.Data()<<" is zombi"<<std::endl;
    assert(0);
  }
  else {
    std::cout<<"  Output Histos file ---> "<<rootFileName.Data()<<std::endl;
  }
  ////////////////////////
  c1->Write();
  c2->Write();
  c3->Write();
  c4->Write();
  c5->Write();
  //
  gr_pa->Write();
  gr_pa_approx01->Write();
  gr_pa_approx02->Write();
  gr_rho->Write();
  gr_t->Write();
  gr_g_cm2_tot->Write();
  gr_g_cm2_tot_vs_impact->Write();
  //
  gr_corsika_rho->Write();
  gr_corsika_thick->Write();
  gr_corsika_n->Write();
  gr_corsika_T->Write();
  gr_corsika_p->Write();
  gr_corsika_pw->Write();
  //
  rootFile->Close();
  return 0;
}

//h Height above Earth-surface m
double get_atm_preasure(double h){
  if(h<0)
    return p_b[0];
  if(h>h_max)
    return 0.0;
  int index = 0;
  if(h>h_b[n-1]){
    index = n-1;
  }
  else{
    for(int i = 1;i<n;i++){
      if(h>=h_b[i-1] && h<h_b[i]){
	index = i-1;
	break;
      }
    }
  }
  double g = g_gravity(Re + h);
  return p_b[index]*TMath::Power((T_b[index] + (h-h_b[index])*L_b[index])/T_b[index],-g*M/R0/L_b[index]);
}

double get_atm_approx01(double h){
  double g = g_gravity(Re + h);
  if((1.0-g*h/cp/T0)<0.0)
    return 0.0;
  return p0*TMath::Power((1.0-g*h/cp/T0),cp*M/R0);
}

double get_atm_approx02(double h){
  double g = g_gravity(Re + h);
  return p0*TMath::Exp(-g*h*M/T0/R0);
}
  
double g_gravity(double R){
  return TMath::G()*Me/R/R;
}

double get_atm_rho(double h){
  if(h<0)
    return rho_b[0];
  if(h>h_max)
    return 0.0;
  int index = 0;
  if(h>h_b[n-1]){
    index = n-1;
  }
  else{
    for(int i = 1;i<n;i++){
      if(h>=h_b[i-1] && h<h_b[i]){
	index = i-1;
	break;
      }
    }
  }
  double g = g_gravity(Re + h);
  return rho_b[index]*TMath::Power(T_b[index]/(T_b[index] + (h-h_b[index])*L_b[index]),(1.0+g*M/R0/L_b[index]));
}

void get_atm_model_from_corsika(TString filen, TGraph *gr_corsika_rho, TGraph *gr_corsika_thick, TGraph *gr_corsika_n, TGraph *gr_corsika_T, TGraph *gr_corsika_p, TGraph *gr_corsika_pw){
  //# Atmospheric Model ECMWF year/month/day   hour h
  //# Col. #1          #2           #3            #4        [ #5 ]        [ #6 ]       [ # 7 ]
  //# Alt [km]    rho [g/cm^3] thick [g/cm^2]    n-1        T [k]       p [mbar]      pw / p
  //0.000     1.21788E-03  1.04061E+03  2.83924E-04  2.91910E+02  1.02049E+03  1.45399E-02
  double alt, rho, thick, n_m1, T, p, pw;
  string mot;
  ifstream filein(filen.Data());
  int i = 0;
  if(filein.is_open()){
    do{
      filein>>mot;
    }while(mot!="[mbar]");
    filein>>mot; filein>>mot; filein>>mot;
    //
    while(filein>>alt>>rho>>thick>>n_m1>>T>>p>>pw){
      //cout<<alt<<endl;
      alt *= 1000;
      gr_corsika_rho->SetPoint(i,alt,rho*1000);
      gr_corsika_thick->SetPoint(i,alt,thick);
      gr_corsika_n->SetPoint(i,alt,n_m1);
      gr_corsika_T->SetPoint(i,alt,T-273.15);
      gr_corsika_p->SetPoint(i,alt,p*100.0);
      gr_corsika_pw->SetPoint(i,alt,pw);
      i++;
    }
    //
    filein.close();
  }
  else cout <<"Unable to open file"<<endl; 
}
