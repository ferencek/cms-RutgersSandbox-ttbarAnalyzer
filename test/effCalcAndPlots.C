#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TF1.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
// #include "exoStyle.C"

using namespace std;


void calc_eff_diLept(const string& fInputFile, const string& fHisto)
{
  TFile *file = new TFile(fInputFile.c_str());

  TH2D *h1_histo = (TH2D*)file->Get(fHisto.c_str());

  double nTotal    = h1_histo->GetBinContent(1);
  double nDiLepton = h1_histo->GetBinContent(3);
  double nZVeto    = h1_histo->GetBinContent(4);
  double n2Jets    = h1_histo->GetBinContent(5);
  double nMET      = h1_histo->GetBinContent(6);
  double n1BTag    = h1_histo->GetBinContent(7);

  cout << nDiLepton/nTotal << std::endl
       << nZVeto/nTotal << std::endl
       << n2Jets/nTotal << std::endl
       << nMET/nTotal << std::endl
       << n1BTag/nTotal << std::endl << std::endl;
}


void calc_eff_semiLept(const string& fInputFile, const string& fHisto)
{
  TFile *file = new TFile(fInputFile.c_str());

  TH2D *h1_histo = (TH2D*)file->Get(fHisto.c_str());

  double nTotal  = h1_histo->GetBinContent(1);
  double nLepton = h1_histo->GetBinContent(2);
  double n4Jets  = h1_histo->GetBinContent(3);
  double n1BTag  = h1_histo->GetBinContent(4);

  cout << nLepton/nTotal << std::endl
       << n4Jets/nTotal << std::endl
       << n1BTag/nTotal << std::endl << std::endl;
}


void plot_dR()
{
  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle();

  TFile *file = new TFile("TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola__Summer12_DR53X-PU_S10_START53_V7C-v1__AODSIM.root");

  TH2D *h1_histo = (TH2D*)file->Get("analyzer/h1_CutFlow_diLept_mumu");
  double nTotal = h1_histo->GetBinContent(1);

  TH2D *h2_dRbJetJet_elel = (TH2D*)file->Get("analyzer/h1_mindRbJetJet_diLept_elel");
  TH2D *h2_dRbJetJet_mumu = (TH2D*)file->Get("analyzer/h1_mindRbJetJet_diLept_mumu");
  TH2D *h2_dRbJetJet_elmu = (TH2D*)file->Get("analyzer/h1_mindRbJetJet_diLept_elmu");

  TH2D *h2_dRbJetbJet_elel = (TH2D*)file->Get("analyzer/h1_mindRbJetbJet_diLept_elel");
  TH2D *h2_dRbJetbJet_mumu = (TH2D*)file->Get("analyzer/h1_mindRbJetbJet_diLept_mumu");
  TH2D *h2_dRbJetbJet_elmu = (TH2D*)file->Get("analyzer/h1_mindRbJetbJet_diLept_elmu");

  h2_dRbJetJet_elel->Add(h2_dRbJetJet_mumu);
  h2_dRbJetJet_elel->Add(h2_dRbJetJet_elmu);

  h2_dRbJetbJet_elel->Add(h2_dRbJetbJet_mumu);
  h2_dRbJetbJet_elel->Add(h2_dRbJetbJet_elmu);

  h2_dRbJetJet_elel->Scale( (227.*20000.)/nTotal );
  h2_dRbJetbJet_elel->Scale( (227.*20000.)/nTotal );

  TCanvas *c = new TCanvas("c", "",1000,800);
  c->cd();

  h2_dRbJetJet_elel->Draw("hist");

  c->SaveAs("h1_dRbJetJet_ttbar_diLept.eps");

  h2_dRbJetbJet_elel->Draw("hist");

  c->SaveAs("h1_dRbJetbJet_ttbar_diLept.eps");

  delete c;
  delete file;
}


void performCalculations()
{
  calc_eff_diLept("TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola__Summer12_DR53X-PU_S10_START53_V7C-v1__AODSIM.root", "analyzer/h1_CutFlow_diLept_mumu");
  calc_eff_diLept("TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola__Summer12_DR53X-PU_S10_START53_V7C-v1__AODSIM.root", "analyzer/h1_CutFlow_diLept_elel");
  calc_eff_diLept("TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola__Summer12_DR53X-PU_S10_START53_V7C-v1__AODSIM.root", "analyzer/h1_CutFlow_diLept_elmu");

  calc_eff_semiLept("TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola__Summer12_DR53X-PU_S10_START53_V7C-v1__AODSIM.root", "analyzer/h1_CutFlow_semiLept_mu");
  calc_eff_semiLept("TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola__Summer12_DR53X-PU_S10_START53_V7C-v1__AODSIM.root", "analyzer/h1_CutFlow_semiLept_el");
}

void makePlots()
{
  plot_dR();
}