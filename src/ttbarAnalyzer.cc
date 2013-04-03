// -*- C++ -*-
//
// Package:    ttbarAnalyzer
// Class:      ttbarAnalyzer
//
/**\class ttbarAnalyzer ttbarAnalyzer.cc RutgersSandbox/ttbarAnalyzer/src/ttbarAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  margaret zientek
//         Created:  Fri Feb  8 12:31:07 CST 2013
// $Id: ttbarAnalyzer.cc,v 1.5 2013/03/30 00:17:35 ferencek Exp $
//
// Analysis of a ttbar MC sample based on di-leptonic (TOP-12-007) and semi-leptonic (TOP-12-006) ttbar cross section measurements
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"

//
// class declaration
//

template <class T> struct ordering_Pt {
    ordering_Pt() {}
    bool operator ()(T const& a, T const& b) { return (a.first->pt() + a.second->pt()) > (b.first->pt() + b.second->pt()); }
};

class ttbarAnalyzer : public edm::EDAnalyzer {
   public:
      explicit ttbarAnalyzer(const edm::ParameterSet&);
      ~ttbarAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------

    const edm::InputTag         genParticleTag;
    const edm::InputTag         genJetTag;
    const edm::InputTag         jetFlavorTag;
    const edm::InputTag         genMetTag;
    const double                elecPtMinDiLept;
    const double                elecAbsEtaMaxDiLept;
    const double                elecPtMinSemiLept;
    const double                elecAbsEtaMaxSemiLept;
    const double                muonPtMinDiLept;
    const double                muonAbsEtaMaxDiLept;
    const double                muonPtMinSemiLept;
    const double                muonAbsEtaMaxSemiLept;
    const double                electronMuonDeltaRDiLept;
    const double                diLeptMassMin;
    const double                zVetoMin;
    const double                zVetoMax;
    const double                jetPtMinDiLept;
    const double                jetAbsEtaMaxDiLept;
    const double                jetLeptonDeltaRDiLept;
    const double                metCutDiLept;
    const double                jetPtMinLowSemiLept;
    const double                jetPtMinHighSemiLept;
    const double                jetAbsEtaMaxSemiLept;
    const double                jetLeptonDeltaRSemiLept;

    edm::Service<TFileService> fs;

    TH1D *h1_CutFlow_diLept_elel;
    TH1D *h1_CutFlow_diLept_mumu;
    TH1D *h1_CutFlow_diLept_elmu;
    TH1D *h1_CutFlow_semiLept_el;
    TH1D *h1_CutFlow_semiLept_mu;

    TH1D *h1_bTagMultiplicity_diLept_elel;
    TH1D *h1_bTagMultiplicity_diLept_mumu;
    TH1D *h1_bTagMultiplicity_diLept_elmu;
    TH1D *h1_bTagMultiplicity_semiLept_el;
    TH1D *h1_bTagMultiplicity_semiLept_mu;

    TH1D *h1_mindR_bJetJet_diLept_elel;
    TH1D *h1_mindR_bJetJet_diLept_mumu;
    TH1D *h1_mindR_bJetJet_diLept_elmu;
    TH1D *h1_mindR_bJetJet_semiLept_el;
    TH1D *h1_mindR_bJetJet_semiLept_mu;

    TH2D *h2_mindR_bJetPt_bJetJet_diLept_elel;
    TH2D *h2_mindR_bJetPt_bJetJet_diLept_mumu;
    TH2D *h2_mindR_bJetPt_bJetJet_diLept_elmu;
    TH2D *h2_mindR_bJetPt_bJetJet_semiLept_el;
    TH2D *h2_mindR_bJetPt_bJetJet_semiLept_mu;

    TH2D *h2_mindR_JetPt_bJetJet_diLept_elel;
    TH2D *h2_mindR_JetPt_bJetJet_diLept_mumu;
    TH2D *h2_mindR_JetPt_bJetJet_diLept_elmu;
    TH2D *h2_mindR_JetPt_bJetJet_semiLept_el;
    TH2D *h2_mindR_JetPt_bJetJet_semiLept_mu;

    TH2D *h2_mindR_FatJetPt_bJetJet_diLept_elel;
    TH2D *h2_mindR_FatJetPt_bJetJet_diLept_mumu;
    TH2D *h2_mindR_FatJetPt_bJetJet_diLept_elmu;
    TH2D *h2_mindR_FatJetPt_bJetJet_semiLept_el;
    TH2D *h2_mindR_FatJetPt_bJetJet_semiLept_mu;

    TH1D *h1_mindR_bJetbJet_diLept_elel;
    TH1D *h1_mindR_bJetbJet_diLept_mumu;
    TH1D *h1_mindR_bJetbJet_diLept_elmu;
    TH1D *h1_mindR_bJetbJet_semiLept_el;
    TH1D *h1_mindR_bJetbJet_semiLept_mu;

    TH2D *h2_mindR_bJet1Pt_bJetbJet_diLept_elel;
    TH2D *h2_mindR_bJet1Pt_bJetbJet_diLept_mumu;
    TH2D *h2_mindR_bJet1Pt_bJetbJet_diLept_elmu;
    TH2D *h2_mindR_bJet1Pt_bJetbJet_semiLept_el;
    TH2D *h2_mindR_bJet1Pt_bJetbJet_semiLept_mu;

    TH2D *h2_mindR_bJet2Pt_bJetbJet_diLept_elel;
    TH2D *h2_mindR_bJet2Pt_bJetbJet_diLept_mumu;
    TH2D *h2_mindR_bJet2Pt_bJetbJet_diLept_elmu;
    TH2D *h2_mindR_bJet2Pt_bJetbJet_semiLept_el;
    TH2D *h2_mindR_bJet2Pt_bJetbJet_semiLept_mu;

    TH2D *h2_mindR_FatJetPt_bJetbJet_diLept_elel;
    TH2D *h2_mindR_FatJetPt_bJetbJet_diLept_mumu;
    TH2D *h2_mindR_FatJetPt_bJetbJet_diLept_elmu;
    TH2D *h2_mindR_FatJetPt_bJetbJet_semiLept_el;
    TH2D *h2_mindR_FatJetPt_bJetbJet_semiLept_mu;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ttbarAnalyzer::ttbarAnalyzer(const edm::ParameterSet& iConfig) :

  genParticleTag(iConfig.getParameter<edm::InputTag>("GenParticleTag")),
  genJetTag(iConfig.getParameter<edm::InputTag>("GenJetTag")),
  jetFlavorTag(iConfig.getParameter<edm::InputTag>("JetFlavorTag")),
  genMetTag(iConfig.getParameter<edm::InputTag>("GenMetTag")),
  elecPtMinDiLept (iConfig.getParameter<double> ("ElecPtMinDiLept")),
  elecAbsEtaMaxDiLept (iConfig.getParameter<double> ("ElecAbsEtaMaxDiLept")),
  elecPtMinSemiLept (iConfig.getParameter<double> ("ElecPtMinSemiLept")),
  elecAbsEtaMaxSemiLept (iConfig.getParameter<double> ("ElecAbsEtaMaxSemiLept")),
  muonPtMinDiLept (iConfig.getParameter<double> ("MuonPtMinDiLept")),
  muonAbsEtaMaxDiLept (iConfig.getParameter<double> ("MuonAbsEtaMaxDiLept")),
  muonPtMinSemiLept (iConfig.getParameter<double> ("MuonPtMinSemiLept")),
  muonAbsEtaMaxSemiLept (iConfig.getParameter<double> ("MuonAbsEtaMaxSemiLept")),
  electronMuonDeltaRDiLept (iConfig.getParameter<double> ("ElectronMuonDeltaRDiLept")),
  diLeptMassMin (iConfig.getParameter<double> ("DiLeptMassMin")),
  zVetoMin (iConfig.getParameter<double> ("ZVetoMin")),
  zVetoMax (iConfig.getParameter<double> ("ZVetoMax")),
  jetPtMinDiLept (iConfig.getParameter<double> ("JetPtMinDiLept")),
  jetAbsEtaMaxDiLept (iConfig.getParameter<double> ("JetAbsEtaMaxDiLept")),
  jetLeptonDeltaRDiLept (iConfig.getParameter<double> ("JetLeptonDeltaRDiLept")),
  metCutDiLept (iConfig.getParameter<double> ("METCutDiLept")),
  jetPtMinLowSemiLept (iConfig.getParameter<double> ("JetPtMinLowSemiLept")),
  jetPtMinHighSemiLept (iConfig.getParameter<double> ("JetPtMinHighSemiLept")),
  jetAbsEtaMaxSemiLept (iConfig.getParameter<double> ("JetAbsEtaMaxSemiLept")),
  jetLeptonDeltaRSemiLept (iConfig.getParameter<double> ("JetLeptonDeltaRSemiLept"))

{
  //now do what ever initialization is needed
  int ncutBins=10;
  double cutmin=-0.5;
  double cutmax=9.5;

  h1_CutFlow_diLept_elel = fs->make<TH1D>("h1_CutFlow_diLept_elel",";Cut number;Number of events",ncutBins,cutmin,cutmax);
  h1_CutFlow_diLept_elel->Sumw2();
  h1_CutFlow_diLept_elel->SetDefaultSumw2(kTRUE); // to have TH1::Sumw2() automatically called for all subsequent histograms
  h1_CutFlow_diLept_mumu = fs->make<TH1D>("h1_CutFlow_diLept_mumu",";Cut number;Number of events",ncutBins,cutmin,cutmax);
  h1_CutFlow_diLept_elmu = fs->make<TH1D>("h1_CutFlow_diLept_elmu",";Cut number;Number of events",ncutBins,cutmin,cutmax);
  h1_CutFlow_semiLept_el = fs->make<TH1D>("h1_CutFlow_semiLept_el",";Cut number;Number of events",ncutBins,cutmin,cutmax);
  h1_CutFlow_semiLept_mu = fs->make<TH1D>("h1_CutFlow_semiLept_mu",";Cut number;Number of events",ncutBins,cutmin,cutmax);

  h1_bTagMultiplicity_diLept_elel = fs->make<TH1D>("h1_bTagMultiplicity_diLept_elel",";Number of b tags;Number of events",6,-0.5,5.5);
  h1_bTagMultiplicity_diLept_mumu = fs->make<TH1D>("h1_bTagMultiplicity_diLept_mumu",";Number of b tags;Number of events",6,-0.5,5.5);
  h1_bTagMultiplicity_diLept_elmu = fs->make<TH1D>("h1_bTagMultiplicity_diLept_elmu",";Number of b tags;Number of events",6,-0.5,5.5);
  h1_bTagMultiplicity_semiLept_el = fs->make<TH1D>("h1_bTagMultiplicity_semiLept_el",";Number of b tags;Number of events",6,-0.5,5.5);
  h1_bTagMultiplicity_semiLept_mu = fs->make<TH1D>("h1_bTagMultiplicity_semiLept_mu",";Number of b tags;Number of events",6,-0.5,5.5);

  h1_mindR_bJetJet_diLept_elel = fs->make<TH1D>("h1_mindR_bJetJet_diLept_elel",";min#DeltaR(b jet,jet);Number of entries",100,0.,5.);
  h1_mindR_bJetJet_diLept_mumu = fs->make<TH1D>("h1_mindR_bJetJet_diLept_mumu",";min#DeltaR(b jet,jet);Number of entries",100,0.,5.);
  h1_mindR_bJetJet_diLept_elmu = fs->make<TH1D>("h1_mindR_bJetJet_diLept_elmu",";min#DeltaR(b jet,jet);Number of entries",100,0.,5.);
  h1_mindR_bJetJet_semiLept_el = fs->make<TH1D>("h1_mindR_bJetJet_semiLept_el",";min#DeltaR(b jet,jet);Number of entries",100,0.,5.);
  h1_mindR_bJetJet_semiLept_mu = fs->make<TH1D>("h1_mindR_bJetJet_semiLept_mu",";min#DeltaR(b jet,jet);Number of entries",100,0.,5.);

  h2_mindR_bJetPt_bJetJet_diLept_elel = fs->make<TH2D>("h1_mindR_bJetPt_bJetJet_diLept_elel",";min#DeltaR(b jet,jet);b jet p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_bJetPt_bJetJet_diLept_mumu = fs->make<TH2D>("h1_mindR_bJetPt_bJetJet_diLept_mumu",";min#DeltaR(b jet,jet);b jet p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_bJetPt_bJetJet_diLept_elmu = fs->make<TH2D>("h1_mindR_bJetPt_bJetJet_diLept_elmu",";min#DeltaR(b jet,jet);b jet p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_bJetPt_bJetJet_semiLept_el = fs->make<TH2D>("h1_mindR_bJetPt_bJetJet_semiLept_el",";min#DeltaR(b jet,jet);b jet p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_bJetPt_bJetJet_semiLept_mu = fs->make<TH2D>("h1_mindR_bJetPt_bJetJet_semiLept_mu",";min#DeltaR(b jet,jet);b jet p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);

  h2_mindR_JetPt_bJetJet_diLept_elel = fs->make<TH2D>("h1_mindR_JetPt_bJetJet_diLept_elel",";min#DeltaR(b jet,jet);Jet p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_JetPt_bJetJet_diLept_mumu = fs->make<TH2D>("h1_mindR_JetPt_bJetJet_diLept_mumu",";min#DeltaR(b jet,jet);Jet p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_JetPt_bJetJet_diLept_elmu = fs->make<TH2D>("h1_mindR_JetPt_bJetJet_diLept_elmu",";min#DeltaR(b jet,jet);Jet p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_JetPt_bJetJet_semiLept_el = fs->make<TH2D>("h1_mindR_JetPt_bJetJet_semiLept_el",";min#DeltaR(b jet,jet);Jet p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_JetPt_bJetJet_semiLept_mu = fs->make<TH2D>("h1_mindR_JetPt_bJetJet_semiLept_mu",";min#DeltaR(b jet,jet);Jet p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);

  h2_mindR_FatJetPt_bJetJet_diLept_elel = fs->make<TH2D>("h1_mindR_FatJetPt_bJetJet_diLept_elel",";min#DeltaR(b jet,jet);Fat jet p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_FatJetPt_bJetJet_diLept_mumu = fs->make<TH2D>("h1_mindR_FatJetPt_bJetJet_diLept_mumu",";min#DeltaR(b jet,jet);Fat jet p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_FatJetPt_bJetJet_diLept_elmu = fs->make<TH2D>("h1_mindR_FatJetPt_bJetJet_diLept_elmu",";min#DeltaR(b jet,jet);Fat jet p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_FatJetPt_bJetJet_semiLept_el = fs->make<TH2D>("h1_mindR_FatJetPt_bJetJet_semiLept_el",";min#DeltaR(b jet,jet);Fat jet p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_FatJetPt_bJetJet_semiLept_mu = fs->make<TH2D>("h1_mindR_FatJetPt_bJetJet_semiLept_mu",";min#DeltaR(b jet,jet);Fat jet p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  
  h1_mindR_bJetbJet_diLept_elel = fs->make<TH1D>("h1_mindR_bJetbJet_diLept_elel",";min#DeltaR(b jet,b jet);Number of entries",100,0.,5.);
  h1_mindR_bJetbJet_diLept_mumu = fs->make<TH1D>("h1_mindR_bJetbJet_diLept_mumu",";min#DeltaR(b jet,b jet);Number of entries",100,0.,5.);
  h1_mindR_bJetbJet_diLept_elmu = fs->make<TH1D>("h1_mindR_bJetbJet_diLept_elmu",";min#DeltaR(b jet,b jet);Number of entries",100,0.,5.);
  h1_mindR_bJetbJet_semiLept_el = fs->make<TH1D>("h1_mindR_bJetbJet_semiLept_el",";min#DeltaR(b jet,b jet);Number of entries",100,0.,5.);
  h1_mindR_bJetbJet_semiLept_mu = fs->make<TH1D>("h1_mindR_bJetbJet_semiLept_mu",";min#DeltaR(b jet,b jet);Number of entries",100,0.,5.);

  h2_mindR_bJet1Pt_bJetbJet_diLept_elel = fs->make<TH2D>("h1_mindR_bJet1Pt_bJetbJet_diLept_elel",";min#DeltaR(b jet,b jet);b jet_{1} p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_bJet1Pt_bJetbJet_diLept_mumu = fs->make<TH2D>("h1_mindR_bJet1Pt_bJetbJet_diLept_mumu",";min#DeltaR(b jet,b jet);b jet_{1} p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_bJet1Pt_bJetbJet_diLept_elmu = fs->make<TH2D>("h1_mindR_bJet1Pt_bJetbJet_diLept_elmu",";min#DeltaR(b jet,b jet);b jet_{1} p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_bJet1Pt_bJetbJet_semiLept_el = fs->make<TH2D>("h1_mindR_bJet1Pt_bJetbJet_semiLept_el",";min#DeltaR(b jet,b jet);b jet_{1} p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_bJet1Pt_bJetbJet_semiLept_mu = fs->make<TH2D>("h1_mindR_bJet1Pt_bJetbJet_semiLept_mu",";min#DeltaR(b jet,b jet);b jet_{1} p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);

  h2_mindR_bJet2Pt_bJetbJet_diLept_elel = fs->make<TH2D>("h1_mindR_bJet2Pt_bJetbJet_diLept_elel",";min#DeltaR(b jet,b jet);b jet_{2} p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_bJet2Pt_bJetbJet_diLept_mumu = fs->make<TH2D>("h1_mindR_bJet2Pt_bJetbJet_diLept_mumu",";min#DeltaR(b jet,b jet);b jet_{2} p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_bJet2Pt_bJetbJet_diLept_elmu = fs->make<TH2D>("h1_mindR_bJet2Pt_bJetbJet_diLept_elmu",";min#DeltaR(b jet,b jet);b jet_{2} p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_bJet2Pt_bJetbJet_semiLept_el = fs->make<TH2D>("h1_mindR_bJet2Pt_bJetbJet_semiLept_el",";min#DeltaR(b jet,b jet);b jet_{2} p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_bJet2Pt_bJetbJet_semiLept_mu = fs->make<TH2D>("h1_mindR_bJet2Pt_bJetbJet_semiLept_mu",";min#DeltaR(b jet,b jet);b jet_{2} p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);

  h2_mindR_FatJetPt_bJetbJet_diLept_elel = fs->make<TH2D>("h1_mindR_FatJetPt_bJetbJet_diLept_elel",";min#DeltaR(b jet,b jet);Fat jet p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_FatJetPt_bJetbJet_diLept_mumu = fs->make<TH2D>("h1_mindR_FatJetPt_bJetbJet_diLept_mumu",";min#DeltaR(b jet,b jet);Fat jet p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_FatJetPt_bJetbJet_diLept_elmu = fs->make<TH2D>("h1_mindR_FatJetPt_bJetbJet_diLept_elmu",";min#DeltaR(b jet,b jet);Fat jet p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_FatJetPt_bJetbJet_semiLept_el = fs->make<TH2D>("h1_mindR_FatJetPt_bJetbJet_semiLept_el",";min#DeltaR(b jet,b jet);Fat jet p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
  h2_mindR_FatJetPt_bJetbJet_semiLept_mu = fs->make<TH2D>("h1_mindR_FatJetPt_bJetbJet_semiLept_mu",";min#DeltaR(b jet,b jet);Fat jet p_{T} [GeV];Number of entries",100,0.,5.,240,0.,600.);
}


ttbarAnalyzer::~ttbarAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ttbarAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  double noCuts=0, leptonCountCuts=1, diLeptInvMassCuts=2, diLeptZVeto=3, diLeptJetCuts=4, diLeptMETCut=5, diLeptBtagCut=6, semiLeptJetCuts=2, semiLeptBtagCut=3;

  h1_CutFlow_diLept_elel->Fill(noCuts);
  h1_CutFlow_diLept_mumu->Fill(noCuts);
  h1_CutFlow_diLept_elmu->Fill(noCuts);
  h1_CutFlow_semiLept_el->Fill(noCuts);
  h1_CutFlow_semiLept_mu->Fill(noCuts);

  bool diElec = false;
  bool diMuon = false;
  bool diElMu = false;
  bool diLept = false;

  bool semiEl = false;
  bool semiMu = false;
  bool semiLept = false;

  std::vector<const reco::GenParticle*> elec_diLept, muon_diLept, lept_diLept, elec_semiLept, muon_semiLept;
  std::vector<const reco::GenJet*> jet_diLept, jet_diLept_btagged, jet_semiLept, jet_semiLept_btagged;
  std::map<const reco::GenJet*, int> jet_flavor;
  std::map<const reco::GenJet*, bool> jet_btag_L, jet_btag_M;

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(genParticleTag, genParticles);

  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByLabel(genJetTag, genJets);

  edm::Handle<reco::JetFlavourMatchingCollection> jetFlavor;
  iEvent.getByLabel(jetFlavorTag, jetFlavor);

  edm::Handle<reco::GenMETCollection> mets;
  iEvent.getByLabel(genMetTag, mets);


  // loop over GenParticles and select electrons and muons
  for( reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it )
  {

    if( !( abs(it->pdgId())==11 || abs(it->pdgId())==13 ) ) continue; // skip if not electron or muon


    if( abs(it->pdgId())==11 && it->status()==1 ) // select stable (status=1) electrons
    {
      if( fabs(it->eta())<elecAbsEtaMaxDiLept && it->pt()>elecPtMinDiLept ) // apply di-leptonic electron selection
        elec_diLept.push_back( &(*it) );

      if( fabs(it->eta())<elecAbsEtaMaxSemiLept && it->pt()>elecPtMinSemiLept ) // apply semi-leptonic electron selection
        elec_semiLept.push_back( &(*it) );
    }


    if( abs(it->pdgId())==13 && it->status()==1 ) // select stable (status=1) muons
    {
      if( fabs(it->eta())<muonAbsEtaMaxDiLept && it->pt()>muonPtMinDiLept ) // apply di-leptonic muon selection
        muon_diLept.push_back( &(*it) );

      if( fabs(it->eta())<muonAbsEtaMaxSemiLept && it->pt()>muonPtMinSemiLept ) // apply semi-leptonic muon selection
        muon_semiLept.push_back( &(*it) );
    }
  }


  // electron-muon cross cleaning
  if( elec_diLept.size()>0 && muon_diLept.size()>0 )
  {
    for(std::vector<const reco::GenParticle*>::iterator eIt = elec_diLept.begin(); eIt != elec_diLept.end(); )
    {
      bool isCloseToMuon = false;
      for(std::vector<const reco::GenParticle*>::const_iterator mIt = muon_diLept.begin(); mIt != muon_diLept.end(); ++mIt )
      {
        if( reco::deltaR( (*eIt)->p4(), (*mIt)->p4() ) < electronMuonDeltaRDiLept ) // remove electrons within dR<0.1 of a muon
        {
          eIt = elec_diLept.erase(eIt);
          isCloseToMuon = true;
          break;
        }
      }
      if( !isCloseToMuon ) ++eIt;
    }
  }


  // form all possible pairs of oppositely changed leptons
  std::vector<std::pair<const reco::GenParticle*,const reco::GenParticle*> > lepton_pairs;

  for(std::vector<const reco::GenParticle*>::const_iterator eIt = elec_diLept.begin(); eIt != elec_diLept.end(); ++eIt )
  {
    for(std::vector<const reco::GenParticle*>::const_iterator e2It = (eIt+1); e2It != elec_diLept.end(); ++e2It )
    {
      if( (*eIt)->charge()==-(*e2It)->charge() ) lepton_pairs.push_back( std::make_pair(*eIt,*e2It) );
    }

    for(std::vector<const reco::GenParticle*>::const_iterator mIt = muon_diLept.begin(); mIt != muon_diLept.end(); ++mIt )
    {
      if( (*eIt)->charge()==-(*mIt)->charge() ) lepton_pairs.push_back( std::make_pair(*eIt,*mIt) );
    }
  }
  for(std::vector<const reco::GenParticle*>::const_iterator mIt = muon_diLept.begin(); mIt != muon_diLept.end(); ++mIt )
  {
    for(std::vector<const reco::GenParticle*>::const_iterator m2It = (mIt+1); m2It != muon_diLept.end(); ++m2It )
    {
      if( (*mIt)->charge()==-(*m2It)->charge() ) lepton_pairs.push_back( std::make_pair(*mIt,*m2It) );
    }
  }

  // sort lepton pairs by their sum Pt
  if( lepton_pairs.size()>1 ) std::sort(lepton_pairs.begin(), lepton_pairs.end(), ordering_Pt<std::pair<const reco::GenParticle*,const reco::GenParticle*> >());

  // check if event is di-leptonic
  if( lepton_pairs.size()>0 )
  {
    lept_diLept.push_back(lepton_pairs.at(0).first);
    lept_diLept.push_back(lepton_pairs.at(0).second);

    int nElecDiLept = 0;
    if( abs(lept_diLept.at(0)->pdgId())==11 ) ++nElecDiLept;
    if( abs(lept_diLept.at(1)->pdgId())==11 ) ++nElecDiLept;

    if( nElecDiLept==2 )       { h1_CutFlow_diLept_elel->Fill(leptonCountCuts); diElec = true; }
    else if ( nElecDiLept==1 ) { h1_CutFlow_diLept_elmu->Fill(leptonCountCuts); diElMu = true; }
    else if ( nElecDiLept==0 ) { h1_CutFlow_diLept_mumu->Fill(leptonCountCuts); diMuon = true; }
  }

  if ( diElec || diMuon || diElMu ) diLept = true;


  // check if event is semi-leptonic (but also not di-leptonic at the same time)
  if ( !diLept && elec_semiLept.size() == 1 ) { h1_CutFlow_semiLept_el->Fill(leptonCountCuts); semiEl = true; }
  else if ( !diLept && muon_semiLept.size() == 1 ) { h1_CutFlow_semiLept_mu->Fill(leptonCountCuts); semiMu = true; }

  if ( semiEl || semiMu) semiLept = true;


  if( !(diLept || semiLept) ) return; // return if event is neither di-leptonic nor semi-leptonic


  TRandom3 rand( iEvent.id().run() + iEvent.id().event() ); // initialize random number generator
  // loop over jets and establish flavor and b-tag information
  for( reco::GenJetCollection::const_iterator it = genJets->begin(); it != genJets->end(); ++it )
  {
    unsigned idx = (it - genJets->begin());
    edm::RefToBase<reco::Jet> jetRef(edm::Ref<reco::GenJetCollection>(genJets, idx));

    jet_flavor[&(*it)] = (*jetFlavor)[jetRef].getFlavour();

    // throw a die
    double rnd = rand.Uniform(0., 1.);

    bool btagged_L = false, btagged_M = false;
    if( abs(jet_flavor[&(*it)])==5 ) // if b jet
    {
      if( rnd < 0.82 ) btagged_L = true;
      if( rnd < 0.67 ) btagged_M = true;
    }
    else if( abs(jet_flavor[&(*it)])==4 ) // if c jet
    {
      if( rnd < 0.38 ) btagged_L = true;
      if( rnd < 0.17 ) btagged_M = true;
    }
    else
    {
      if( rnd < 0.1 )  btagged_L = true;
      if( rnd < 0.01 ) btagged_M = true;
    }

    jet_btag_L[&(*it)] = btagged_L;
    jet_btag_M[&(*it)] = btagged_M;
  }


  // if event is di-leptonic
  if( diLept )
  {
    bool passDiLeptInvMassCuts = false;
    bool passZVeto = false;

    // di-lepton invariant mass cuts
    if( diElec )
    {
      lept_diLept.push_back(elec_diLept.at(0));
      lept_diLept.push_back(elec_diLept.at(1));

      double diLeptMass = (lept_diLept.at(0)->p4() + lept_diLept.at(1)->p4()).mass();

      if( diLeptMass > diLeptMassMin ) passDiLeptInvMassCuts = true;
      if( diLeptMass < zVetoMin || diLeptMass > zVetoMax ) passZVeto = true;
    }
    else if( diMuon )
    {
      lept_diLept.push_back(muon_diLept.at(0));
      lept_diLept.push_back(muon_diLept.at(1));

      double diLeptMass = (lept_diLept.at(0)->p4() + lept_diLept.at(1)->p4()).mass();

      if( diLeptMass > diLeptMassMin ) passDiLeptInvMassCuts = true;
      if( diLeptMass < zVetoMin || diLeptMass > zVetoMax ) passZVeto = true;
    }
    else if( diElMu )
    {
      lept_diLept.push_back(elec_diLept.at(0));
      lept_diLept.push_back(muon_diLept.at(0));

      double diLeptMass = (lept_diLept.at(0)->p4() + lept_diLept.at(1)->p4()).mass();

      if( diLeptMass > diLeptMassMin ) passDiLeptInvMassCuts = true;
      passZVeto = true; // always true for el-mu channel
    }

    if( passDiLeptInvMassCuts )
    {
      if( diElec )      h1_CutFlow_diLept_elel->Fill(diLeptInvMassCuts);
      else if( diMuon ) h1_CutFlow_diLept_mumu->Fill(diLeptInvMassCuts);
      else if( diElMu ) h1_CutFlow_diLept_elmu->Fill(diLeptInvMassCuts);
    }
    else
      return;

    if( passZVeto )
    {
      if( diElec )      h1_CutFlow_diLept_elel->Fill(diLeptZVeto);
      else if( diMuon ) h1_CutFlow_diLept_mumu->Fill(diLeptZVeto);
      else if( diElMu ) h1_CutFlow_diLept_elmu->Fill(diLeptZVeto);
    }
    else
      return;


    // loop over jets
    for( reco::GenJetCollection::const_iterator it = genJets->begin(); it != genJets->end(); ++it )
    {
      if( it->pt()>jetPtMinDiLept && fabs(it->eta())<jetAbsEtaMaxDiLept )
      {
        if( reco::deltaR( it->p4(), lept_diLept.at(0)->p4() )>jetLeptonDeltaRDiLept &&
            reco::deltaR( it->p4(), lept_diLept.at(1)->p4() )>jetLeptonDeltaRDiLept )
          jet_diLept.push_back( &(*it) );
      }
    }

    if( jet_diLept.size()>=2 )
    {
      if( diElec )      h1_CutFlow_diLept_elel->Fill(diLeptJetCuts);
      else if( diMuon ) h1_CutFlow_diLept_mumu->Fill(diLeptJetCuts);
      else if( diElMu ) h1_CutFlow_diLept_elmu->Fill(diLeptJetCuts);
    }
    else
      return;


    // MET cut
    if( diElec && mets->front().pt()>metCutDiLept )      h1_CutFlow_diLept_elel->Fill(diLeptMETCut);
    else if( diMuon && mets->front().pt()>metCutDiLept ) h1_CutFlow_diLept_mumu->Fill(diLeptMETCut);
    else if( diElMu )                                    h1_CutFlow_diLept_elmu->Fill(diLeptMETCut);
    else return;


    // b-tagged jets
    for( std::vector<const reco::GenJet*>::const_iterator it = jet_diLept.begin(); it != jet_diLept.end(); ++it )
      if( jet_btag_L[*it] ) jet_diLept_btagged.push_back(*it);

    if( jet_diLept_btagged.size()>=1 )
    {
      if( diElec )
      {
        h1_CutFlow_diLept_elel->Fill(diLeptBtagCut);
        h1_bTagMultiplicity_diLept_elel->Fill(jet_diLept_btagged.size());

        std::vector<std::pair<const reco::GenJet*,const reco::GenJet*> > filledbJetJetPairs, filledbJetbJetPairs;
        // loop over b-tagged selected jets
        for( std::vector<const reco::GenJet*>::const_iterator it = jet_diLept_btagged.begin(); it != jet_diLept_btagged.end(); ++it )
        {
          double mindR_bJetJet = 1e6;
          const reco::GenJet* mindR_Jet = 0;
          // loop over all selected jets
          for( std::vector<const reco::GenJet*>::const_iterator jt = jet_diLept.begin(); jt != jet_diLept.end(); ++jt )
          {
            if( *it == *jt ) continue; // skip the b-tagged jet itself
            double tempdR = reco::deltaR( (*it)->p4(), (*jt)->p4() );
            if( tempdR < mindR_bJetJet ) { mindR_bJetJet = tempdR; mindR_Jet = *jt; }
          }
          if( mindR_Jet != 0 )
          {
            bool alreadyFilled = false;
            for( std::vector<std::pair<const reco::GenJet*,const reco::GenJet*> >::const_iterator pIt = filledbJetJetPairs.begin(); pIt != filledbJetJetPairs.end(); ++pIt )
            {
              if( (pIt->first==(*it) && pIt->second==mindR_Jet) || (pIt->first==mindR_Jet && pIt->second==(*it)) ) // if already filled
              {
                alreadyFilled = true;
                break;
              }
            }
            if( !alreadyFilled ) // if not already filled
            {
              filledbJetJetPairs.push_back( std::make_pair( *it, mindR_Jet ) );
              h1_mindR_bJetJet_diLept_elel->Fill(mindR_bJetJet);
              h2_mindR_bJetPt_bJetJet_diLept_elel->Fill(mindR_bJetJet,(*it)->pt());
              h2_mindR_JetPt_bJetJet_diLept_elel->Fill(mindR_bJetJet,mindR_Jet->pt());
              h2_mindR_FatJetPt_bJetJet_diLept_elel->Fill(mindR_bJetJet,((*it)->p4() + mindR_Jet->p4()).pt());
            }
          }

          double mindR_bJetbJet = 1e6;
          const reco::GenJet* mindR_bJet = 0;
          // loop over all b-tagged selected jets
          for( std::vector<const reco::GenJet*>::const_iterator jt = jet_diLept_btagged.begin(); jt != jet_diLept_btagged.end(); ++jt )
          {
            if( *it == *jt ) continue; // skip the b-tagged jet itself
            double tempdR = reco::deltaR( (*it)->p4(), (*jt)->p4() );
            if( tempdR < mindR_bJetbJet ) { mindR_bJetbJet = tempdR; mindR_bJet = *jt; }
          }
          if( mindR_bJet != 0 )
          {
            bool alreadyFilled = false;
            for( std::vector<std::pair<const reco::GenJet*,const reco::GenJet*> >::const_iterator pIt = filledbJetbJetPairs.begin(); pIt != filledbJetbJetPairs.end(); ++pIt )
            {
              if( (pIt->first==(*it) && pIt->second==mindR_bJet) || (pIt->first==mindR_bJet && pIt->second==(*it)) ) // if already filled
              {
                alreadyFilled = true;
                break;
              }
            }
            if( !alreadyFilled ) // if not already filled
            {
              filledbJetbJetPairs.push_back( std::make_pair( *it, mindR_bJet ) );
              h1_mindR_bJetbJet_diLept_elel->Fill(mindR_bJetbJet);
              h2_mindR_bJet1Pt_bJetbJet_diLept_elel->Fill(mindR_bJetbJet,(*it)->pt());
              h2_mindR_bJet2Pt_bJetbJet_diLept_elel->Fill(mindR_bJetbJet,mindR_bJet->pt());
              h2_mindR_FatJetPt_bJetbJet_diLept_elel->Fill(mindR_bJetbJet,((*it)->p4() + mindR_bJet->p4()).pt());
            }
          }
        }
      }
      else if( diMuon )
      {
        h1_CutFlow_diLept_mumu->Fill(diLeptBtagCut);
        h1_bTagMultiplicity_diLept_mumu->Fill(jet_diLept_btagged.size());

        std::vector<std::pair<const reco::GenJet*,const reco::GenJet*> > filledbJetJetPairs, filledbJetbJetPairs;
        // loop over b-tagged selected jets
        for( std::vector<const reco::GenJet*>::const_iterator it = jet_diLept_btagged.begin(); it != jet_diLept_btagged.end(); ++it )
        {
          double mindR_bJetJet = 1e6;
          const reco::GenJet* mindR_Jet = 0;
          // loop over all selected jets
          for( std::vector<const reco::GenJet*>::const_iterator jt = jet_diLept.begin(); jt != jet_diLept.end(); ++jt )
          {
            if( *it == *jt ) continue; // skip the b-tagged jet itself
            double tempdR = reco::deltaR( (*it)->p4(), (*jt)->p4() );
            if( tempdR < mindR_bJetJet ) { mindR_bJetJet = tempdR; mindR_Jet = *jt; }
          }
          if( mindR_Jet != 0 )
          {
            bool alreadyFilled = false;
            for( std::vector<std::pair<const reco::GenJet*,const reco::GenJet*> >::const_iterator pIt = filledbJetJetPairs.begin(); pIt != filledbJetJetPairs.end(); ++pIt )
            {
              if( (pIt->first==(*it) && pIt->second==mindR_Jet) || (pIt->first==mindR_Jet && pIt->second==(*it)) ) // if already filled
              {
                alreadyFilled = true;
                break;
              }
            }
            if( !alreadyFilled ) // if not already filled
            {
              filledbJetJetPairs.push_back( std::make_pair( *it, mindR_Jet ) );
              h1_mindR_bJetJet_diLept_mumu->Fill(mindR_bJetJet);
              h2_mindR_bJetPt_bJetJet_diLept_mumu->Fill(mindR_bJetJet,(*it)->pt());
              h2_mindR_JetPt_bJetJet_diLept_mumu->Fill(mindR_bJetJet,mindR_Jet->pt());
              h2_mindR_FatJetPt_bJetJet_diLept_mumu->Fill(mindR_bJetJet,((*it)->p4() + mindR_Jet->p4()).pt());
            }
          }

          double mindR_bJetbJet = 1e6;
          const reco::GenJet* mindR_bJet = 0;
          // loop over all b-tagged selected jets
          for( std::vector<const reco::GenJet*>::const_iterator jt = jet_diLept_btagged.begin(); jt != jet_diLept_btagged.end(); ++jt )
          {
            if( *it == *jt ) continue; // skip the b-tagged jet itself
            double tempdR = reco::deltaR( (*it)->p4(), (*jt)->p4() );
            if( tempdR < mindR_bJetbJet ) { mindR_bJetbJet = tempdR; mindR_bJet = *jt; }
          }
          if( mindR_bJet != 0 )
          {
            bool alreadyFilled = false;
            for( std::vector<std::pair<const reco::GenJet*,const reco::GenJet*> >::const_iterator pIt = filledbJetbJetPairs.begin(); pIt != filledbJetbJetPairs.end(); ++pIt )
            {
              if( (pIt->first==(*it) && pIt->second==mindR_bJet) || (pIt->first==mindR_bJet && pIt->second==(*it)) ) // if already filled
              {
                alreadyFilled = true;
                break;
              }
            }
            if( !alreadyFilled ) // if not already filled
            {
              filledbJetbJetPairs.push_back( std::make_pair( *it, mindR_bJet ) );
              h1_mindR_bJetbJet_diLept_mumu->Fill(mindR_bJetbJet);
              h2_mindR_bJet1Pt_bJetbJet_diLept_mumu->Fill(mindR_bJetbJet,(*it)->pt());
              h2_mindR_bJet2Pt_bJetbJet_diLept_mumu->Fill(mindR_bJetbJet,mindR_bJet->pt());
              h2_mindR_FatJetPt_bJetbJet_diLept_mumu->Fill(mindR_bJetbJet,((*it)->p4() + mindR_bJet->p4()).pt());
            }
          }
        }
      }
      else if( diElMu )
      {
        h1_CutFlow_diLept_elmu->Fill(diLeptBtagCut);
        h1_bTagMultiplicity_diLept_elmu->Fill(jet_diLept_btagged.size());

        std::vector<std::pair<const reco::GenJet*,const reco::GenJet*> > filledbJetJetPairs, filledbJetbJetPairs;
        // loop over b-tagged selected jets
        for( std::vector<const reco::GenJet*>::const_iterator it = jet_diLept_btagged.begin(); it != jet_diLept_btagged.end(); ++it )
        {
          double mindR_bJetJet = 1e6;
          const reco::GenJet* mindR_Jet = 0;
          // loop over all selected jets
          for( std::vector<const reco::GenJet*>::const_iterator jt = jet_diLept.begin(); jt != jet_diLept.end(); ++jt )
          {
            if( *it == *jt ) continue; // skip the b-tagged jet itself
            double tempdR = reco::deltaR( (*it)->p4(), (*jt)->p4() );
            if( tempdR < mindR_bJetJet ) { mindR_bJetJet = tempdR; mindR_Jet = *jt; }
          }
          if( mindR_Jet != 0 )
          {
            bool alreadyFilled = false;
            for( std::vector<std::pair<const reco::GenJet*,const reco::GenJet*> >::const_iterator pIt = filledbJetJetPairs.begin(); pIt != filledbJetJetPairs.end(); ++pIt )
            {
              if( (pIt->first==(*it) && pIt->second==mindR_Jet) || (pIt->first==mindR_Jet && pIt->second==(*it)) ) // if already filled
              {
                alreadyFilled = true;
                break;
              }
            }
            if( !alreadyFilled ) // if not already filled
            {
              filledbJetJetPairs.push_back( std::make_pair( *it, mindR_Jet ) );
              h1_mindR_bJetJet_diLept_elmu->Fill(mindR_bJetJet);
              h2_mindR_bJetPt_bJetJet_diLept_elmu->Fill(mindR_bJetJet,(*it)->pt());
              h2_mindR_JetPt_bJetJet_diLept_elmu->Fill(mindR_bJetJet,mindR_Jet->pt());
              h2_mindR_FatJetPt_bJetJet_diLept_elmu->Fill(mindR_bJetJet,((*it)->p4() + mindR_Jet->p4()).pt());
            }
          }

          double mindR_bJetbJet = 1e6;
          const reco::GenJet* mindR_bJet = 0;
          // loop over all b-tagged selected jets
          for( std::vector<const reco::GenJet*>::const_iterator jt = jet_diLept_btagged.begin(); jt != jet_diLept_btagged.end(); ++jt )
          {
            if( *it == *jt ) continue; // skip the b-tagged jet itself
            double tempdR = reco::deltaR( (*it)->p4(), (*jt)->p4() );
            if( tempdR < mindR_bJetbJet ) { mindR_bJetbJet = tempdR; mindR_bJet = *jt; }
          }
          if( mindR_bJet != 0 )
          {
            bool alreadyFilled = false;
            for( std::vector<std::pair<const reco::GenJet*,const reco::GenJet*> >::const_iterator pIt = filledbJetbJetPairs.begin(); pIt != filledbJetbJetPairs.end(); ++pIt )
            {
              if( (pIt->first==(*it) && pIt->second==mindR_bJet) || (pIt->first==mindR_bJet && pIt->second==(*it)) ) // if already filled
              {
                alreadyFilled = true;
                break;
              }
            }
            if( !alreadyFilled ) // if not already filled
            {
              filledbJetbJetPairs.push_back( std::make_pair( *it, mindR_bJet ) );
              h1_mindR_bJetbJet_diLept_elmu->Fill(mindR_bJetbJet);
              h2_mindR_bJet1Pt_bJetbJet_diLept_elmu->Fill(mindR_bJetbJet,(*it)->pt());
              h2_mindR_bJet2Pt_bJetbJet_diLept_elmu->Fill(mindR_bJetbJet,mindR_bJet->pt());
              h2_mindR_FatJetPt_bJetbJet_diLept_elmu->Fill(mindR_bJetbJet,((*it)->p4() + mindR_bJet->p4()).pt());
            }
          }
        }
      }
    }
    else
      return;
  }
  // if event is semi-leptonic
  else if( semiLept )
  {
    // loop over jets
    for( reco::GenJetCollection::const_iterator it = genJets->begin(); it != genJets->end(); ++it )
    {
      if( it->pt()>jetPtMinLowSemiLept && fabs(it->eta())<jetAbsEtaMaxSemiLept )
      {
        if( reco::deltaR( it->p4(), ( semiEl ? elec_semiLept.at(0)->p4() : muon_semiLept.at(0)->p4() ) )>jetLeptonDeltaRSemiLept )
          jet_semiLept.push_back( &(*it) );
      }
    }

    if( jet_semiLept.size()>=4 &&
        jet_semiLept.at(0)->pt()>jetPtMinHighSemiLept &&
        jet_semiLept.at(1)->pt()>jetPtMinHighSemiLept &&
        jet_semiLept.at(2)->pt()>jetPtMinLowSemiLept &&
        jet_semiLept.at(3)->pt()>jetPtMinLowSemiLept )
    {
      if( semiEl )      h1_CutFlow_semiLept_el->Fill(semiLeptJetCuts);
      else if( semiMu ) h1_CutFlow_semiLept_mu->Fill(semiLeptJetCuts);
    }
    else
      return;


    // b-tagged jets
    for( std::vector<const reco::GenJet*>::const_iterator it = jet_semiLept.begin(); it != jet_semiLept.end(); ++it )
      if( jet_btag_M[*it] ) jet_semiLept_btagged.push_back(*it);

    if( jet_semiLept_btagged.size()>=1 )
    {
      if( semiEl )
      {
        h1_CutFlow_semiLept_el->Fill(semiLeptBtagCut);
        h1_bTagMultiplicity_semiLept_el->Fill(jet_semiLept_btagged.size());

        std::vector<std::pair<const reco::GenJet*,const reco::GenJet*> > filledbJetJetPairs, filledbJetbJetPairs;
        // loop over b-tagged selected jets
        for( std::vector<const reco::GenJet*>::const_iterator it = jet_semiLept_btagged.begin(); it != jet_semiLept_btagged.end(); ++it )
        {
          double mindR_bJetJet = 1e6;
          const reco::GenJet* mindR_Jet = 0;
          // loop over all selected jets
          for( std::vector<const reco::GenJet*>::const_iterator jt = jet_semiLept.begin(); jt != jet_semiLept.end(); ++jt )
          {
            if( *it == *jt ) continue; // skip the b-tagged jet itself
            double tempdR = reco::deltaR( (*it)->p4(), (*jt)->p4() );
            if( tempdR < mindR_bJetJet ) { mindR_bJetJet = tempdR; mindR_Jet = *jt; }
          }
          if( mindR_Jet != 0 )
          {
            bool alreadyFilled = false;
            for( std::vector<std::pair<const reco::GenJet*,const reco::GenJet*> >::const_iterator pIt = filledbJetJetPairs.begin(); pIt != filledbJetJetPairs.end(); ++pIt )
            {
              if( (pIt->first==(*it) && pIt->second==mindR_Jet) || (pIt->first==mindR_Jet && pIt->second==(*it)) ) // if already filled
              {
                alreadyFilled = true;
                break;
              }
            }
            if( !alreadyFilled ) // if not already filled
            {
              filledbJetJetPairs.push_back( std::make_pair( *it, mindR_Jet ) );
              h1_mindR_bJetJet_semiLept_el->Fill(mindR_bJetJet);
              h2_mindR_bJetPt_bJetJet_semiLept_el->Fill(mindR_bJetJet,(*it)->pt());
              h2_mindR_JetPt_bJetJet_semiLept_el->Fill(mindR_bJetJet,mindR_Jet->pt());
              h2_mindR_FatJetPt_bJetJet_semiLept_el->Fill(mindR_bJetJet,((*it)->p4() + mindR_Jet->p4()).pt());
            }
          }

          double mindR_bJetbJet = 1e6;
          const reco::GenJet* mindR_bJet = 0;
          // loop over all b-tagged selected jets
          for( std::vector<const reco::GenJet*>::const_iterator jt = jet_semiLept_btagged.begin(); jt != jet_semiLept_btagged.end(); ++jt )
          {
            if( *it == *jt ) continue; // skip the b-tagged jet itself
            double tempdR = reco::deltaR( (*it)->p4(), (*jt)->p4() );
            if( tempdR < mindR_bJetbJet ) { mindR_bJetbJet = tempdR; mindR_bJet = *jt; }
          }
          if( mindR_bJet != 0 )
          {
            bool alreadyFilled = false;
            for( std::vector<std::pair<const reco::GenJet*,const reco::GenJet*> >::const_iterator pIt = filledbJetbJetPairs.begin(); pIt != filledbJetbJetPairs.end(); ++pIt )
            {
              if( (pIt->first==(*it) && pIt->second==mindR_bJet) || (pIt->first==mindR_bJet && pIt->second==(*it)) ) // if already filled
              {
                alreadyFilled = true;
                break;
              }
            }
            if( !alreadyFilled ) // if not already filled
            {
              filledbJetbJetPairs.push_back( std::make_pair( *it, mindR_bJet ) );
              h1_mindR_bJetbJet_semiLept_el->Fill(mindR_bJetbJet);
              h2_mindR_bJet1Pt_bJetbJet_semiLept_el->Fill(mindR_bJetbJet,(*it)->pt());
              h2_mindR_bJet2Pt_bJetbJet_semiLept_el->Fill(mindR_bJetbJet,mindR_bJet->pt());
              h2_mindR_FatJetPt_bJetbJet_semiLept_el->Fill(mindR_bJetbJet,((*it)->p4() + mindR_bJet->p4()).pt());
            }
          }
        }
      }
      else if( semiMu )
      {
        h1_CutFlow_semiLept_mu->Fill(semiLeptBtagCut);
        h1_bTagMultiplicity_semiLept_mu->Fill(jet_semiLept_btagged.size());

        std::vector<std::pair<const reco::GenJet*,const reco::GenJet*> > filledbJetJetPairs, filledbJetbJetPairs;
        // loop over b-tagged selected jets
        for( std::vector<const reco::GenJet*>::const_iterator it = jet_semiLept_btagged.begin(); it != jet_semiLept_btagged.end(); ++it )
        {
          double mindR_bJetJet = 1e6;
          const reco::GenJet* mindR_Jet = 0;
          // loop over all selected jets
          for( std::vector<const reco::GenJet*>::const_iterator jt = jet_semiLept.begin(); jt != jet_semiLept.end(); ++jt )
          {
            if( *it == *jt ) continue; // skip the b-tagged jet itself
            double tempdR = reco::deltaR( (*it)->p4(), (*jt)->p4() );
            if( tempdR < mindR_bJetJet ) { mindR_bJetJet = tempdR; mindR_Jet = *jt; }
          }
          if( mindR_Jet != 0 )
          {
            bool alreadyFilled = false;
            for( std::vector<std::pair<const reco::GenJet*,const reco::GenJet*> >::const_iterator pIt = filledbJetJetPairs.begin(); pIt != filledbJetJetPairs.end(); ++pIt )
            {
              if( (pIt->first==(*it) && pIt->second==mindR_Jet) || (pIt->first==mindR_Jet && pIt->second==(*it)) ) // if already filled
              {
                alreadyFilled = true;
                break;
              }
            }
            if( !alreadyFilled ) // if not already filled
            {
              filledbJetJetPairs.push_back( std::make_pair( *it, mindR_Jet ) );
              h1_mindR_bJetJet_semiLept_mu->Fill(mindR_bJetJet);
              h2_mindR_bJetPt_bJetJet_semiLept_mu->Fill(mindR_bJetJet,(*it)->pt());
              h2_mindR_JetPt_bJetJet_semiLept_mu->Fill(mindR_bJetJet,mindR_Jet->pt());
              h2_mindR_FatJetPt_bJetJet_semiLept_mu->Fill(mindR_bJetJet,((*it)->p4() + mindR_Jet->p4()).pt());
            }
          }

          double mindR_bJetbJet = 1e6;
          const reco::GenJet* mindR_bJet = 0;
          // loop over all b-tagged selected jets
          for( std::vector<const reco::GenJet*>::const_iterator jt = jet_semiLept_btagged.begin(); jt != jet_semiLept_btagged.end(); ++jt )
          {
            if( *it == *jt ) continue; // skip the b-tagged jet itself
            double tempdR = reco::deltaR( (*it)->p4(), (*jt)->p4() );
            if( tempdR < mindR_bJetbJet ) { mindR_bJetbJet = tempdR; mindR_bJet = *jt; }
          }
          if( mindR_bJet != 0 )
          {
            bool alreadyFilled = false;
            for( std::vector<std::pair<const reco::GenJet*,const reco::GenJet*> >::const_iterator pIt = filledbJetbJetPairs.begin(); pIt != filledbJetbJetPairs.end(); ++pIt )
            {
              if( (pIt->first==(*it) && pIt->second==mindR_bJet) || (pIt->first==mindR_bJet && pIt->second==(*it)) ) // if already filled
              {
                alreadyFilled = true;
                break;
              }
            }
            if( !alreadyFilled ) // if not already filled
            {
              filledbJetbJetPairs.push_back( std::make_pair( *it, mindR_bJet ) );
              h1_mindR_bJetbJet_semiLept_mu->Fill(mindR_bJetbJet);
              h2_mindR_bJet1Pt_bJetbJet_semiLept_mu->Fill(mindR_bJetbJet,(*it)->pt());
              h2_mindR_bJet2Pt_bJetbJet_semiLept_mu->Fill(mindR_bJetbJet,mindR_bJet->pt());
              h2_mindR_FatJetPt_bJetbJet_semiLept_mu->Fill(mindR_bJetbJet,((*it)->p4() + mindR_bJet->p4()).pt());
            }
          }
        }
      }
    }
    else
      return;
  }
  else
    return;
}



// ------------ method called once each job just before starting event loop  ------------
void
ttbarAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
ttbarAnalyzer::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void
ttbarAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
ttbarAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
ttbarAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
ttbarAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ttbarAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ttbarAnalyzer);
