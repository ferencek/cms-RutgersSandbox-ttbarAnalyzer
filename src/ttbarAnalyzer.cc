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
// $Id$
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
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TH1D.h"
#include "TH2D.h"

//
// class declaration
//

template <class T> struct ordering_Pt {
    ordering_Pt() {}
    bool operator ()(const T* const& a, const T* const& b) { return a->pt() > b->pt(); }
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
  h1_CutFlow_diLept_mumu = fs->make<TH1D>("h1_CutFlow_diLept_mumu",";Cut number;Number of events",ncutBins,cutmin,cutmax);
  h1_CutFlow_diLept_elmu = fs->make<TH1D>("h1_CutFlow_diLept_elmu",";Cut number;Number of events",ncutBins,cutmin,cutmax);
  h1_CutFlow_semiLept_el = fs->make<TH1D>("h1_CutFlow_semiLept_el",";Cut number;Number of events",ncutBins,cutmin,cutmax);
  h1_CutFlow_semiLept_mu = fs->make<TH1D>("h1_CutFlow_semiLept_mu",";Cut number;Number of events",ncutBins,cutmin,cutmax);
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
  double noCuts=0, leptonCuts=1, diLeptInvMassCuts=2, diLeptJetCuts=3, diLeptMETCut=4, semiLeptJetCuts=2;

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

  int nElecDiLept=0, nMuonDiLept=0, nElecSemiLept=0, nMuonSemiLept=0;
  std::vector<const reco::GenParticle*> elec_diLept, muon_diLept, elec_semiLept, muon_semiLept;
  std::vector<const reco::GenJet*> jet_diLept, jet_semiLept;

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(genParticleTag, genParticles);

  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByLabel(genJetTag, genJets);

  edm::Handle<reco::GenMETCollection> mets;
  iEvent.getByLabel(genMetTag, mets);


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


  // count leptons
  nElecDiLept = elec_diLept.size();
  nMuonDiLept = muon_diLept.size();
  nElecSemiLept = elec_semiLept.size();
  nMuonSemiLept = muon_semiLept.size();


  // sort leptons by Pt (if more than 1 lepton found)
  if( nElecDiLept>1 )   std::sort(elec_diLept.begin(), elec_diLept.end(), ordering_Pt<reco::GenParticle>());
  if( nMuonDiLept>1 )   std::sort(muon_diLept.begin(), muon_diLept.end(), ordering_Pt<reco::GenParticle>());
  //if( nElecSemiLept>1 ) std::sort(elec_semiLept.begin(), elec_semiLept.end(), ordering_Pt<reco::GenParticle>());
  //if( nMuonSemiLept>1 ) std::sort(muon_semiLept.begin(), muon_semiLept.end(), ordering_Pt<reco::GenParticle>());


  // check if event is di-leptonic (here we are ignoring lepton charge)
  if ( nElecDiLept >= 2 && nMuonDiLept == 0 ) { h1_CutFlow_diLept_elel->Fill(leptonCuts); diElec = true; }
  else if ( nMuonDiLept >= 2 && nElecDiLept == 0 ) { h1_CutFlow_diLept_mumu->Fill(leptonCuts); diMuon = true; }
  else if ( nElecDiLept>=1 && nMuonDiLept>=1 ) // this is a more complicated case and we need to select the two leading leptons
  {
    int nElecAboveMuon1 = 0, nElecAboveMuon2 = 0;
    for(std::vector<const reco::GenParticle*>::const_iterator it = elec_diLept.begin(); it != elec_diLept.end(); ++it )
    {
      if( (*it)->pt() > muon_diLept.at(0)->pt() ) nElecAboveMuon1++;
      if( muon_diLept.size()>1 )
        if( (*it)->pt() > muon_diLept.at(1)->pt() ) nElecAboveMuon2++;
    }

    if( nElecAboveMuon1>=2 ) { h1_CutFlow_diLept_elel->Fill(leptonCuts); diElec = true; }
    else if ( nElecAboveMuon1==1 ) { h1_CutFlow_diLept_elmu->Fill(leptonCuts); diElMu = true; }
    else if ( nElecAboveMuon1==0 && muon_diLept.size()==1 ) { h1_CutFlow_diLept_elmu->Fill(leptonCuts); diElMu = true; }
    else if ( nElecAboveMuon1==0 && muon_diLept.size()>1 && nElecAboveMuon2==0 ) { h1_CutFlow_diLept_mumu->Fill(leptonCuts); diMuon = true; }
    else if ( nElecAboveMuon1==0 && muon_diLept.size()>1 && nElecAboveMuon2>0 ) { h1_CutFlow_diLept_elmu->Fill(leptonCuts); diElMu = true; }
  }

  if ( diElec || diMuon || diElMu ) diLept = true;


  // check if event is semi-leptonic (but also not di-leptonic at the same time)
  if ( !diLept && nElecSemiLept == 1 ) { h1_CutFlow_semiLept_el->Fill(leptonCuts); semiEl = true; }
  else if ( !diLept && nMuonSemiLept == 1 ) { h1_CutFlow_semiLept_mu->Fill(leptonCuts); semiMu = true; }

  if ( semiEl || semiMu) semiLept = true;


  if( !(diLept || semiLept) ) return; // return if event is neither di-leptonic nor semi-leptonic


  // if event is di-leptonic
  if( diLept )
  {
    std::vector<const reco::GenParticle*> lept_diLept;
    bool passDiLeptInvMassCuts = false;

    // di-lepton invariant mass cuts
    if( diElec )
    {
      lept_diLept.push_back(elec_diLept.at(0));
      lept_diLept.push_back(elec_diLept.at(1));

      double diLeptMass = (lept_diLept.at(0)->p4() + lept_diLept.at(1)->p4()).mass();

      if( diLeptMass>diLeptMassMin && (diLeptMass<zVetoMin || diLeptMass>zVetoMax) ) passDiLeptInvMassCuts = true;
    }
    else if( diMuon )
    {
      lept_diLept.push_back(muon_diLept.at(0));
      lept_diLept.push_back(muon_diLept.at(1));

      double diLeptMass = (lept_diLept.at(0)->p4() + lept_diLept.at(1)->p4()).mass();

      if( diLeptMass>diLeptMassMin && (diLeptMass<zVetoMin || diLeptMass>zVetoMax) ) passDiLeptInvMassCuts = true;
    }
    else if( diElMu )
    {
      lept_diLept.push_back(elec_diLept.at(0));
      lept_diLept.push_back(muon_diLept.at(0));

      double diLeptMass = (lept_diLept.at(0)->p4() + lept_diLept.at(1)->p4()).mass();

      if( diLeptMass>diLeptMassMin ) passDiLeptInvMassCuts = true;
    }

    if( passDiLeptInvMassCuts )
    {
      if( diElec )      h1_CutFlow_diLept_elel->Fill(diLeptInvMassCuts);
      else if( diMuon ) h1_CutFlow_diLept_mumu->Fill(diLeptInvMassCuts);
      else if( diElMu ) h1_CutFlow_diLept_elmu->Fill(diLeptInvMassCuts);
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
    if( diElec && mets->front().pt()>metCutDiLept )     h1_CutFlow_diLept_elel->Fill(diLeptMETCut);
    else if( diMuon && mets->front().pt()>metCutDiLept) h1_CutFlow_diLept_mumu->Fill(diLeptMETCut);
    else if( diElMu )                                   h1_CutFlow_diLept_elmu->Fill(diLeptMETCut);
  }
  // if event is semi-leptonic
  if( semiLept )
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
