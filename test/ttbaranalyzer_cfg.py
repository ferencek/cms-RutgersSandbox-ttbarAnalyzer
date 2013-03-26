import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import copy

options = VarParsing ('python')

options.register('reportEvery',
    1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=1)"
)
options.register('wantSummary',
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Print out trigger and timing summary"
)
## 'maxEvents' is already registered by the Framework, just changing default value
options.setDefault('maxEvents', 10)

options.parseArguments()

process = cms.Process("USER")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(options.wantSummary) )

## Input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       '/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/04E4CE54-CFE1-E111-9801-003048D439C6.root',
       '/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/06035EC2-E3E1-E111-804B-003048C69408.root',
       '/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/FED775BD-B8E1-E111-8ED5-003048C69036.root',
       '/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/0680BB8B-A4E1-E111-958C-003048F02CB2.root',
       '/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/FEC90779-74E1-E111-A432-0025901D493E.root',
       '/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/06C55C32-C7E1-E111-AC3F-0030487F16BF.root',
       '/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/FE4C2F81-D0E1-E111-9080-0030487E0A2D.root',
       '/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/0076C8E3-9AE1-E111-917C-003048D439AA.root',
       '/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/FCE9FDAC-59E1-E111-82BB-0030487E5247.root',
       '/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/02A13705-B0E1-E111-8248-0030487E4EB5.root',
       '/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/FCB7FB42-ACE1-E111-B8AB-0025901D4936.root',
       '/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/02F5A838-8FE1-E111-B0C8-00266CFFA654.root',
       '/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/FC59218E-8AE1-E111-9934-0030487E52A3.root',
    )
)

## Output file
process.TFileService = cms.Service('TFileService',
   fileName = cms.string('output.root')
)

process.myPartons = cms.EDProducer("PartonSelector",
    src = cms.InputTag("genParticles"),
    withLeptons = cms.bool(False)
)

# Flavour byReference
process.AK5byRef = cms.EDProducer("JetPartonMatcher",
    jets = cms.InputTag("ak5GenJets"),
    coneSizeToAssociate = cms.double(0.3),
    partons = cms.InputTag("myPartons")
)

# Flavour byValue AlgoDef
process.AK5byValAlgo = cms.EDProducer("JetFlavourIdentifier",
    srcByReference = cms.InputTag("AK5byRef"),
    physicsDefinition = cms.bool(False),
    leptonInfo = cms.bool(True)
)

## EDAnalyzer
process.analyzer = cms.EDAnalyzer('ttbarAnalyzer',
    GenParticleTag            = cms.InputTag('genParticles'),
    GenJetTag                 = cms.InputTag('ak5GenJets'),
    JetFlavorTag              = cms.InputTag('AK5byValAlgo'),
    GenMetTag                 = cms.InputTag('genMetTrue'),
    ElecPtMinDiLept           = cms.double(20),
    ElecAbsEtaMaxDiLept       = cms.double(2.5),
    ElecPtMinSemiLept         = cms.double(30),
    ElecAbsEtaMaxSemiLept     = cms.double(2.5),
    MuonPtMinDiLept           = cms.double(20),
    MuonAbsEtaMaxDiLept       = cms.double(2.4),
    MuonPtMinSemiLept         = cms.double(26),
    MuonAbsEtaMaxSemiLept     = cms.double(2.1),
    ElectronMuonDeltaRDiLept  = cms.double(0.1),
    DiLeptMassMin             = cms.double(20.),
    ZVetoMin                  = cms.double(76.),
    ZVetoMax                  = cms.double(106.),
    JetPtMinDiLept            = cms.double(30.),
    JetAbsEtaMaxDiLept        = cms.double(2.5),
    JetLeptonDeltaRDiLept     = cms.double(0.5),
    METCutDiLept              = cms.double(40.),
    JetPtMinLowSemiLept       = cms.double(35.),
    JetPtMinHighSemiLept      = cms.double(45.),
    JetAbsEtaMaxSemiLept      = cms.double(2.5),
    JetLeptonDeltaRSemiLept   = cms.double(0.5)
)

process.p = cms.Path(
  process.myPartons
  * process.AK5byRef
  * process.AK5byValAlgo
  * process.analyzer
)
