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
       # Di-leptonic ttbar
       #'/store/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v2/00001/FEAC1583-031B-E211-AF54-00215E2283FA.root',
       #'/store/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v2/00001/FE1884B5-231B-E211-8822-00215E21D57C.root',
       #'/store/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v2/00001/FCFD5C18-2D1B-E211-9512-001A645C1FFC.root',
       #'/store/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v2/00001/FC2A0CAC-0D1B-E211-9017-00215E221EEA.root',
       #'/store/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v2/00001/FAAC928B-091B-E211-8858-00215E21D690.root',
       #'/store/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v2/00001/F8FEE31B-231B-E211-996B-00215E21DD0E.root',
       #'/store/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v2/00001/F8B70B6D-2D1B-E211-88B1-00215E21DAB0.root',
       #'/store/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v2/00001/F8805736-EB1A-E211-AA09-00215E2222DA.root',
       #'/store/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v2/00001/F87CDCEE-ED1A-E211-9AB9-00215E222256.root'
       # Semi-leptonic ttbar
       #'/store/mc/Summer12_DR53X/TTJets_SemiLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A_ext-v1/00001/E8B4F7A1-1025-E211-BB45-002590200B40.root',
       #'/store/mc/Summer12_DR53X/TTJets_SemiLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A_ext-v1/00001/E8324A46-0525-E211-BD58-001E6739751C.root',
       #'/store/mc/Summer12_DR53X/TTJets_SemiLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A_ext-v1/00001/E870615A-2825-E211-A499-002590200B60.root',
       #'/store/mc/Summer12_DR53X/TTJets_SemiLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A_ext-v1/00001/EA0D1319-F624-E211-A6EA-001E67396BAD.root',
       #'/store/mc/Summer12_DR53X/TTJets_SemiLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A_ext-v1/00001/E885EA9C-D924-E211-A1F5-001E673981C4.root',
       #'/store/mc/Summer12_DR53X/TTJets_SemiLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A_ext-v1/00001/EA4BE856-D724-E211-ABBD-002590200B14.root',
       #'/store/mc/Summer12_DR53X/TTJets_SemiLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A_ext-v1/00001/E89F6565-4925-E211-A8B7-001E67397D7D.root',
       #'/store/mc/Summer12_DR53X/TTJets_SemiLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A_ext-v1/00001/EA583F0B-0F25-E211-9369-001E67396C9D.root',
       #'/store/mc/Summer12_DR53X/TTJets_SemiLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A_ext-v1/00001/EA6146BB-EA24-E211-9F28-001E67396577.root',
       #'/store/mc/Summer12_DR53X/TTJets_SemiLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A_ext-v1/00001/EA6B1BAB-9D25-E211-B921-0025902008E4.root'
       # Inclusive ttbar
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

## EDAnalyzer
process.analyzer = cms.EDAnalyzer('ttbarAnalyzer',
    GenParticleTag            = cms.InputTag('genParticles'),
    GenJetTag                 = cms.InputTag('ak5GenJets'),
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

process.p = cms.Path(process.analyzer)
