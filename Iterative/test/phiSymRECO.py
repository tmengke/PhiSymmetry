import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
process = cms.Process("PHASEHFX",eras.Run2_2018)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32 (1000)
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.cerr.default.limit = 10

process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.GeometryDB_cff")


#--- Global Tag conditions                                                                                                                                                      
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = '101X_upgrade2018_realistic_v3'
process.GlobalTag.globaltag = '101X_dataRun2_Prompt_v9'


process.phaseHF = cms.EDAnalyzer ("phiSym",
                                  textFile = cms.untracked.string('yhisto_6.txt'),
                                  #   gainFile = cms.untracked.string('gain_6.txt'),                                                                                                                     
                                  #   triggerResultsLabel = cms.InputTag("TriggerResults::HLT"),                                                                                                         
                                  # corrFile = cms.untracked.string('CorrFactor_rechit_HE_SeptChByCh.txt'),                                                                                              
                                  
                                  hfreco =  cms.InputTag("hfreco"),
                                  hbhereco = cms.InputTag("hbheprereco")
 )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('phisym.root'),
)

# -------------------------  HBHEHF RECO  

process.load("EventFilter.HcalRawToDigi.HcalRawToDigi_cfi")
process.load("RecoLocalCalo.Configuration.hcalLocalReco_cff")
#from RecoLocalCalo.HcalRecProducers.HBHEIsolatedNoiseReflagger_cfi import *

# Load revenant hcal noise modules
process.load("RecoMET.METProducers.hcalnoiseinfoproducer_cfi")
process.load("CommonTools.RecoAlgos.HBHENoiseFilter_cfi")
process.load("CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi")

# Configuration parameters (especially if a customization is needed):
process.HBHENoiseFilterResultProducer = cms.EDProducer(
    'HBHENoiseFilterResultProducer',
    noiselabel = cms.InputTag('hcalnoise'),
    useBunchSpacingProducer = cms.bool(True), # set to False to disable
    minHPDHits = cms.int32(17),
    minHPDNoOtherHits = cms.int32(10),
    minZeros = cms.int32(9999), #disables the filter
    IgnoreTS4TS5ifJetInLowBVRegion = cms.bool(False),  # set to True if
    defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose"), # set
    minNumIsolatedNoiseChannels = cms.int32(10),
    minIsolatedNoiseSumE = cms.double(50.0),
    minIsolatedNoiseSumEt = cms.double(25.0)
)
# To apply filter decision in CMSSW as an EDFilter:
process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
    inputLabel =
cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
    reverseDecision = cms.bool(False)
)

#----------------------------
# Paths/Sequences Definitions
#----------------------------
process.digiPath = cms.Path(
    process.hcalDigis
)

process.recoPath = cms.Path(
    process.hfprereco
    *process.hfreco
    *process.hbheprereco
)

process.p = cms.Path(
 process.hcalnoise
 *process.HBHENoiseFilterResultProducer
 *process.ApplyBaselineHBHENoiseFilter
 *process.phaseHF
)

process.source = cms.Source ("PoolSource" ,
                             fileNames=cms.untracked.vstring(
			     '/store/data/Run2018A/Commissioning/RAW/v1/000/315/252/00000/6A8B727E-9F49-E811-9014-FA163E8D185B.root'
				)
)
