
import FWCore.ParameterSet.Config as cms

process = cms.Process("PHASEHFX")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32 (30)
)
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.cerr.default.limit = 10

process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")



#--- Global Tag conditions                                                                                                                                                      
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = '92X_dataRun2_Prompt_v9'

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                                      vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                                      minimumNDOF = cms.uint32(4) ,
                                                      maxAbsZ = cms.double(24),
                                                      maxd0 = cms.double(2)
                                                      )

process.phaseHF = cms.EDAnalyzer ("phiSym",
         textFile = cms.untracked.string('yhisto_6.txt'),
         #   gainFile = cms.untracked.string('gain_6.txt'),                                                                                                                     
         #   triggerResultsLabel = cms.InputTag("TriggerResults::HLT"),                                                                                                         
         # corrFile = cms.untracked.string('CorrFactor_rechit_HE_SeptChByCh.txt'),                                                                                              

hfreco =  cms.InputTag("hfreco"),
hbhereco = cms.InputTag("hbhereco")

 )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('yhisto_6.root'),
)

# -------------------------  HBHEHF RECO  

from EventFilter.HcalRawToDigi.HcalRawToDigi_cfi import *
from RecoLocalCalo.HcalRecProducers.HcalHitReconstructor_hbhe_cfi import *
from RecoLocalCalo.HcalRecProducers.HcalHitReconstructor_hf_cfi import *
from RecoLocalCalo.HcalRecProducers.HBHEIsolatedNoiseReflagger_cfi import *
process.newhbheprereco = hbheprereco.clone()
process.newhfreco      = hfreco.clone()
process.newhcalLocalRecoSequence = cms.Sequence(process.newhbheprereco+process.newhfreco)

process.p = cms.Path(
 process.phaseHF
)

process.source = cms.Source ("PoolSource" ,
fileNames=cms.untracked.vstring(

'/store/data/Run2017C/SingleMuon/ALCARECO/HcalCalIterativePhiSym-PromptReco-v1/000/299/368/00000/28A87F56-8A6D-E711-84F0-02163E011A5A.root'

)
)
import FWCore.PythonUtilities.LumiList as LumiList
