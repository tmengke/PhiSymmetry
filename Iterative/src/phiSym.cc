// include files                                                                                                                                                        
#include <memory>

// user include files                                                                                                                                                   
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "CondFormats/HcalObjects/interface/HcalQIEShape.h"
#include "CondFormats/HcalObjects/interface/HcalQIECoder.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CondFormats/HcalObjects/interface/HcalPedestals.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//                                            

#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include <TROOT.h>
#include <TSystem.h>
#include "TFile.h"
#include <TCanvas.h>
#include <cmath>
#include "TMath.h"

using namespace cms;
using namespace edm;
using namespace std;
using namespace reco;

//                                                                                                                                                                      
// class declaration                                                                                                                                                    
//                                                                                                                                                                      

class phiSym : public edm::EDAnalyzer {
 public:
  explicit phiSym (const edm::ParameterSet&);
  ~phiSym();

 private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  edm::Service<TFileService> fs;

  std::string histfile;
  std::string textfile;

  TFile* mFile;
  FILE* tFile;

  Int_t runNumb, EventN;

  double e_RecHitHE[14][72], e_RecHitHBp15[72], e_RecHitHBm15[72];
  TH1F *hcounter,*herun, *hlumi, *heventn, *hBX, *hvertex;
  TH1F *hen[26][36][2];
  TH1F *henhbp[16][72][2], *henhbm[16][72][2];
  TH1F *henhep[14][72][3], *henhem[14][72][3];

  InputTag HBHENoiseFilterResultLabel_;
  HcalCalibrations calibs_;

  edm::EDGetTokenT<HFRecHitCollection> mhfreco;
  edm::EDGetTokenT<HBHERecHitCollection> mhbhereco;

};

// constructors and destructor                                                                                                                                                  
//                                                                                                                                                                              
phiSym::phiSym(const edm::ParameterSet& iConfig)
{
  mhfreco  = consumes<HFRecHitCollection>(iConfig.getParameter<edm::InputTag>("hfreco"));//for RECO                                                                             
  mhbhereco = consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("hbhereco"));//for RECO                                                                        
  textfile = iConfig.getUntrackedParameter<string>("textFile");
}

phiSym::~phiSym()
{
}

// member functions                                                                                                                                                             

// ------------ method called to for each event  ------------                                                                                                                   
 void phiSym::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   printf("Starting :\n");
   edm::EventID eventId = iEvent.id();
   int runNumber = eventId.run ();
   int eventNumber = eventId.event ();
   int lumi = eventId.luminosityBlock();

   float eventWeight = 1;

   runNumb=runNumber;
   EventN++;
 cout <<"event number--- "<<EventN<<endl;


   hcounter->Fill(0); // total                                                                                                                                                  

   Int_t nBX=0;// iBX=1, nORBIT=0;                                                                                                                                              
   nBX = iEvent.bunchCrossing();

   hcounter->Fill(1);

   //hf rechits                                                                                                                                                                 
   Handle<HFRecHitCollection> hf_hits_h;
   iEvent.getByToken(mhfreco, hf_hits_h);
   const HFRecHitCollection* hf_hits = hf_hits_h.failedToGet () ? 0 : &*hf_hits_h;

   //hcal rechits                                                                                                                                                               
   Handle<HBHERecHitCollection> hbhe_hits_h;
   iEvent.getByToken(mhbhereco, hbhe_hits_h);
   const HBHERecHitCollection* hbhe_hits = hbhe_hits_h.failedToGet () ? 0 : &*hbhe_hits_h;

   double Etot=0;//, enen = 0;                                                                                                                                                  
   int jeta;//,jphi,jdepth;                                                                                                                                                     


   // ------------ HF -----------                                                                                                                                               
   if (hf_hits_h.isValid()) {
     cout<<"hf_hits...."<<endl;
     for (HFRecHitCollection::const_iterator hfhit=hf_hits->begin(); hfhit!=hf_hits->end(); hfhit++) {
       int ieta = hfhit->id().ieta();
       int iphi = hfhit->id().iphi();
       int depth = hfhit->id().depth();

        double energy = hfhit->energy();

       if (ieta>0) jeta = ieta-29;
       else jeta = 13-ieta-29;
       hen[jeta][(iphi-1)/2][depth-1]->Fill(energy, eventWeight );//****                                                                                                        
       if (energy>10) Etot += energy;
     }
   }
   else printf("No HF RecHits: run= %d  ev= %d :\n",runNumber,eventNumber); // ------------                                                                                     

   // ------------ HBHE -----------                                                                                                                                             


   if (hbhe_hits_h.isValid()) {
     cout<<"hbhe_hits...."<<endl;
     hcounter->Fill(3);

     for (HBHERecHitCollection::const_iterator hbhehit=hbhe_hits->begin(); hbhehit!=hbhe_hits->end(); hbhehit++) {
       int ieta = hbhehit->id().ieta();
       int iphi = hbhehit->id().iphi();
       int depth = hbhehit->id().depth();
        double energy = hbhehit->energy();

       if (abs(ieta)>16 || (abs(ieta)==16 && depth==3)) { // HE                                                                                                                 
         if (ieta>0) henhep[ieta-16][iphi-1][depth-1]->Fill(energy, eventWeight);
         else  henhem[-ieta-16][iphi-1][depth-1]->Fill(energy, eventWeight);
         if (energy>4) Etot += energy;
       }
       else { // HB                                                                                                                                                             

         if (ieta>0) {
           if (ieta!=15) henhbp[ieta-1][iphi-1][depth-1]->Fill(energy, eventWeight);
              else  e_RecHitHBp15[iphi-1]+= energy;
         }
         else if (ieta!=-15) henhbm[-ieta-1][iphi-1][depth-1]->Fill(energy, eventWeight);
         else  e_RecHitHBm15[iphi-1]+= energy;
          if (energy>4) Etot += energy;
       }
     }

       for (int j = 0; j<72; j++){   // iphi                                                                                                                                    
         if (e_RecHitHBp15[j]>0.001)     henhbp[14][j][0]->Fill(e_RecHitHBp15[j],eventWeight);
         if (e_RecHitHBm15[j]>0.001)     henhbm[14][j][0]->Fill(e_RecHitHBm15[j],eventWeight);//}                                                                               

       }
     }

   else printf("No RecHits: run= %d  ev= %d :\n",runNumber,eventNumber); // ------------                                                                                        



   herun->Fill(runNumber);
   heventn->Fill(EventN);
   hlumi->Fill(lumi);
   hBX->Fill(nBX,Etot);
   return;
}

void phiSym::beginJob() {

  EventN=0;
  char htit[64];
  if ((tFile = fopen(textfile.c_str(),"w"))==NULL) printf("\nNo textfile open\n\n");

  herun = fs->make< TH1F>("herun","E(HF) vs Nrun;Nrun;GeV",5000,246000,251000);
  heventn = fs->make< TH1F>("heventn","E(HF) vs EventN;EventN;GeV",100000,0,100000);
  hlumi = fs->make< TH1F>("hlumi","E(HF) vs Lumi;EventN;GeV",100000,0,100000);
  hBX = fs->make< TH1F>("hBX","E(HF) vs nBX;BX;GeV",4096,0,4096);
  hcounter = fs->make< TH1F>("hcounter","hcounter",101,-0.5,100.5);
  hvertex = fs->make< TH1F>("hvertex","hvertex",100, 0, 99.5);

  TFileDirectory ESpec = fs->mkdir( "espec" );
  for (int i=0;i<13;i++) for (int j=0;j<36;j++) for (int k=0;k<2;k++) {
    if (i>10 && j%2==0) continue;
    sprintf(htit,"E_+%d_%d_%d",i+29,j*2+1,k+1);
    hen[i][j][k] = ESpec.make< TH1F>(htit,htit,1000,0,250); // E Rec                                                                                                            
    sprintf(htit,"E_-%d_%d_%d",i+29,j*2+1,k+1);
    hen[i+13][j][k] = ESpec.make< TH1F>(htit,htit,1000,0,250);
  }

  TFileDirectory EBSpec = fs->mkdir( "eHBspec" );
  for (int i=0;i<16;i++) for (int j=0;j<72;j++) for (int k=0;k<2;k++) {
    if (i+1<15 && k==1) continue;
    sprintf(htit,"E_+%d_%d_%d",i+1,j+1,k+1);
    henhbp[i][j][k] = EBSpec.make< TH1F>(htit,htit,1000,0,250); // E Rec                                                                                                        
    sprintf(htit,"E_-%d_%d_%d",i+1,j+1,k+1);
    henhbm[i][j][k] = EBSpec.make< TH1F>(htit,htit,1000,0,250);
  }

  TFileDirectory EESpec = fs->mkdir( "eHEspec" );
  for (int i=0;i<14;i++) for (int j=0;j<72;j++) for (int k=0;k<3;k++) {
    if (i+16==16 && k<2) continue;
    if (i+16==17 && k>0) continue;
    if (i+16>17 && i+16<27 && k==2) continue;
    if (i+16==29 && k==2) continue;
    if (i+16>20 && (j+1)%2==0) continue;
    sprintf(htit,"E_+%d_%d_%d",i+16,j+1,k+1);
    henhep[i][j][k] = EESpec.make< TH1F>(htit,htit,1000,0,250); // E Rec                                                                                                        
    sprintf(htit,"E_-%d_%d_%d",i+16,j+1,k+1);
    henhem[i][j][k] = EESpec.make< TH1F>(htit,htit,1000,0,250);
  }

  std::cout<<std::endl<<"beginJob: histfile="<<histfile.c_str()<<"  textfile="<<textfile.c_str()<<std::endl;
  return;

}

void phiSym::endJob() {

  fprintf(tFile,"#RunN %d   Events processed %d\n",runNumb,EventN);
  std::cout<<"endJob: histos processing..."<<std::endl;
  std::cout<<"RunN= "<<runNumb<<"  Events processed= "<<EventN<<std::endl;

  fclose(tFile);
  std::cout<<std::endl<<" --endJob-- done"<<std::endl;
  return;
}

//define this as a plug-in                                                                                                                                                      
DEFINE_FWK_MODULE(phiSym);


