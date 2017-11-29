#include <vector>
#include <exception>
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "Riostream.h"
#include "TLatex.h"
#include "TColor.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream.h>
#include "histspec.C"
#include "histStat.C"

void cetaflatHFP12(int nIterN=1, double Ethr1=10, double Ethr2=150) {
  
  gStyle->SetOptLogz(0);
  gStyle->SetMarkerSize(0.7);
  gStyle->SetMarkerStyle(20);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetTitleOffset(1.7,"Y");
  gStyle->SetTitleOffset(0.9,"X");
  //gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadLeftMargin(0.18);
  //gStyle->SetNdivisions(516);
  gStyle->SetStatH(0.025);
  gStyle->SetStatW(0.3);
  gStyle->SetTitleW(0.4);
  gStyle->SetTitleX(0.28);
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();

  char ctit[245],ftit[245];
  float etaBounds[14] = {2.853,2.964,3.139,3.314,3.489,3.664,3.839,4.013,4.191,4.363,4.538,4.716,4.889,5.205};

  // ------Histos input: spectra of all channels-----------------------------------

  	sprintf(ctit, "Run2017B.root"); 


 TFile *fila = new TFile (ctit);
  cout<<"File= "<<ctit<<endl;

  TH1F *hcounter =   new TH1F(*((TH1F*)fila->Get("phaseHF/hcounter")));
  cout<<"Stat= "<<hcounter->GetBinContent(2)<<endl;
  cout<<"E within: "<<Ethr1<<" - "<<Ethr2<<endl;

  TH2F* hLmapP = new TH2F("hLmapP","E L HFP;i#eta;i#phi",13,28.5,41.5,36,0,72);
  TH2F* hSmapP = new TH2F("hSmapP","E S HFP;i#eta;i#phi",13,28.5,41.5,36,0,72);
  TH2F* hLmapP0 = new TH2F("hLmapP0","E0 L HFP;i#eta;i#phi",13,28.5,41.5,36,0,72);
  TH2F* hSmapP0 = new TH2F("hSmapP0","E0 S HFP;i#eta;i#phi",13,28.5,41.5,36,0,72);
  TH2F* hLmapPc = new TH2F("hLmapPc","corr L HFP;i#eta;i#phi",13,28.5,41.5,36,0,72);
  TH2F* hSmapPc = new TH2F("hSmapPc","corr S HFP;i#eta;i#phi",13,28.5,41.5,36,0,72);
  hLmapPc->Sumw2(); hSmapPc->Sumw2();
   TH1F *hLcorr1D = new TH1F("hLcorr1Derr","Err corr L",180,0.7,1.5);//300,0.,0.025);
  TH1F *hScorr1D = new TH1F("hScorr1Derr","Err corr S",180,0.7,1.5);//300,0.,0.025);
  TH1F *hLdatP[13][36], *hSdatP[13][36], *hLdatPx[13][36], *hSdatPx[13][36];
  for (int ii=0;ii<13;ii++) for (int jj=0;jj<36;jj++) {
    sprintf(ctit,"hL%d_%d",ii+29,2*jj+1);
    hLdatP[ii][jj] = new TH1F(ctit,ctit,8000,0,250);
    sprintf(ctit,"hS%d_%d",ii+29,2*jj+1);
    hSdatP[ii][jj] = new TH1F(ctit,ctit,8000,0,250);
  }

  TH1F *hEg = new TH1F("hEg","hEg",1000,0,250);
  TH1F *hEb = new TH1F("hEb","hEb",1000,0,250);

  TH1F *htL = new TH1F("htL","htL",20000,0,7e8/3.);
  TH1F *htS = new TH1F("htS","htS",20000,0,5e8/3.);
  TH1F *hLdatPx[13][36], *hSdatPx[13][36];

  TCanvas *cLx[200],*cSx[200];
  TSpline5 *ttL,*ttS;

  Double_t x,y,rPL,rPS,drPL,drPS,mLE,mSE,ermean,rms;
  Double_t xxL[1000],yyL[1000];
  Double_t xxS[1000],yyS[1000];
  Int_t nELP, nESP, nIter=0;
  Double_t mcorrL,scorrL,mcorrS,scorrS,erLP,erSP,rLP,drLP,rSP,drSP,corrL,corrS,dcorrL,dcorrS;
  double mLEphi[13],mSEphi[13], mmmLEphi[13], mmmSEphi[13],dmLEphi[13],dmSEphi[13];
  Int_t ier;
  TCanvas *ccxx = new TCanvas("ccxx","ccxx",100,300,900,500);
  ccxx->Divide(2,1);

  for (int ii=0;ii<13;ii++) {
  //for (int ii=1;ii<2;ii++) {
    int ieta=ii+29;

    mLE=mSE=0;   // ------------------for initial condition
    int nmLE=0, nmSE=0;
    ier = 0;
    htL->Reset(); htS->Reset();
    for (int ll=1;ll<=72;ll+=2) {

      int iphi=ll;
      if ((ieta==30) && (iphi==39)) ier++;
      if ((ieta==34) && (iphi==39)) ier++;
     if (abs(ieta)==40||abs(ieta)==29) continue;
 
      if (abs(ieta)==41 && (iphi-1)%4==0) continue;
       
      hSmapPc->SetBinContent(ii+1,ll/2+1,1);
      hLmapPc->SetBinContent(ii+1,ll/2+1,1);
      hSmapPc->SetBinError(ii+1,ll/2+1,1.e-6);
      hLmapPc->SetBinError(ii+1,ll/2+1,1.e-6);
      
      sprintf(ctit,"phaseHF/espec/E_+%d_%d_1",ieta,iphi);
      
      hLdatPx[ii][ll/2]  =   new TH1F(*((TH1F*)fila->Get(ctit)));
      
      hLdatPx[ii][ll/2]->SetAxisRange(Ethr1,Ethr2);
     
      rLP = hLdatPx[ii][ll/2]->Integral()*hLdatPx[ii][ll/2]->GetMean();
      hLmapP0->SetBinContent(ii+1,ll/2+1,rLP);
      sprintf(ctit,"phaseHF/espec/E_+%d_%d_2",ieta,iphi);
      hSdatPx[ii][ll/2]  =   new TH1F(*((TH1F*)fila->Get(ctit)));
      hSdatPx[ii][ll/2]->SetAxisRange(Ethr1,Ethr2);
      rSP = hSdatPx[ii][ll/2]->Integral()*hSdatPx[ii][ll/2]->GetMean();
      hSmapP0->SetBinContent(ii+1,ll/2+1,rSP);
     
      //if (ieta<=32 && iphi==67) continue;
      if (rLP>0) {
	htL->Fill(rLP);
	mLE += rLP;
	nmLE++;
      }
      if (rSP>0) {
	htS->Fill(rSP);
	mSE += rSP;
	nmSE++;
      }
    }
   
    // }
    if (nmLE>0) mLE /= nmLE; 
    else continue;
    if (nmSE>0) mSE /= nmSE; 
    else continue;
    ccxx->cd(1); htL->Draw("hist");
    ccxx->cd(2); htS->Draw("hist");
    ccxx->Update();
    //histspec(htL,mLE,ermean,rms,4,3);
    //histspec(htS,mSE,ermean,rms,4,3);
    mLEphi[ii]=mLE;
    mSEphi[ii]=mSE;
    mmmLEphi[ii]=mLE;
    mmmSEphi[ii]=mSE;
    dmLEphi[ii]=htL->GetRMS();
    dmSEphi[ii]=htS->GetRMS();
    printf("ieta %2d :  <E>L= %8.1f (%6.1f) x %d    <E>S= %8.1f (%6.1f) x %d \n",
	   ieta,mLE,dmLEphi[ii],nmLE,mSE,dmSEphi[ii],nmSE);
    
    for (int ll=1;ll<=72;ll+=2) {
      if (ieta==40||ieta==29) continue;
      //if (ieta==29) continue;   
	if (abs(ieta)==41 && (ll-1)%4==0) continue;
      int iphi=ll;
      if ((ieta==30) && (iphi==39)) ier++;
      if ((ieta==34) && (iphi==39)) ier++;
      //      if (abs(ieta)>39 && (iphi-1)%4==0) continue;
      //if (ier ==0) {                                                                                                                                                                         

      sprintf(ctit,"phaseHF/espec/E_+%d_%d_1",ieta,iphi);
      hLdatPx[ii][ll/2]  =   new TH1F(*((TH1F*)fila->Get(ctit)));
      hLdatPx[ii][ll/2]->SetAxisRange(Ethr1,Ethr2);
      rLP = hLdatPx[ii][ll/2]->Integral()*hLdatPx[ii][ll/2]->GetMean();
      hLmapP0->SetBinContent(ii+1,ll/2+1,rLP);
      //cout<<mLE<< "!"<<rLP/mLE<<endl;
       sprintf(ctit,"phaseHF/espec/E_+%d_%d_2",ieta,iphi);
      hSdatPx[ii][ll/2]  =   new TH1F(*((TH1F*)fila->Get(ctit)));
      hSdatPx[ii][ll/2]->SetAxisRange(Ethr1,Ethr2);
      rSP = hSdatPx[ii][ll/2]->Integral()*hSdatPx[ii][ll/2]->GetMean();
      hSmapP0->SetBinContent(ii+1,ll/2+1,rSP);
      //      if (iphi == 39) cout<<"ieta = "<<ieta<<"; iphi = "<<iphi<<"; ratio1 = "<<rLP/mLE<<"; ratio2 = "<<rSP/mSE<<endl; 
      if (ieta == 39) cout<<"ieta = "<<ieta<<"  "<<ii<<"  "<<"; iphi = "<<iphi<<"; <E> = "<<hLdatPx[ii][ll/2]->GetMean()<<endl;
    }
    

    // }
    for (int jj=1;jj<=72;jj+=2) {
      int iphi=jj;
      if (abs(ieta)>39 && (iphi-1)%4==0) continue;
      // if (ieta<=32 && iphi==67) {
      //	hLmapP->SetBinContent(ii+1,jj/2+1,hLmapP0->GetBinContent(ii+1,jj/2+1));
      //	hSmapP->SetBinContent(ii+1,jj/2+1,hSmapP0->GetBinContent(ii+1,jj/2+1));
	//	continue;
	// }

      for (nIter=1;nIter<nIterN;nIter++) { //cout<<nIter<<" |  ";
	corrL=hLmapPc->GetBinContent(ii+1,jj/2+1);
	hLdatP[ii][jj/2]->Reset();

	for (int kk=1;kk<=hLdatPx[ii][jj/2]->GetNbinsX();kk++) {
	  xxL[kk-1]=hLdatPx[ii][jj/2]->GetBinCenter(kk);
	  yyL[kk-1]=hLdatPx[ii][jj/2]->GetBinContent(kk);
	}
	ttL = new TSpline5("tt",xxL,yyL,1000,"",10,20);

	for (int kk=1;kk<=hLdatP[ii][jj/2]->GetNbinsX();kk++) {
	  x=hLdatP[ii][jj/2]->GetBinCenter(kk);
	  y=hLdatP[ii][jj/2]->GetBinContent(kk);
	  hLdatP[ii][jj/2]->Fill(x*corrL,ttL->Eval(x)/8.0);
	}
	ttL->Delete();

	hLdatP[ii][jj/2]->SetAxisRange(Ethr1,Ethr2);
	rLP = hLdatP[ii][jj/2]->Integral()*hLdatP[ii][jj/2]->GetMean();
	dcorrL=(rLP-mLE)/mLE;
	if (rLP>0) drLP=
	      sqrt(pow(hLdatP[ii][jj/2]->GetMeanError()/hLdatP[ii][jj/2]->GetMean(),2)+
		   1.f/hLdatP[ii][jj/2]->Integral()+
		   pow(dcorrL/(1.0+sqrt((float) nIter)),2));
	else drLP=1.e-6;
	if (fabs(dcorrL)>0.001) { 
	  corrL*=1-dcorrL/(1.0+sqrt((float) nIter));
	  //printf("%2d : %2d / %2d / 1 %7.3f %7.3f\n",nIter,ieta,iphi,dcorrL,corrL);
	  hLmapPc->SetBinContent(ii+1,jj/2+1,corrL);
	  hLmapPc->SetBinError(ii+1,jj/2+1,corrL*drLP);
	  hLmapP->SetBinContent(ii+1,jj/2+1,rLP);
	}
	else {
	  //  printf("%2d : %2d / %2d / 1 %7.3f %8.4f %8.4f\n",nIter,ieta,iphi,dcorrL,corrL,corrL*drLP);
	  hLmapP->SetBinContent(ii+1,jj/2+1,rLP);
	  hLmapPc->SetBinError(ii+1,jj/2+1,corrL*drLP);
	  break;
	}
	if (nIter==nIterN-1) {
	  //printf("%2d : %2d / %2d / 1 %7.3f %8.4f %8.4f\n",nIter,ieta,iphi,dcorrL,corrL,corrL*drLP);
	}
      }

      for (nIter=1;nIter<nIterN;nIter++) { //cout<<nIter<<" |  ";
	corrS=hSmapPc->GetBinContent(ii+1,jj/2+1);
	hSdatP[ii][jj/2]->Reset();

	for (int kk=1;kk<=hSdatPx[ii][jj/2]->GetNbinsX();kk++) {
	  xxS[kk-1]=hSdatPx[ii][jj/2]->GetBinCenter(kk);
	  yyS[kk-1]=hSdatPx[ii][jj/2]->GetBinContent(kk);
	}
	ttS = new TSpline5("tt",xxS,yyS,1000,"",10,20);

	for (int kk=1;kk<=hSdatP[ii][jj/2]->GetNbinsX();kk++) {
	  x=hSdatP[ii][jj/2]->GetBinCenter(kk);
	  y=hSdatP[ii][jj/2]->GetBinContent(kk);
	  hSdatP[ii][jj/2]->Fill(x*corrS,ttS->Eval(x)/8.0);
	}
	ttS->Delete();

	hSdatP[ii][jj/2]->SetAxisRange(Ethr1,Ethr2);
	rSP = hSdatP[ii][jj/2]->Integral()*hSdatP[ii][jj/2]->GetMean();
	dcorrS=(rSP-mSE)/mSE;
	if (rSP>0) drSP=sqrt(pow(hSdatP[ii][jj/2]->GetMeanError()/hSdatP[ii][jj/2]->GetMean(),2)+
			     1.f/hSdatP[ii][jj/2]->Integral()+
			     pow(dcorrS/(1.0+sqrt((float) nIter)),2));
	else drSP=1.e-6;
	if (fabs(dcorrS)>0.001) { 
	  corrS*=1-dcorrS/(1.0+sqrt((float) nIter));
	  //printf("%2d : %2d / %2d / 1 %7.3f %7.3f\n",nIter,ieta,iphi,dcorrS,corrS);
	  hSmapPc->SetBinContent(ii+1,jj/2+1,corrS);
	  hSmapPc->SetBinError(ii+1,jj/2+1,corrS*drSP);
	  hSmapP->SetBinContent(ii+1,jj/2+1,rSP);
	}
	else {
	  // printf("%2d : %2d / %2d / 2 %7.3f %8.4f %8.4f\n",nIter,ieta,iphi,dcorrS,corrS,corrS*drSP);
	  hSmapP->SetBinContent(ii+1,jj/2+1,rSP);
	  hSmapPc->SetBinError(ii+1,jj/2+1,corrS*drSP);
	  break;
	}
	if (nIter==nIterN-1) {
	  // printf("%2d : %2d / %2d / 2 %7.3f %8.4f %8.4f\n",nIter,ieta,iphi,dcorrS,corrS,corrS*drSP);
	}
      }
    }
  }
  //fila->Close();

  cout<<endl<<"Rings :  "<<endl;
  cout<<"       E L        "<<"E S        "<<"eta     "<<"delta eta"<<endl;
  double xeta[13], weta[13], reta[13];
  for (int i=29;i<=41;i++) {
    xeta[i-29]=(etaBounds[i-28]+etaBounds[i-29])/2;
    weta[i-29]=(etaBounds[i-28]-etaBounds[i-29]);
    mLEphi[i-29]=mLEphi[i-29]*36/weta[i-29];
    mSEphi[i-29]=mSEphi[i-29]*36/weta[i-29];
    dmLEphi[i-29]=dmLEphi[i-29]*36/weta[i-29];
    dmSEphi[i-29]=dmSEphi[i-29]*36/weta[i-29];
    if (i>39) {  mLEphi[i-29]/=2; mSEphi[i-29]/=2; dmLEphi[i-29]/=2; dmSEphi[i-29]/=2; }
    reta[i-29] = mSEphi[i-29]/mLEphi[i-29];
    cout<<i<<" :  "<<mLEphi[i-29]<<"    "<<mSEphi[i-29]<<"    "<<xeta[i-29]<<"   "<<weta[i-29]<<endl;
  }
  TCanvas *cgL = new TCanvas("cgL","cgL",300,300,600,600);
  TGraphErrors *grL = new TGraphErrors(13,xeta,mLEphi,0,dmLEphi);
  grL->SetTitle("HFP L;#eta;E_{Ring} / #Delta#eta_{Ring} ,  GeV");
  grL->SetMinimum(0);
  grL->SetMarkerStyle(20);
  grL->Draw("1+PAl");
  cgL->Print("pictHFplot/etaProfHFPL.gif");
  mSEphi[12]/=2; mSEphi[11]/=2;
  TCanvas *cgS = new TCanvas("cgS","cgS",300,300,600,600);
  TGraphErrors *grS = new TGraphErrors(13,xeta,mSEphi,0,dmSEphi);
  grS->SetTitle("HFP S;#eta;E_{Ring} / #Delta#eta_{Ring} ,  GeV");
  grS->SetMinimum(0);
  grS->SetMarkerStyle(20);
  grS->Draw("1+PAl");
  cgS->Print("pictHFplot/etaProfHFPS.gif");
  TCanvas *crg = new TCanvas("crg","crg",300,300,600,600);
  TGraphErrors *rg = new TGraphErrors(13,xeta,reta,0,0);
  rg->SetTitle("HFP;#eta;E(S) / E(L)");
  rg->SetMinimum(0);
  rg->Draw("1+PAl");
  crg->Print("pictHFplot/SoverLetaHFP.gif");

  TCanvas *cL0 = new TCanvas("cL0","cL0",0,0,650,600);
  hLmapP0->Draw("colz");
  cL0->Update();
  TCanvas *cS = new TCanvas("cS0","cS0",1000,0,650,600);
  hSmapP0->Draw("colz");
  cS0->Update();

  //TFile *histf = new TFile("HFPmc.root","RECREATE");

  FILE *ft1;
  sprintf(ctit,"Run2017B_C/corrHFP.txt",ftit,((int) Ethr1),((int) Ethr2));

  //sprintf(ctit,"corrHFPmc_%d_%d.txt",((int) Ethr1),((int) Ethr2));
  //  sprintf(ctit,"C.txt",ftit,((int) Ethr1),((int) Ethr2));
  if ((ft1 = fopen(ctit,"w"))==NULL){               // Open new file
    printf("\nNo file %s open => EXIT\n\n",file);
    return;
  }
  printf("\n\n File '%s' open \n\n",ctit);

  TH1D *hprL[13],*hprS[13],*hprL0[13],*hprS0[13];
  TH1D *hprcL[13],*hprcS[13];
  TCanvas *cpr[13],*ccc[13];
  TLine *lin1 = new TLine(0,1,71,1); lin1->SetLineWidth(1);

  int noff=0;
  for (int ii=0;ii<13;ii++) {

    sprintf(ctit,"HFPcorr_%d_L",ii+29);  // draw corrections
    hprcL[ii] = hLmapPc->ProjectionY(ctit,ii+1,ii+1);
    hprcL[ii]->SetTitle(ctit);
    sprintf(ctit,"HFPcorr_%d_S",ii+29);
    hprcS[ii] = hSmapPc->ProjectionY(ctit,ii+1,ii+1);
    hprcS[ii]->SetTitle(ctit);
    ccc[ii] = new TCanvas(ctit,ctit,800,100,500,900);
    ccc[ii]->Divide(1,2);
    ccc[ii]->cd(1);
    if (ii+29>39) {
      hprcL[ii]->Rebin(2);
      hprcS[ii]->Rebin(2);
    }
    hprcL[ii]->SetMinimum(0);
    hprcL[ii]->SetTitleOffset(0.9,"X");
    hprcL[ii]->Draw("e");
    lin1->Draw();
    ccc[ii]->cd(2);
    hprcS[ii]->SetMinimum(0);
    hprcS[ii]->SetTitleOffset(0.9,"X");
    hprcS[ii]->Draw("e");
    lin1->Draw();
    sprintf(ctit,"pictHFplot/HFPcorr_%d.gif",ii+29);
    ccc[ii]->Update();
    ccc[ii]->Print(ctit);
    //hprcL[ii]->Write();
    //hprcS[ii]->Write();

    sprintf(ctit,"HFP_%d_L",ii+29);  //  draw E depositions
    hprL0[ii] = hLmapP0->ProjectionY(ctit,ii+1,ii+1);//////mmmLEphi[ii];

    sprintf(ctit,"HFP_%d_L;i#phi;GeV;",29+ii);  //  draw E depositions
    hprL0[ii]->SetTitle(ctit);
    sprintf(ctit,"HFP_L_%d",ii+29);
    hprL[ii] = hLmapP->ProjectionY(ctit,ii+1,ii+1);
    sprintf(ctit,"HFP_%d_S",ii+29);
    hprS0[ii] = hSmapP0->ProjectionY(ctit,ii+1,ii+1);//////mmmSEphi[ii];
    sprintf(ctit,"HFP_%d_S;i#phi;GeV;",29+ii);  //  draw E depositions
    hprS0[ii]->SetTitle(ctit);
    sprintf(ctit,"HFP_S_%d",ii+29);
    hprS[ii] = hSmapP->ProjectionY(ctit,ii+1,ii+1);

    cpr[ii] = new TCanvas(ctit,ctit,800,100,500,900);
    cpr[ii]->Divide(1,2);
    cpr[ii]->cd(1);
    if (ii+29>39) {
      hprL0[ii]->Rebin(2);
      hprL[ii]->Rebin(2);
      hprS0[ii]->Rebin(2);
      hprS[ii]->Rebin(2);
    }
    hprL0[ii]->SetFillColor(3);hprL0[ii]->SetLineColor(3);hprL0[ii]->SetLineWidth(3);
    hprL0[ii]->SetMinimum(0);
    hprL0[ii]->SetTitleOffset(0.9,"X");
    hprL0[ii]->Draw("");
    hprL[ii]->Draw("samehist");//////
    cpr[ii]->cd(2);
    hprS0[ii]->SetMinimum(0);
    hprS0[ii]->SetTitleOffset(0.9,"X");
    hprS0[ii]->SetFillColor(3);hprS0[ii]->SetLineColor(3);hprS0[ii]->SetLineWidth(3);
    hprS0[ii]->Draw("");
    hprS[ii]->Draw("samehist");//////
    sprintf(ctit,"pictHFplot/HFP_%d.gif",ii+29);
    cpr[ii]->Print(ctit);
    //hprS0[ii]->Write();
    //hprL0[ii]->Write();

    cout<<"Results : "<<endl;
    for (int jj=1;jj<=72;jj+=2) {
      int ieta=ii+29;
      int iphi=jj;
      if (abs(ieta)>39 && (iphi-1)%4==0) continue;
      if (ieta==29) continue;
       corrL=hLmapPc->GetBinContent(ii+1,jj/2+1);
      if (fabs(1-corrL)>0.05) {
	cout <<":>5%: eta = "<<ieta<<"; phi = "<<iphi<<endl;
	//cout <<":>5%: eta = "<<ii+1<<"; phi = "<<jj/2+1<<endl;
		hEb->Add(hLdatPx[ii][jj/2]); 
      }
      if (fabs(1-corrL)<0.05) hEg->Add(hLdatPx[ii][jj/2]);

      corrS=hSmapPc->GetBinContent(ii+1,jj/2+1);
      dcorrL=hLmapPc->GetBinError(ii+1,jj/2+1);
      dcorrS=hSmapPc->GetBinError(ii+1,jj/2+1);
      hLcorr1D->Fill(corrL); hScorr1D->Fill(corrS);
      noff++;
      //printf("%2d : %2d / %2d / 1 %9.4f %9.4f\n",noff,ieta,iphi,corrL,dcorrL);
      if ((corrL<0.85)||(corrL>1.15)) fprintf(ft1,"%2d   %2d   1 %9.4f %9.4f\n",ieta,iphi,corrL,dcorrL);
      noff++;
      //printf("%2d : %2d / %2d / 2 %9.4f %9.4f\n",noff,ieta,iphi,corrS,dcorrS);
      if ((corrS<0.85)||(corrS>1.15)) fprintf(ft1,"%2d   %2d   2 %9.4f %9.4f\n",ieta,iphi,corrS,dcorrS);
    }
  }
  fclose(ft1);

  for (int ii=0;ii<13;ii++) for (int jj=1;jj<=72;jj+=2) {
      int ieta=ii+29;
      int iphi=jj;
      if (abs(ieta)>39 && (iphi-1)%4==0) continue;
      if (ieta==29) continue; //&& iphi==67) continue;
      corrL=hLmapPc->GetBinContent(ii+1,jj/2+1);
      if (fabs(corrL-1)>0.16) printf("%2d / %2d / 1 %9.4f %9.4f\n",ieta,iphi,corrL,dcorrL);
      corrS=hSmapPc->GetBinContent(ii+1,jj/2+1);
      if (fabs(corrS-1)>0.16) printf("%2d / %2d / 2 %9.4f %9.4f\n",ieta,iphi,corrS,dcorrS);
  }

  TCanvas *cLcorr =new TCanvas("cLcorr","cLcorr",30,30,600,600);
  cLcorr->SetRightMargin(0.12);
  hLmapPc->SetAxisRange(0.6,1.6,"Z");
  hLmapPc->Draw("colz");
  TCanvas *cScorr =new TCanvas("cScorr","cScorr",30,300,600,600);
  cScorr->SetRightMargin(0.12);
  hSmapPc->SetAxisRange(0.6,1.6,"Z");
  hSmapPc->Draw("colz");

  TCanvas *cL = new TCanvas("cL","cL",0,0,650,600);
  hLmapP->Draw("colz");
  cL->Update();
  TCanvas *cS = new TCanvas("cS","cS",1000,0,650,600);
  hSmapP->Draw("colz");
  cS->Update();

  TCanvas *c1corr =new TCanvas("c1corr","c1corr",30,30,900,500);
  c1corr->Divide(2,1);
  c1corr->cd(1);  hLcorr1D->Draw("hist");  histStat(hLcorr1D,1);
  c1corr->cd(2);  hScorr1D->Draw("hist");  histStat(hScorr1D,1);
  //  c1corr->cd(1);  hEg->Draw("hist");  histStat(hEg,1);                                                                         
  //c1corr->cd(2);  hEb->Draw("hist");  histStat(hEb,1);     




//hLcorr1D->Write(); hScorr1D->Write();  

  c1corr->Print("pictHFplot/corrHFP.gif");
  //c1corr->Print("pictHFmc/corrHFP.gif");
  c1corr->Update();
  
  //fila->Close();
  //histf->Close();

  sprintf(ctit,"HFP_%s_%d_%d.root",ftit,((int) Ethr1),((int) Ethr2));
  TFile *histf = new TFile(ctit,"RECREATE");
  hLcorr1D->Write(); 
  hScorr1D->Write();  
  hLmapP->Write(); 
  hLmapP0->Write(); 
  hLmapPc->Write(); 
  hSmapP->Write(); 
  hSmapP0->Write(); 
  hSmapPc->Write(); 
  grL->Write();
  grS->Write();
  histf->Close();
}

