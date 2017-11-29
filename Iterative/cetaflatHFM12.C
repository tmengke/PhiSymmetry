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
//#include <math.h>
//#include <iostream.h>
//#include "histspec.C"
#include "histStat.C"

void cetaflatHFM12(int nIterN=1, double Ethr1=10, double Ethr2=150) { // for HFM, L and S separately
  // nIterN - number of iterations
  // Ethr1 and Ethr2 - E thresholds within which E is estimated
  
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


    sprintf(ctit,"Run2017B.root");



TFile *fila = new TFile (ctit);
  cout<<"File= "<<ctit<<endl;

  TH1F *hcounter =   new TH1F(*((TH1F*)fila->Get("phaseHF/hcounter")));
  cout<<"Stat= "<<hcounter->GetBinContent(2)<<endl;
  cout<<"E within: "<<Ethr1<<" - "<<Ethr2<<endl;

  TH2F* hLmapP = new TH2F("hLmapP","E L HFM;i#eta;i#phi",13,-41.5,-28.5,36,0,72);
  TH2F* hSmapP = new TH2F("hSmapP","E S HFM;i#eta;i#phi",13,-41.5,-28.5,36,0,72);
  TH2F* hLmapM0 = new TH2F("hLmapM0","E0 L HFM;i#eta;i#phi",13,-41.5,-28.5,36,0,72);
  TH2F* hSmapM0 = new TH2F("hSmapM0","E0 S HFM;i#eta;i#phi",13,-41.5,-28.5,36,0,72);
  TH2F* hLmapPc = new TH2F("hLmapPc","corr L HFM;i#eta;i#phi",13,-41.5,-28.5,36,0,72);
  TH2F* hSmapPc = new TH2F("hSmapPc","corr S HFM;i#eta;i#phi",13,-41.5,-28.5,36,0,72);
  hLmapPc->Sumw2(); hSmapPc->Sumw2();
  //TH1F *hLcorr1D = new TH1F("hLcorr1D","Corr L",300,0.5,2);
  //TH1F *hScorr1D = new TH1F("hScorr1D","Corr S",300,0.5,2);
  TH1F *hLcorr1D = new TH1F("hLcorr1D","Corr L",180,0.7,1.5);
  TH1F *hScorr1D = new TH1F("hScorr1D","Corr S",180,0.7,1.5);
  TH1F *hLcorr1Derr = new TH1F("hLcorr1Derr","Err corr L",300,0.,0.025);
  TH1F *hScorr1Derr = new TH1F("hScorr1Derr","Err corr S",300,0.,0.025);
  TH1F *hLdatP[13][36], *hSdatP[13][36], *hLdatPx[13][36], *hSdatPx[13][36];
  for (int ii=0;ii<13;ii++) for (int jj=0;jj<36;jj++) {
    sprintf(ctit,"hL%d_%d",-29-ii,2*jj+1);
    hLdatP[ii][jj] = new TH1F(ctit,ctit,8000,0,250);
    sprintf(ctit,"hS%d_%d",-29-ii,2*jj+1);
    hSdatP[ii][jj] = new TH1F(ctit,ctit,8000,0,250);
  }
  TH1F *htL = new TH1F("htL","htL",20000,0,7e8/3.);
  TH1F *htS = new TH1F("htS","htS",20000,0,5e8/3.);
  //TH1F *htL = new TH1F("htL","htL",20000,0,4e8/40);
  //TH1F *htS = new TH1F("htS","htS",20000,0,2e8/40);
  //TH1F *hLdatPx[13][36], *hSdatPx[13][36];

  TCanvas *cLx[200],*cSx[200];
  TSpline5 *ttL,*ttS;

  Double_t x,y,rPL,rPS,drPL,drPS,mLE,mSE,ermean,rms;
  Double_t xxL[1000],yyL[1000];
  Double_t xxS[1000],yyS[1000];
  Int_t nELP, nESP, nIter=0;
  Double_t mcorrL,scorrL,mcorrS,scorrS,erLP,erSP,rLP,drLP,rSP,drSP,corrL,corrS,dcorrL,dcorrS;
  double mLEphi[13],mSEphi[13],dmLEphi[13],dmSEphi[13];

  TCanvas *ccxx = new TCanvas("ccxx","ccxx",100,300,900,500);
  ccxx->Divide(2,1);

  // loop over channels, 
  for (int ii=0;ii<13;ii++) {  // loop over HFM ieta (rings)
  //for (int ii=1;ii<2;ii++) {
    int ieta=-ii-29;

    mLE=mSE=0;   // ------------------for initial condition
    int nmLE=0, nmSE=0;
    htL->Reset(); htS->Reset();
    for (int ll=1;ll<=72;ll+=2) { // loop over channels within HFM ieta (ring) to get data
      int iphi=ll;
      if (abs(ieta)>39 && (iphi-1)%4==0) continue;

      if (abs(ieta)==29) continue;// && iphi==67) continue;
      hSmapPc->SetBinContent(13-ii,ll/2+1,1);  //  initial correction set to 1
      hLmapPc->SetBinContent(13-ii,ll/2+1,1);
      hSmapPc->SetBinError(13-ii,ll/2+1,1.e-6);  //  non-zero err required
      hLmapPc->SetBinError(13-ii,ll/2+1,1.e-6);
      sprintf(ctit,"phaseHF/espec/E_%d_%d_1",ieta,iphi);
      hLdatPx[ii][ll/2]  =   new TH1F(*((TH1F*)fila->Get(ctit)));  // read spectrum from root file
      hLdatPx[ii][ll/2]->SetAxisRange(Ethr1,Ethr2);
      rLP = hLdatPx[ii][ll/2]->Integral()*hLdatPx[ii][ll/2]->GetMean();  // channel's total E between thresholds
      if (rLP>0) {
	htL->Fill(rLP);
	mLE += rLP;  // ring's total E between thresholds
	nmLE++;
      }
      sprintf(ctit,"phaseHF/espec/E_%d_%d_2",ieta,iphi);
      hSdatPx[ii][ll/2]  =   new TH1F(*((TH1F*)fila->Get(ctit)));
      hSdatPx[ii][ll/2]->SetAxisRange(Ethr1,Ethr2);
      rSP = hSdatPx[ii][ll/2]->Integral()*hSdatPx[ii][ll/2]->GetMean();
      if (rSP>0) {
	htS->Fill(rSP);
	mSE += rSP;
	nmSE++;
      }
      hLmapM0->SetBinContent(13-ii,ll/2+1,rLP);
      hSmapM0->SetBinContent(13-ii,ll/2+1,rSP);
      hSdatP[ii][ll/2]->SetAxisRange(Ethr1,Ethr2);
    }
    if (nmLE>0) mLE /= nmLE; // ring's <E> 
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
    dmLEphi[ii]=htL->GetRMS();
    dmSEphi[ii]=htS->GetRMS();
    printf("ieta %2d :  <E>L= %8.1f (%6.1f) x %d    <E>S= %8.1f (%6.1f) x %d \n",
	   ieta,mLE,dmLEphi[ii],nmLE,mSE,dmSEphi[ii],nmSE);
    
    for (int jj=1;jj<=72;jj+=2) { // loop over channels within HFM ieta (ring) to get correction
      int iphi=jj;
      if (abs(ieta)>39 && (iphi-1)%4==0) continue;
      if (abs(ieta)==20) continue;
      for (nIter=1;nIter<nIterN;nIter++) { // loop over iterations, L channels
	corrL=hLmapPc->GetBinContent(13-ii,jj/2+1);  // correction factor
	hLdatP[ii][jj/2]->Reset();

	for (int kk=1;kk<=hLdatPx[ii][jj/2]->GetNbinsX();kk++) {
	  xxL[kk-1]=hLdatPx[ii][jj/2]->GetBinCenter(kk);
	  yyL[kk-1]=hLdatPx[ii][jj/2]->GetBinContent(kk);
	}
	ttL = new TSpline5("tt",xxL,yyL,1000,"",10,20);  // spline function over spectrum to get
	// smoother image of spectrum to improve convergence 

	for (int kk=1;kk<=hLdatP[ii][jj/2]->GetNbinsX();kk++) {
	  x=hLdatP[ii][jj/2]->GetBinCenter(kk);
	  y=hLdatP[ii][jj/2]->GetBinContent(kk);
	  hLdatP[ii][jj/2]->Fill(x*corrL,ttL->Eval(x)/8.0);  // smoother image of spectrum,
	  // with E shifted (corrected) by correction factor;
	  // factor 8 due to difference in bin size of spectrum and smoother image of spectrum.
	}
	ttL->Delete();

	hLdatP[ii][jj/2]->SetAxisRange(Ethr1,Ethr2);
	rLP = hLdatP[ii][jj/2]->Integral()*hLdatP[ii][jj/2]->GetMean(); // E of corrected spectrum 
	dcorrL=(rLP-mLE)/mLE;  // estimator of difference between <E> and E of corrected spectrum
	if (rLP>0) drLP=
	      sqrt(pow(hLdatP[ii][jj/2]->GetMeanError()/hLdatP[ii][jj/2]->GetMean(),2)+
		   1.f/hLdatP[ii][jj/2]->Integral()+
		   pow(dcorrL/(1.0+sqrt((float) nIter)),2));
	else drLP=1.e-6;
	if (fabs(dcorrL)>0.001) { 
	  corrL*=1-dcorrL/(1.0+sqrt((float) nIter));  // new correction factor
	  //printf("%2d : %2d / %2d / 1 %7.3f %7.3f\n",nIter,ieta,iphi,dcorrL,corrL);
	  hLmapPc->SetBinContent(13-ii,jj/2+1,corrL);
	  hLmapPc->SetBinError(13-ii,jj/2+1,corrL*drLP);
	  hLmapP->SetBinContent(13-ii,jj/2+1,rLP);
	}
	else {
	  printf("%2d : %2d / %2d / 1 %7.3f %8.4f %8.4f\n",nIter,ieta,iphi,dcorrL,corrL,corrL*drLP);
	  hLmapP->SetBinContent(13-ii,jj/2+1,rLP);
	  hLmapPc->SetBinError(13-ii,jj/2+1,corrL*drLP);
	  break;
	}
	if (nIter==nIterN-1) {
	  printf("%2d : %2d / %2d / 1 %7.3f %8.4f %8.4f\n",nIter,ieta,iphi,dcorrL,corrL,corrL*drLP);
	}
      }

      for (nIter=1;nIter<nIterN;nIter++) { 
	corrS=hSmapPc->GetBinContent(13-ii,jj/2+1);
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
	  hSmapPc->SetBinContent(13-ii,jj/2+1,corrS);
	  hSmapPc->SetBinError(13-ii,jj/2+1,corrS*drSP);
	  hSmapP->SetBinContent(13-ii,jj/2+1,rSP);
	}
	else {
	  printf("%2d : %2d / %2d / 2 %7.3f %8.4f %8.4f\n",nIter,ieta,iphi,dcorrS,corrS,corrS*drSP);
	  hSmapP->SetBinContent(13-ii,jj/2+1,rSP);
	  hSmapPc->SetBinError(13-ii,jj/2+1,corrS*drSP);
	  break;
	}
	if (nIter==nIterN-1) {
	  printf("%2d : %2d / %2d / 2 %7.3f %8.4f %8.4f\n",nIter,ieta,iphi,dcorrS,corrS,corrS*drSP);
	}
      }
    }
  }
  //fila->Close();

  cout<<endl<<"Rings :  "<<endl;
  cout<<"       E L        "<<"E S        "<<"eta     "<<"delta eta"<<endl;
  double xeta[13], weta[13], reta[13];
  for (int i=0;i<13;i++) {
    xeta[i]=-(etaBounds[i+1]+etaBounds[i])/2;
    weta[i]=(etaBounds[i+1]-etaBounds[i]);
    mLEphi[i]=mLEphi[i]*36/weta[i];
    mSEphi[i]=mSEphi[i]*36/weta[i];
    dmLEphi[i]=dmLEphi[i]*36/weta[i];
    dmSEphi[i]=dmSEphi[i]*36/weta[i];
    if (i>10) {  mLEphi[i]/=2; mSEphi[i]/=2; dmLEphi[i]/=2; dmSEphi[i]/=2; }
    reta[i] = mSEphi[i]/mLEphi[i];
    cout<<i<<" :  "<<mLEphi[i]<<"    "<<mSEphi[i]<<"    "<<xeta[i]<<"   "<<weta[i]<<endl;
  }
  TCanvas *cgL = new TCanvas("cgL","cgL",300,300,600,600);
  TGraphErrors *grL = new TGraphErrors(13,xeta,mLEphi,0,dmLEphi);
  grL->SetName("grL");
  grL->SetTitle("HFM L;#eta;E_{Ring} / #Delta#eta_{Ring} ,  GeV");
  grL->SetMinimum(0);
  grL->SetMarkerStyle(20);
  grL->Draw("1+PAl");
  cgL->Print("pictHFplot/etaProfHFML.gif");
  //cgL->Print("pictHFmc/etaProfHFML.gif");
  mSEphi[12]/=2; mSEphi[11]/=2;
  TCanvas *cgS = new TCanvas("cgS","cgS",300,300,600,600);
  TGraphErrors *grS = new TGraphErrors(13,xeta,mSEphi,0,dmSEphi);
  grS->SetName("grS");
  grS->SetTitle("HFM S;#eta;E_{Ring} / #Delta#eta_{Ring} ,  GeV");
  grS->SetMinimum(0);
  grS->SetMarkerStyle(20);
  grS->Draw("1+PAl");
  cgS->Print("pictHFplot/etaProfHFMS.gif");
  //cgS->Print("pictHFmc/etaProfHFMS.gif");
  TCanvas *crg = new TCanvas("crg","crg",300,300,600,600);
  TGraphErrors *rg = new TGraphErrors(13,xeta,reta,0,0);
  rg->SetTitle("HFM;#eta;E(S) / E(L)");
  rg->SetMinimum(0);
  rg->Draw("1+PAl");
  crg->Print("pictHFplot/SoverLetaHFM.gif");
  //crg->Print("pictHFmc/SoverLetaHFM.gif");

  TCanvas *cL0 = new TCanvas("cL0","cL0",0,0,650,600);
  hLmapM0->Draw("colz");
  cL0->Update();
  TCanvas *cS0 = new TCanvas("cS0","cS0",1000,0,650,600);
  hSmapM0->Draw("colz");
  cS0->Update();

  //TFile *histf = new TFile("HFMmc.root","RECREATE");

  FILE *ft1;
  //sprintf(ctit,"corrHFMmc_%d_%d.txt",((int) Ethr1),((int) Ethr2));
  sprintf(ctit,"Run2017B_C/corrHFM.txt",ftit,((int) Ethr1),((int) Ethr2));
  if ((ft1 = fopen(ctit,"w"))==NULL){               // Open new file
    //printf("\nNo file %s open => EXIT\n\n",file);
    return;
  }
  printf("\n\n File '%s' open \n\n",ctit);

  TH1D *hprL[13],*hprS[13],*hprL0[13],*hprS0[13];
  TH1D *hprcL[13],*hprcS[13];
  TCanvas *cpr[13],*ccc[13];
  TLine *lin1 = new TLine(0,1,71,1); lin1->SetLineWidth(1);

  int noff=0;
  for (int ii=0;ii<13;ii++) {
    int ieta=-ii-29;

    sprintf(ctit,"HFMcorr_%d_L",-ii-29);  // draw corrections
    hprcL[ii] = hLmapPc->ProjectionY(ctit,13-ii,13-ii);
    hprcL[ii]->SetTitle(ctit);
    sprintf(ctit,"HFMcorr_%d_S",-ii-29);
    hprcS[ii] = hSmapPc->ProjectionY(ctit,13-ii,13-ii);
    hprcS[ii]->SetTitle(ctit);
    ccc[ii] = new TCanvas(ctit,ctit,800,100,500,900);
    ccc[ii]->Divide(1,2);
    ccc[ii]->cd(1);
    if (abs(ieta)>39) {
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
    sprintf(ctit,"pictHFplot/HFMcorr_%d.gif",ii+29);
    //sprintf(ctit,"pictHFmc/HFMcorr_%d.gif",ii+29);
    ccc[ii]->Update();
    ccc[ii]->Print(ctit);
    //hprcL[ii]->Write();
    //hprcS[ii]->Write();

    sprintf(ctit,"HFM_%d_L",-29-ii);  //  draw E depositions
    hprL0[ii] = hLmapM0->ProjectionY(ctit,13-ii,13-ii);
    sprintf(ctit,"HFM_%d_L;i#phi;GeV;",-29-ii);  //  draw E depositions
    hprL0[ii]->SetTitle(ctit);
    sprintf(ctit,"HFM_L_%d",-29-ii);
    hprL[ii] = hLmapP->ProjectionY(ctit,13-ii,13-ii);
    sprintf(ctit,"HFM_%d_S",-29-ii);
    hprS0[ii] = hSmapM0->ProjectionY(ctit,13-ii,13-ii);
    sprintf(ctit,"HFM_%d_S;i#phi;GeV;",-29-ii);  //  draw E depositions
    hprS0[ii]->SetTitle(ctit);
    sprintf(ctit,"HFM_S_%d",-29-ii);
    hprS[ii] = hSmapP->ProjectionY(ctit,13-ii,13-ii);

    cpr[ii] = new TCanvas(ctit,ctit,800,100,500,900);
    cpr[ii]->Divide(1,2);
    cpr[ii]->cd(1);
    if (abs(ieta)>39) {
      hprL0[ii]->Rebin(2);
      hprL[ii]->Rebin(2);
      hprS0[ii]->Rebin(2);
      hprS[ii]->Rebin(2);
    }
    hprL0[ii]->SetFillColor(3);hprL0[ii]->SetLineColor(3);hprL0[ii]->SetLineWidth(3);
    hprL0[ii]->SetMinimum(0);
    hprL0[ii]->SetTitleOffset(0.9,"X");
    hprL0[ii]->Draw("hist");
    hprL[ii]->Draw("samehist");
    cpr[ii]->cd(2);
    hprS0[ii]->SetMinimum(0);
    hprS0[ii]->SetTitleOffset(0.9,"X");
    hprS0[ii]->SetFillColor(3);hprS0[ii]->SetLineColor(3);hprS0[ii]->SetLineWidth(3);
    hprS0[ii]->Draw("hist");
    hprS[ii]->Draw("samehist");
    sprintf(ctit,"pictHFplot/HFM_%d.gif",ii+29);
    //sprintf(ctit,"pictHFmc/HFM_10_100G_%d.gif",ii+29);
    cpr[ii]->Update();
    cpr[ii]->Print(ctit);
    //hprS0[ii]->Write();
    //hprL0[ii]->Write();

    cout<<"Results : "<<endl;
    for (int jj=1;jj<=72;jj+=2) {
      int ieta=-ii-29;
      int iphi=jj;
      if (abs(ieta)>39 && (iphi-1)%4==0) continue;
      if (abs(ieta)==29) continue;// && iphi==67) continue;
      corrL=hLmapPc->GetBinContent(13-ii,jj/2+1);
      corrS=hSmapPc->GetBinContent(13-ii,jj/2+1);
      dcorrL=hLmapPc->GetBinError(13-ii,jj/2+1);
      dcorrS=hSmapPc->GetBinError(13-ii,jj/2+1);
      // -------------------------------------------------------------- !!!
      //if (abs(ieta)>29 && abs(ieta)<41) 
      hLcorr1D->Fill(corrL); 
      hScorr1D->Fill(corrS);
      hLcorr1Derr->Fill(dcorrL);
      hScorr1Derr->Fill(dcorrS);

      noff++;
      //printf("%2d : %2d / %2d / 1 %9.4f %9.4f\n",noff,ieta,iphi,corrL,dcorrL);
      if ((corrL<0.85)||(corrL>1.15))    fprintf(ft1,"%2d   %2d   1 %9.4f %9.4f\n",ieta,iphi,corrL,dcorrL);
      noff++;
      //printf("%2d : %2d / %2d / 2 %9.4f %9.4f\n",noff,ieta,iphi,corrS,dcorrS);
      if ((corrS<0.85)||(corrS>1.15))   fprintf(ft1,"%2d   %2d   2 %9.4f %9.4f\n",ieta,iphi,corrS,dcorrS);
    }
  }
  fclose(ft1);

  for (int ii=0;ii<13;ii++) for (int jj=1;jj<=72;jj+=2) {
      int ieta=-ii-29;
      int iphi=jj;
      if (abs(ieta)>39 && (iphi-1)%4==0) continue;
      if (abs(ieta)==29) continue;//&& iphi==67) continue;
      corrL=hLmapPc->GetBinContent(13-ii,jj/2+1);
      if (fabs(corrL-1)>0.16) printf("%2d / %2d / 1 %9.4f %9.4f\n",ieta,iphi,corrL,dcorrL);
      corrS=hSmapPc->GetBinContent(13-ii,jj/2+1);
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
  c1corr->Print("pictHFplot/corrHFM.gif");
  //c1corr->Print("pictHFmc/corrM.gif");
  c1corr->Update();

  TCanvas *c1correrr =new TCanvas("c1correrr","c1correrr",30,30,900,500);
  c1correrr->Divide(2,1);
  c1correrr->cd(1);  hLcorr1Derr->Draw("hist");  histStat(hLcorr1Derr,1);
  c1correrr->cd(2);  hScorr1Derr->Draw("hist");  histStat(hScorr1Derr,1);
  c1correrr->Print("pictHFplot/corrHFMerr.gif");
  //c1corr->Print("pictHFmc/corrM.gif");                                                                                                                       
  c1corr->Update();
  
  //fila->Close();
  //histf->Close();

  sprintf(ctit,"HFMo_%s_%d_%d.root",ftit,((int) Ethr1),((int) Ethr2));
  TFile *histf = new TFile(ctit,"RECREATE");
  hLcorr1D->Write(); 
  hScorr1D->Write();  
  hLmapP->Write(); 
  hLmapM0->Write(); 
  hLmapPc->Write(); 
  hSmapP->Write(); 
  hSmapM0->Write(); 
  hSmapPc->Write(); 
  grL->Write();
  grS->Write();
  histf->Close();
}

