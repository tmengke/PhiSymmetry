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
#include "skipHBChannel.C"

void cetaflatHBMt(int nIterN=20, double Ethr1=4, double Ethr2=100) {
  
  static const double theHBHEEtaBounds[] = { 0.000, 0.087, 0.087*2, 0.087*3, 0.087*4, 
					     0.087*5, 0.087*6, 0.087*7, 0.087*8, 0.087*9,   
					     0.087*10, 0.087*11, 0.087*12, 0.087*13, 0.087*14,   
					     0.087*15, 0.087*16, 0.087*17, 0.087*18, 0.087*19,   
					     1.74, 1.83, 1.93, 2.043, 2.172, 
					     2.332, 2.5, 2.65, 2.868, 3.000 };
 
  static const double theHFEtaBounds[] = { 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 
					   4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191 };


  gStyle->SetOptLogz(0);
  gStyle->SetMarkerSize(0.7);
  gStyle->SetMarkerStyle(20);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetTitleOffset(1.7,"Y");
  gStyle->SetTitleOffset(0.9,"X");
  //gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadLeftMargin(0.2);
  //gStyle->SetNdivisions(516);
  gStyle->SetStatH(0.09);
  gStyle->SetStatW(0.3);
  gStyle->SetTitleW(0.4);
  gStyle->SetTitleX(0.25);
  gStyle->SetOptStat(0);
  //gStyle->SetOptStat(1111111);
  gROOT->ForceStyle();

  double ahmax;
  TLatex tt1, tt2, tt3, tt4, tt5;
  tt3.SetTextFont(42);
  tt3.SetTextAlign(11);
  tt3.SetTextColor(4);
  tt3.SetTextSize(0.045);
  tt3.SetTextAngle(0);
  tt2.SetTextFont(42);
  tt2.SetTextAlign(11);
  tt2.SetTextColor(2);
  tt2.SetTextSize(0.045);
  tt2.SetTextAngle(0);
  tt1.SetTextFont(42);
  tt1.SetTextAlign(11);
  tt1.SetTextColor(1);
  tt1.SetTextSize(0.045);
  tt1.SetTextAngle(0);
  tt4.SetTextFont(42);
  tt4.SetTextAlign(11);
  tt4.SetTextColor(6);
  tt4.SetTextSize(0.045);
  tt4.SetTextAngle(0);
  tt5.SetTextFont(42);
  tt5.SetTextAlign(11);
  tt5.SetTextColor(3);
  tt5.SetTextSize(0.045);
  tt5.SetTextAngle(0);

  char ctit[145];

  // ---------------- Histos input --------------------------------------

  char ftit[145];
  sprintf(ftit,"%s","phiSymCorr_2017");

  sprintf(ctit,"Run2017B.root");
  
  TFile *fila = new TFile (ctit);
  cout<<"File= "<<ctit<<endl;

  TH1F *hcounter =   new TH1F(*((TH1F*)fila->Get("phaseHF/hcounter")));
  cout<<"Stat= "<<hcounter->GetBinContent(2)<<endl;
  cout<<"E within: "<<Ethr1<<" - "<<Ethr2<<endl;

  hprof  = new TProfile("hprof","Profile",72,0,72,0.5,1.5);

  TH2F* hmapP = new TH2F("hmapP","E  HBM;i#eta;i#phi",16,-16.5,-0.5,72,0,72);
  TH2F* hmapP0 = new TH2F("hmapP0","E0  HBM;i#eta;i#phi",16,-16.5,-0.5,72,0,72);
  TH2F* hmapPc = new TH2F("hmapPc","corr  HBM;i#eta;i#phi",16,-16.5,-0.5,72,0,72);
  hmapPc->Sumw2();
  TH1F *hcorr1D_1 = new TH1F("hcorr1D_1","Corr",300,0.,2.);
  TH1F *hcorr1D_5 = new TH1F("hcorr1D_5","Corr",300,0.,2.);
  TH1F *hcorr1D_11 = new TH1F("hcorr1D_11","Corr",300,0.,2.);
  TH1F *hcorr1Derr_1 = new TH1F("hcorr1Derr_1","Corr",300,0.,0.04);
  TH1F *hcorr1Derr_5 = new TH1F("hcorr1Derr_5","Corr",300,0.,2.);
  TH1F *hcorr1Derr_11 = new TH1F("hcorr1Derr_11","Corr",300,0.,2.);

  TH1F *ht = new TH1F("ht","ht",20000,0,1e6);
  TH1F *htx = new TH1F("htx","htx",20000,0,1e4);
  //TH1F *ht = new TH1F("ht","ht",40000,0,7e4);
  //TH1F *ht = new TH1F("ht","ht",40000,0,7e3);
  TH1F *htr = new TH1F("htr","htr",5000,0,3);

  TH1F *hdatP[16][72], *hdatPx[16][72];
  for (int ii=0;ii<16;ii++) for (int jj=0;jj<72;jj++) {
    sprintf(ctit,"h%d_%d",ii+1,jj+1);
    hdatP[ii][jj] = new TH1F(ctit,ctit,10000,0,250);
  }

  TCanvas *cx[400];
  TSpline5 *tt;

  Double_t x,y,rPL,rPS,mLE,mSE,ermean,rms;
  Double_t xx[4000],yy[4000];
  Int_t nELP, nESP, nIter=0;
  Double_t mcorrL,scorrL,mcorrS,scorrS,erLP,erSP,rLP,drLP,rSP,corrL,corrS,dcorrL,dcorrS;
  double mLEphi[16];
  int ieta;
  TCanvas *ccxx = new TCanvas("ccxx","ccxx",0,400,700,400);
  ccxx->Divide(2,1);
  for (int ii=0;ii<16;ii++) {
    ieta=-ii-1;

    mLE=mSE=0;   // ------------------for initial condition
    int nmLE=0, nmSE=0;
    ht->Reset(); htx->Reset();
    for (int ll=0;ll<72;ll++) {
      int iphi=ll+1;
      hmapPc->SetBinContent(16-ii,ll+1,1);
      hmapPc->SetBinError(16-ii,ll+1,1.e-6);
      sprintf(ctit,"phaseHF/eHBspec/E_-%d_%d_1",-ieta,iphi);
      hdatPx[ii][ll]  =   new TH1F(*((TH1F*)fila->Get(ctit)));
      hdatPx[ii][ll]->SetAxisRange(Ethr1,Ethr2);
      rLP = hdatPx[ii][ll]->Integral()*hdatPx[ii][ll]->GetMean();
      hmapP0->SetBinContent(16-ii,ll+1,rLP);
      if (skipHBChannel(iphi,ieta)) continue;
      if (rLP>0) {
	ht->Fill(rLP); htx->Fill(rLP);
	mLE += rLP;
	nmLE++;
	drLP=rLP*sqrt(pow(1./hdatPx[ii][ll]->Integral(),2)+
		      pow(hdatPx[ii][ll]->GetMeanError()/hdatPx[ii][ll]->GetMean(),2));
	hmapP0->SetBinError(16-ii,ll+1,drLP);
      }
      else hmapP0->SetBinError(16-ii,ll+1,0);
    }
    if (nmLE>0) mLE /= nmLE; 
    else mLE=0;
    ccxx->cd(1); ht->Draw("hist");
    ccxx->cd(2); htx->Draw("hist");
    ccxx->Update();
    // if (htx->GetBinContent(20001)>1) histspec(ht,mLE,ermean,rms,4,-5);
    //else histspec(htx,mLE,ermean,rms,4,-5);
    mLEphi[ii]=mLE;
    printf("ieta %2d :  <E>= %8.1f \n",ieta,mLE);
    if (ht->GetMean()>0) htr->Fill(ht->GetRMS()/ht->GetMean());
    
    for (int jj=0;jj<72;jj++) {
      int iphi=jj+1;
      if (skipHBChannel(iphi,ieta)) {
      hmapP->SetBinContent(16-ii,jj+1,hmapP0->GetBinContent(16-ii,jj+1,rLP));
	 	continue;
       }

      for (nIter=1;nIter<nIterN;nIter++) { //cout<<nIter<<" |  ";
	if (hmapP0->GetBinContent(16-ii,jj+1)<=0) continue;
	corrL=hmapPc->GetBinContent(16-ii,jj+1);
	hdatP[ii][jj]->Reset();

	for (int kk=1;kk<=hdatPx[ii][jj]->GetNbinsX();kk++) {
	  xx[kk-1]=hdatPx[ii][jj]->GetBinCenter(kk);
	  yy[kk-1]=hdatPx[ii][jj]->GetBinContent(kk);
	}
	tt = new TSpline5("tt",xx,yy,1000,"",10,20);

	for (int kk=1;kk<=hdatP[ii][jj]->GetNbinsX();kk++) {
	  x=hdatP[ii][jj]->GetBinCenter(kk);
	  y=hdatP[ii][jj]->GetBinContent(kk);
	  hdatP[ii][jj]->Fill(x*corrL,tt->Eval(x)/10.0);
	}
	tt->Delete();

	hdatP[ii][jj]->SetAxisRange(Ethr1,Ethr2);
	rLP = hdatP[ii][jj]->Integral()*hdatP[ii][jj]->GetMean();
	//dcorrL=(rLP-mLE)/mLE; if (fabs(dcorrL)>0.1) dcorrL=0.1*dcorrL/fabs(dcorrL);
	dcorrL=(rLP-mLE)/mLE; if (fabs(dcorrL)>0.7) dcorrL=0.7*dcorrL/fabs(dcorrL);
	if (rLP>0) drLP=
	      sqrt(pow(hdatP[ii][jj]->GetMeanError()/hdatP[ii][jj]->GetMean(),2)+
		   1.f/hdatP[ii][jj]->Integral()+
		   pow(dcorrL/(1.0+sqrt((float) nIter)),2));
	else drLP=1.e-6;
	if (fabs(dcorrL)>0.001) { 
	  //corrL*=1-20*dcorrL/(40+nIter*nIter);
	  //corrL*=1-dcorrL/(2+nIter);
	  corrL*=1-dcorrL/(1.0+sqrt((float) nIter));
	  // cout<<"!!!!!!!!corrL = "<<corrL<<endl;
	  printf("%2d : %2d / %2d / 1 %7.3f %7.3f\n",nIter,ieta,iphi,dcorrL,corrL);
	  hmapPc->SetBinContent(16-ii,jj+1,corrL);
	  hmapPc->SetBinError(16-ii,jj+1,corrL*drLP);
	  hmapP->SetBinContent(16-ii,jj+1,rLP);
	}
	else {
	  printf("%2d : %2d / %2d / 1 %7.3f %8.4f %8.4f\n",
		 nIter,ieta,iphi,dcorrL,corrL,corrL*drLP);
	  hmapP->SetBinContent(16-ii,jj+1,rLP);
	  hmapPc->SetBinError(16-ii,jj+1,corrL*drLP);
	  break;
	}
	if (nIter==nIterN-1) 
	  printf("%2d : %2d / %2d / 1 %8.4f %8.4f %8.4f\n",nIter,ieta,iphi,dcorrL,corrL,corrL*drLP);
      }
    }
  }

  printf("\nieta      eta  width  dE/dPhidEta\n");
  double xeta[16], weta[16];
  for (int i=0;i<16;i++) {
    xeta[i]=-(theHBHEEtaBounds[i]+theHBHEEtaBounds[i+1])/2;
    weta[i]=(theHBHEEtaBounds[i+1]-theHBHEEtaBounds[i]);
    mLEphi[i]*=72./weta[i];
    printf("%3d    1 %7.3f%7.3f   %g\n",-(i+1),xeta[i],weta[i],mLEphi[i]);
  }
  TCanvas *cgL = new TCanvas("cgL","cgL",300,300,600,600);
  TGraphErrors *grL = new TGraphErrors(16,xeta,mLEphi,0,0);
  grL->SetTitle("HBM ;#eta;E / #Delta#eta ,  GeV");
  grL->Draw("1+PAl");
  cgL->Print("pictHBplot_DY2/phiProfHBM1.gif");

  FILE *ft1;
  sprintf(ctit,"Run2017B_C/corrHBM1.txt",ftit,((int) Ethr1),((int) Ethr2));
  if ((ft1 = fopen(ctit,"w"))==NULL){               // Open new file
    // printf("\nNo file %s open => EXIT\n\n",file);
  return;
 }
  printf("\n\n File '%s' open \n\n",ctit);

  TH1D *hprL[16],*hprL0[16],*hprcL[16];
  TCanvas *cpr[16],*ccc[16];
  TLine *lin1 = new TLine(0,1,71,1); lin1->SetLineWidth(1);

  int noff=0;
  for (int ii=0;ii<16;ii++) {
    ieta=-ii-1;

    sprintf(ctit,"HBMcorr_%d_1",ieta);  // draw corrections
    hprcL[ii] = hmapPc->ProjectionY(ctit,16-ii,16-ii);
    hprcL[ii]->SetTitle(ctit);
    ccc[ii] = new TCanvas(ctit,ctit,800,100,500,500);
    hprcL[ii]->SetMinimum(0.41);
    hprcL[ii]->SetMaximum(hprcL[ii]->GetMaximum()*1.1);
    hprcL[ii]->SetTitleOffset(0.9,"X");
    hprcL[ii]->Draw("e");
    lin1->Draw();
    sprintf(ctit,"pictHBplot_DY2/HBMcorr_%d.gif",ieta);
    ccc[ii]->Print(ctit);

    sprintf(ctit,"HBM_E_%d_1;i#phi;GeV",ieta);
    hprL0[ii] = hmapP0->ProjectionY(ctit,16-ii,16-ii);
    hprL0[ii]->SetTitle(ctit);
    sprintf(ctit,"HBM__%d",ieta);
    hprL[ii] = hmapP->ProjectionY(ctit,16-ii,16-ii);

    cpr[ii] = new TCanvas(ctit,ctit,800,100,500,500);
    hprL0[ii]->SetFillColor(3);hprL0[ii]->SetLineColor(3);hprL0[ii]->SetLineWidth(1);
    hprL0[ii]->SetTitleOffset(0.9,"X");
    hprL0[ii]->SetMinimum(0);
    hprL0[ii]->Draw("hist");
    hprL[ii]->SetLineWidth(2);
    hprL[ii]->Draw("samehist");
    sprintf(ctit,"pictHBplot_DY2/HBM_E_%d_1.gif",ieta);
    cpr[ii]->Print(ctit);

    for (int jj=0;jj<72;jj++) {
      int ieta=-ii-1;
      int iphi=jj+1;
      corrL=hmapPc->GetBinContent(16-ii,jj+1);
      dcorrL=hmapPc->GetBinError(16-ii,jj+1);

      hcorr1D_1->Fill(corrL);
      hcorr1Derr_1->Fill(dcorrL);
      hprof->Fill(jj+1,corrL);

      noff++;

      if ((corrL<0.85)||(corrL>1.15))   fprintf(ft1,"%2d   %2d   1 %8.4f %8.4f\n",ieta,iphi,corrL,dcorrL);
    }
  }
  fclose(ft1);

  TCanvas *c1prof =new TCanvas("c1prof","c1prof",30,30,600,600);
  hprof->Draw();

  TCanvas *c1corr_1 =new TCanvas("c1corr_1","c1corr_1",30,30,600,600);
  hcorr1D_1->Draw("hist"); histStat(hcorr1D_1,1);

  TCanvas *c1corr_err_1 =new TCanvas("c1corr_err_1","c1corr_err_1",30,30,600,600);

  hcorr1Derr_1->Draw("hist");
  histStat(hcorr1Derr_1,1);

  TCanvas *c1corr_5 =new TCanvas("c1corr_5","c1corr_5",30,30,600,600);
  hcorr1D_5->Draw("hist"); histStat(hcorr1D_5,1);

  TCanvas *c1corr_err_5 =new TCanvas("c1corr_err_5","c1corr_err_5",30,30,600,600);
  hcorr1Derr_5->Draw("hist"); histStat(hcorr1Derr_5,1);

  TCanvas *c1corr_11 =new TCanvas("c1corr_11","c1corr_11",30,30,600,600);

  hcorr1D_11->Draw("hist"); histStat(hcorr1D_11,1);

  TCanvas *c1corr_err_11 =new TCanvas("c1corr_err_11","c1corr_err_11",30,30,600,600);
  hcorr1Derr_11->Draw("hist"); histStat(hcorr1Derr_11,1);


  TCanvas *ctr = new TCanvas("ctr","ctr",0,0,650,600);
  htr->Draw("hist");
  ctr->Update();

  TCanvas *chmapP = new TCanvas("chmapP","chmapP",0,0,650,600);
  chmapP->cd(); chmapP->SetRightMargin(0.12); chmapP->cd();
  hmapP->SetAxisRange(hmapP->GetBinContent(1,1)/2,-1111,"Z");
  hmapP->Draw("colz");
  chmapP->Print("pictHBplot_DY2/hmapHBM1.gif");
  chmapP->Update();

  TCanvas *chmapP0 = new TCanvas("chmapP0","chmapP0",0,0,650,600);
  chmapP0->cd(); chmapP0->SetRightMargin(0.12); chmapP0->cd(); 
  hmapP0->SetAxisRange(hmapP0->GetBinContent(1,1)/2,-1111,"Z");
  hmapP0->Draw("colz");
  chmapP0->Print("pictHBplot_DY2/hmap0HBM1.gif");
  chmapP0->Update();

  TCanvas *chmapPc = new TCanvas("chmapPc","chmapPc",0,0,650,600);
  chmapPc->cd(); chmapPc->SetRightMargin(0.12);
  hmapPc->SetAxisRange(0.6,2,"Z");
  hmapPc->Draw("colz");
  chmapPc->Print("pictHBplot_DY2/hmapcHBM1.gif");
  chmapPc->Update();

  sprintf(ctit,"HBM1o_%s_%d_%d.root",ftit,((int) Ethr1),((int) Ethr2));
  TFile *histf = new TFile(ctit,"RECREATE");
  hmapP->Write(); 
  hmapP0->Write(); 
  hmapPc->Write(); 
  histf->Close();

}

