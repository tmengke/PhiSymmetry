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
//#include <iostream.h>
//#include "histspec.C"
#include "histStat.C"

void cetaflatHEP2t(int nIterN=1, double Ethr1=4, double Ethr2=150) {
  
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
  gStyle->SetTitleX(0.3);
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

    sprintf(ctit, "Run2017B.root");

  TFile *fila = new TFile (ctit);
  cout<<"File= "<<ctit<<endl;

  TH1F *hcounter =   new TH1F(*((TH1F*)fila->Get("phaseHF/hcounter")));
  cout<<"Stat= "<<hcounter->GetBinContent(2)<<endl;
  cout<<"E within: "<<Ethr1<<" - "<<Ethr2<<endl;

  TH2F* hmapP = new TH2F("hmapP","E  HEP;i#eta;i#phi",14,15.5,29.5,72,0,72);
  TH2F* hmapP0 = new TH2F("hmapP0","E0  HEP;i#eta;i#phi",14,15.5,29.5,72,0,72);
  TH2F* hmapPc = new TH2F("hmapPc","corr  HEP;i#eta;i#phi",14,15.5,29.5,72,0,72);
  hmapPc->Sumw2();
  TH1F *hcorr1D = new TH1F("hcorr1D","Corr",150,0.5,2);
  TH1F *ht = new TH1F("ht","ht",20000,0,1e7);
  TH1F *htx = new TH1F("htx","htx",20000,0,1e5);
  TH1F *htr = new TH1F("htr","htr",5000,0,3);

  

  TH1F *hdatP[14][72], *hdatPx[14][72];
  for (int ii=0;ii<14;ii++) for (int jj=0;jj<72;jj++) {
    sprintf(ctit,"h%d_%d",ii+16,jj+1);
    hdatP[ii][jj] = new TH1F(ctit,ctit,10000,0,250);
  }

  TCanvas *cx[400];
  TSpline5 *tt;

  Double_t x,y,rPL,rPS,mLE,mSE,ermean,rms;
  Double_t xx[4000],yy[4000];
  Int_t nELP, nESP, nIter=0;
  Double_t mcorrL,scorrL,mcorrS,scorrS,erLP,erSP,rLP,drLP,rSP,corrL,corrS,dcorrL,dcorrS;
  double mLEphi[14];

  TH1F *spectra26 = new TH1F("spectra26","HEP17",1000,0.,250.);
  TH1F *spectra26n = new TH1F("spectra26n","!HEP17",1000,0.,250.);
  TH1F *spectra27 = new TH1F("spectra27","HEP17",1000,0.,250.);
  TH1F *spectra27n = new TH1F("spectra27n","!HEP17",1000,0.,250.);
  TH1F *spectra28 = new TH1F("spectra28","HEP17",1000,0.,250.);
  TH1F *spectra28n = new TH1F("spectra28n","!HEP17",1000,0.,250.);



  TCanvas *ccxx = new TCanvas("ccxx","ccxx",0,400,700,400);
  ccxx->Divide(2,1);
  for (int ii=0;ii<14;ii++) {
    int ieta=ii+16;

    mLE=mSE=0;   // ------------------for initial condition
    int nmLE=0, nmSE=0;
    ht->Reset(); htx->Reset();
    for (int ll=0;ll<72;ll++) {
      int iphi=ll+1;
      if (abs(ieta)==29) continue;
      if (abs(ieta)==16) continue;

      if (abs(ieta)<18) continue;
      if (abs(ieta)>20 && iphi%2==0) continue;  

      hmapPc->SetBinContent(ii+1,ll+1,1);
      hmapPc->SetBinError(ii+1,ll+1,1.e-6);
      sprintf(ctit,"phaseHF/eHEspec/E_+%d_%d_2",ieta,iphi);
      hdatPx[ii][ll]  =   new TH1F(*((TH1F*)fila->Get(ctit)));



      if ((ieta==26)&&(iphi >= 63) && (iphi <= 66)) spectra26->Add( hdatPx[ii][ll]);
      if ((ieta==27)&&(iphi >= 63) && (iphi <= 66)) spectra27->Add( hdatPx[ii][ll]);
      if ((ieta==28)&&(iphi >= 63) && (iphi <= 66)) spectra28->Add( hdatPx[ii][ll]);
      if ((ieta==26)&&((iphi < 63) || (iphi > 66))) spectra26n->Add( hdatPx[ii][ll]);
      if ((ieta==27)&&((iphi < 63) || (iphi > 66))) spectra27n->Add( hdatPx[ii][ll]);
      if ((ieta==28)&&((iphi < 63) || (iphi > 66))) spectra28n->Add( hdatPx[ii][ll]);



      hdatPx[ii][ll]->SetAxisRange(Ethr1,Ethr2);



      rLP = hdatPx[ii][ll]->Integral()*hdatPx[ii][ll]->GetMean();
      hmapP0->SetBinContent(ii+1,ll+1,rLP);
      if (rLP>0) {
	ht->Fill(rLP); htx->Fill(rLP);
	mLE += rLP;
	nmLE++;
	drLP=rLP*sqrt(pow(1./hdatPx[ii][ll]->Integral(),2)+
		      pow(hdatPx[ii][ll]->GetMeanError()/hdatPx[ii][ll]->GetMean(),2));
	hmapP0->SetBinError(ii+1,ll+1,drLP);
      }
      else hmapP0->SetBinError(ii+1,ll+1,0);
    }
    if (nmLE>0) mLE /= nmLE; 
    else mLE=0;
    ccxx->cd(1); ht->Draw("hist");
    ccxx->cd(2); htx->Draw("hist");
    ccxx->Update();
    if (htx->GetBinContent(20001)>1) histspec(ht,mLE,ermean,rms,4,-5);
    else histspec(htx,mLE,ermean,rms,4,-5);
    mLEphi[ii]=mLE;
    printf("ieta %2d :  <E>= %8.1f \n",ieta,mLE);
    if (ht->GetMean()>0) htr->Fill(ht->GetRMS()/ht->GetMean());
    
    for (int jj=0;jj<72;jj++) {
      int iphi=jj+1;

      if (abs(ieta)<18) continue;
      if (abs(ieta)>20 && iphi%2==0) continue;  

      for (nIter=1;nIter<nIterN;nIter++) { //cout<<nIter<<" |  ";
	if (hmapP0->GetBinContent(ii+1,jj+1)<=0) continue;
	corrL=hmapPc->GetBinContent(ii+1,jj+1);
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
	dcorrL=(rLP-mLE)/mLE; if (fabs(dcorrL)>0.5) dcorrL=0.5*dcorrL/fabs(dcorrL);
	if (rLP>0) drLP=
	      sqrt(pow(hdatP[ii][jj]->GetMeanError()/hdatP[ii][jj]->GetMean(),2)+
		   1.f/hdatP[ii][jj]->Integral()+
		   pow(dcorrL/(1.0+sqrt((float) nIter)),2));
	else drLP=1.e-6;
	if (fabs(dcorrL)>0.001) { 
	  //corrL*=1-20*dcorrL/(40+nIter*nIter);
	  //corrL*=1-dcorrL/(2+nIter);
	  corrL*=1-dcorrL/(1.0+sqrt((float) nIter));
	  //printf("%2d : %2d / %2d / 1 %7.3f %7.3f\n",nIter,ieta,iphi,dcorrL,corrL);
	  hmapPc->SetBinContent(ii+1,jj+1,corrL);
	  hmapPc->SetBinError(ii+1,jj+1,corrL*drLP);
	  hmapP->SetBinContent(ii+1,jj+1,rLP);
	}
	else {
	  printf("%2d : %2d / %2d / 2 %7.3f %8.4f %8.4f\n",
		 nIter,ieta,iphi,dcorrL,corrL,corrL*drLP);
	  hmapP->SetBinContent(ii+1,jj+1,rLP);
	  hmapPc->SetBinError(ii+1,jj+1,corrL*drLP);
	  break;
	}
	if (nIter==nIterN-1) 
	  printf("%2d : %2d / %2d / 2 %8.4f %8.4f %8.4f\n",nIter,ieta,iphi,dcorrL,corrL,corrL*drLP);
      }
    }
  }

  printf("\nieta      eta  width  dE/dPhidEta\n");
  double xeta[14], weta[14], yield[14];
  int ind=0;
  for (int i=0;i<14;i++) {
    int ieta=i+16;
    if (abs(ieta)<18) continue;
    xeta[ind]=(theHBHEEtaBounds[i+15]+theHBHEEtaBounds[i+16])/2;
    weta[ind]=(theHBHEEtaBounds[i+16]-theHBHEEtaBounds[i+15]);
    yield[ind]=mLEphi[i];
    if (abs(ieta)<21) yield[ind]*=72/weta[ind];
    else yield[ind]*=36/weta[ind];
    printf("%3d   2 %7.3f%7.3f   %g\n",ieta,xeta[ind],weta[ind],yield[ind]);
    ind++;
  }
  TCanvas *cgL = new TCanvas("cgL","cgL",300,300,600,600);
  TGraphErrors *grL = new TGraphErrors(ind,xeta,yield,0,0);
  grL->SetTitle("HEP ;#eta;E / #Delta#eta ,  GeV");
  grL->Draw("1+PAl");
  cgL->Print("pictHEplot/phiProfHEP2.gif");

  FILE *ft1;
  sprintf(ctit,"TEST/HEP2.txt",ftit,((int) Ethr1),((int) Ethr2));
  if ((ft1 = fopen(ctit,"w"))==NULL){               // Open new file
    //printf("\nNo file %s open => EXIT\n\n",file);
    return;
  }
  printf("\n\n File '%s' open \n\n",ctit);

  TH1D *hprL[14],*hprL0[14],*hprcL[16];
  TCanvas *cpr[14],*ccc[16];
  TLine *lin1 = new TLine(0,1,71,1); lin1->SetLineWidth(1);

  int noff=0;
  for (int ii=0;ii<14;ii++) {

    int ieta=ii+16;
    if (abs(ieta)<18) continue;

    sprintf(ctit,"HEPcorr_%d_2",ieta);  // draw corrections
    hprcL[ii] = hmapPc->ProjectionY(ctit,ii+1,ii+1);
    hprcL[ii]->SetTitle(ctit);
    ccc[ii] = new TCanvas(ctit,ctit,800,100,500,500);
    hprcL[ii]->SetMinimum(0.41);
    hprcL[ii]->SetMaximum(hprcL[ii]->GetMaximum()*1.1);
    hprcL[ii]->SetTitleOffset(0.9,"X");
    hprcL[ii]->Draw("e");
    lin1->Draw();
    sprintf(ctit,"pictHEplot/HEP2corr_%d.gif",ieta);
    ccc[ii]->Print(ctit);

    sprintf(ctit,"HEP_E_%d_2;i#phi;GeV",ieta);
    hprL0[ii] = hmapP0->ProjectionY(ctit,ii+1,ii+1);
    hprL0[ii]->SetTitle(ctit);
    sprintf(ctit,"HEP__%d",ieta);
    hprL[ii] = hmapP->ProjectionY(ctit,ii+1,ii+1);
    if (ieta>20) {
      hprL[ii]->Rebin();
      hprL0[ii]->Rebin();
    }
    cpr[ii] = new TCanvas(ctit,ctit,800,100,500,500);
    hprL0[ii]->SetFillColor(3);hprL0[ii]->SetLineColor(3);hprL0[ii]->SetLineWidth(1);
    hprL0[ii]->SetTitleOffset(0.9,"X");
    hprL0[ii]->SetMinimum(0);
    hprL0[ii]->Draw("hist");
    hprL[ii]->SetLineWidth(2);
    hprL[ii]->Draw("samehist");
    sprintf(ctit,"pictHEplot/HEP_E_%d_2.gif",ieta);
    cpr[ii]->Print(ctit);

    for (int jj=0;jj<72;jj++) {
      int ieta=ii+16;
      int iphi=jj+1;

      if (abs(ieta)>20 && iphi%2==0) continue; 
 
      corrL=hmapPc->GetBinContent(ii+1,jj+1);
      dcorrL=hmapPc->GetBinError(ii+1,jj+1);

      if ((ieta>25)||(ieta<18)){

              if ((ieta!=16)&&(ieta!=29)){                                                                                                                                                     
			if ((iphi == 63)||(iphi==64)||(iphi==65)||(iphi==66)) hcorr1D->Fill(corrL);
	}                                                                                                                                                                            
      }


      noff++;
      //printf("%2d : %2d / %2d / 2 %8.4f %8.4f\n",noff,ieta,iphi,corrL,dcorrL);
      //      if ((corrL<0.85)||(corrL>1.15))
 fprintf(ft1,"%2d   %2d   2 %8.4f %8.4f\n",ieta,iphi,corrL,dcorrL);
    }
  }
  fclose(ft1);

  TCanvas *c1corr =new TCanvas("c1corr","c1corr",30,30,600,600);
  hcorr1D->Draw("hist");  histStat(hcorr1D,1);
  c1corr->Print("pictHEplot/corrHEP2.gif");
  //hcorr1D->Write(); 

  TCanvas *ctr = new TCanvas("ctr","ctr",0,0,650,600);
  htr->Draw("hist");
  ctr->Update();

  TCanvas *chmapP = new TCanvas("chmapP","chmapP",0,0,650,600);
  chmapP->cd(); chmapP->SetRightMargin(0.12); chmapP->SetLogz();
  hmapP->SetAxisRange(hmapP->GetBinContent(3,1)/2,-1111,"Z");
  hmapP->Draw("colz");
  chmapP->Print("pictHEplot/hmapHEP2.gif");
  chmapP->Update();

  TCanvas *chmapP0 = new TCanvas("chmapP0","chmapP0",0,0,650,600);
  chmapP0->cd(); chmapP0->SetRightMargin(0.12); chmapP0->SetLogz();
  hmapP0->SetAxisRange(hmapP0->GetBinContent(3,1)/2,-1111,"Z");
  hmapP0->Draw("colz");
  chmapP0->Print("pictHEplot/hmap0HEP2.gif");
  chmapP0->Update();

  TCanvas *chmapPc = new TCanvas("chmapPc","chmapPc",0,0,650,600);
  chmapPc->cd(); chmapPc->SetRightMargin(0.12);
  hmapPc->SetAxisRange(0.6,2,"Z");
  hmapPc->Draw("colz");
  chmapPc->Print("pictHEplot/hmapcHEP2.gif");
  chmapPc->Update();

  sprintf(ctit,"HEP2o_%s_%d_%d.root",ftit,((int) Ethr1),((int) Ethr2));
  TFile *histf = new TFile(ctit,"RECREATE");
  hmapP->Write(); 
  hmapP0->Write(); 
  hmapPc->Write(); 
  histf->Close();


  TCanvas *chmapppp26 = new TCanvas("chmapppp26","chmapppp26",0,0,650,600);
  chmapppp26->cd(); chmapppp26->SetRightMargin(0.12);
  TLegend *legend26 = new TLegend(0.4,0.6,0.89,0.89);

  Double_t norm = 1;
  Double_t scale = norm/(spectra26->Integral());
  // spectra26->Scale(scale);                                                                                                                                                                             
  spectra26->Rebin(20);
  spectra26->Draw();

  norm = 1;
  scale = norm/(17*spectra26n->Integral());
  spectra26n->Scale(1/17.);
  spectra26n->Rebin(20);

  spectra26n->Draw("same");

  legend26->AddEntry(spectra26,"HEP17", "l");
  legend26->AddEntry(spectra26n,"!HEP17", "l");
  legend26->Draw();


  TCanvas *chmapppp27 = new TCanvas("chmapppp27","chmapppp27",0,0,650,600);
  chmapppp27->cd(); chmapppp27->SetRightMargin(0.12);

  TLegend *legend27 = new TLegend(0.4,0.6,0.89,0.89);

  Double_t norm = 1;
  Double_t scale = norm/(spectra27->Integral());
  // spectra27->Scale(scale);                                                                                                                                                                             
  spectra27->Rebin(20);
  spectra27->Draw();

  norm = 1;
  scale = norm/(17*spectra27n->Integral());
  spectra27n->Scale(1/17.);
  spectra27n->Rebin(20);
  spectra27n->Draw("same");


  legend27->AddEntry(spectra27,"HEP17", "l");
  legend27->AddEntry(spectra27n,"!HEP17", "l");
  legend27->Draw();


  TCanvas *chmapppp28 = new TCanvas("chmapppp28","chmapppp28",0,0,650,600);
  chmapppp28->cd(); chmapppp28->SetRightMargin(0.12);
  TLegend *legend28 = new TLegend(0.4,0.6,0.89,0.89);

  Double_t norm = 1;
  Double_t scale = norm/(spectra28->Integral());
  // spectra28->Scale(scale);                                                                                                                                                                             
  spectra28->Rebin(20);
  spectra28->Draw();

  norm = 1;
  scale = norm/(17*spectra28n->Integral());
  spectra28n->Scale(1/17.);
  spectra28n->Rebin(20);
  spectra28n->Draw("same");
  legend28->AddEntry(spectra28,"HEP17", "l");
  legend28->AddEntry(spectra28n,"!HEP17", "l");

  legend28->Draw();

}

