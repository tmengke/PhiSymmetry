#include "TH1F.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TSpline.h"
#include "TLine.h"
#include "TKey.h"
#include "TClass.h"
#include "TStyle.h"
#include <map>
#include <stdio.h>
#include <sstream>
#include <iostream>
#include <string>

//g++ `root-config --cflags` -o Iterative Iterative.cc  `root-config --glibs`

using namespace std;


int main(){

  gStyle->SetOptLogz(0);
  gStyle->SetMarkerSize(0.7);
  gStyle->SetMarkerStyle(20);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetTitleOffset(1.7,"Y");
  gStyle->SetTitleOffset(0.9,"X");
  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetStatH(0.09);
  gStyle->SetStatW(0.3);
  gStyle->SetTitleW(0.4);
  gStyle->SetTitleX(0.3);
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();


int nIterN=5;
double Ethr1=4;
double Ethr2=150;

char ctit[145];

std::map<std::string, TH1F*> histo1F;
std::map<std::string, TH2F*> histo2F;
  
std::string det[] = {"HEP","HEM","HBP","HBM","HFP","HFM"};
const double xl[] = {15.5,-29.5,0.5,-16.5,28.5,-41.5};
const double xu[] = {29.5,-15.5,16.5,-0.5,41.5,-28.5};
const int ietaN[] = {14,14,16,16,13,13};
const int sieta[] = {16,16,1,1,29,29};
const int iphiN[] = {72,72,72,72,36,36};

sprintf(ctit, "phisym.root");
TFile *fila = new TFile (ctit);
cout<<"File= "<<ctit<<endl;
TH1F *hcounter =   new TH1F(*((TH1F*)fila->Get("phaseHF/hcounter")));
cout<<"Stat= "<<hcounter->GetBinContent(2)<<endl;
cout<<"E within: "<<Ethr1<<" - "<<Ethr2<<endl;

TFile *histf = new TFile("itercali.root","recreate");

//booking histgrams
std::map<std::string, TDirectory*> cdDet;

for(int i=0;i<6;i++){
	for(int j=1;j<4;j++){
		if(i>1&&j==3)continue;
	    	std::string hname=det[i]+"_depth"+std::to_string(j);
	  	cdDet[hname] = histf->mkdir(hname.c_str());
        	cdDet[hname]->cd();
		hname="E_"+det[i]+std::to_string(j);  
		histo2F[hname] = new TH2F(hname.c_str(),hname.c_str(),ietaN[i],xl[i],xu[i],iphiN[i],0,72);
		hname="E0_"+det[i]+std::to_string(j);
		histo2F[hname] = new TH2F(hname.c_str(),hname.c_str(),ietaN[i],xl[i],xu[i],iphiN[i],0,72);
		hname="corr_"+det[i]+std::to_string(j);
		histo2F[hname] = new TH2F(hname.c_str(),hname.c_str(),ietaN[i],xl[i],xu[i],iphiN[i],0,72);
		histo2F[hname]->Sumw2();
		hname="hcorr_"+det[i]+std::to_string(j);
		histo1F[hname]=new TH1F(hname.c_str(),hname.c_str(),300,0.,2.);
	}
}

//loop over all subdetectors and depths

for(int i=0;i<1;i++){

	for(int j=1;j<4;j++){

		if(i>1&&j==3)continue;
		std::map<std::string, TH1F*> hdatP;
		std::map<std::string, TH1F*> hdatPx;		

		std::string dir;
	
		std::string hname = "E0_"+det[i]+std::to_string(j);
		std::string cname = "corr_"+det[i]+std::to_string(j);
                std::string ename = "E_"+det[i]+std::to_string(j);

		
		switch (i) {
			case 0:
			case 1: dir="phaseHF/HEcollapsed";
				break;
			case 2:
			case 3: dir="phaseHF/eHBspec";
				break;
			case 4:
			case 5: dir="phaseHF/espec";
				break;
		}	

		TDirectory *fs=(TDirectory*)fila->Get(dir.c_str());
		TIter next(fs->GetListOfKeys());
		TKey *key;
		while ((key = (TKey*)next())) {
			TClass *cl = gROOT->GetClass(key->GetClassName());
			/*if (!cl->InheritsFrom("TH1F")) {
				std::cout<<"Not a histogram"<<std::endl;;
				continue;
			}*/
			string pname = key->GetName();
			/*string pname_tmp = pname;
			for ( std::string::iterator it=pname_tmp.begin(); it!=pname_tmp.end(); ++it){
  			 //cout << *it<<endl;
   				 if(!isdigit(*it)){
       				 *it =' ';
       				 }
   			}	

  			Int_t ieta,iphi,depth;
			stringstream ss(pname_tmp);
			ss>>ieta>>iphi>>depth;	

			cout<<"ieta  "<<ieta<<"iphi   "<<iphi<<"depth   "<<depth<<endl;
			
			if(i==0||i==2||i==4){if(pname[2]!='+')continue;}
			else {if(pname[2]!='-')continue;}			
			if(j!=depth)continue;  //depths
			if(ieta==29)continue;// no calibration for ieta 29
			*/


			if(pname=="E_+17_58_1")continue;
			hdatP[pname] = new TH1F(pname.c_str(),pname.c_str(),10000,0,250);
			hdatPx[pname] =	(TH1F*)key->ReadObj();
			hdatPx[pname] -> SetAxisRange(Ethr1,Ethr2);

				
		}

		for(int ii=0;ii<ietaN[i];ii++){

		   	Int_t ieta=ii+sieta[i];
			if(ieta==29)continue;
			if(i==1||i==3||i==5)ieta*=-1;	
			Double_t x,mLE=0;
			Double_t xx[4000],yy[4000];		
			Double_t rLP,drLP,corrL,dcorrL;
			Int_t nIter=0;
	
			Int_t nmLE=0;
			TSpline5 *tt;

			for(int jj=0;jj<72;jj++){
				 Int_t iphi=jj+1;
				 Int_t ll=jj;

				 if(i>3){
					if(iphi%2==0)continue;
					ll = iphi/2;
				 }

				 std::string pname;
				 if(ieta>0)pname="E_+"+std::to_string(abs(ieta))+"_"+std::to_string(iphi)+"_"+std::to_string(j);
				 else pname="E_-"+std::to_string(abs(ieta))+"_"+std::to_string(iphi)+"_"+std::to_string(j);
				 
				if(hdatPx.count(pname.c_str())==0)continue;
			
				 if(ieta>0){
				 	histo2F[cname]->SetBinContent(ii+1,ll+1,1);
                        	 	histo2F[cname]->SetBinError(ii+1,ll+1,1.e-6);
                        	 }
                        	 else{ 
                        	 	histo2F[cname]->SetBinContent(ietaN[i]-ii,ll+1,1);
                        		histo2F[cname]->SetBinError(ietaN[i]-ii,ll+1,1.e-6);
				 }
				 rLP = hdatPx[pname]->Integral()*hdatPx[pname]->GetMean();

				if(ieta>0) histo2F[hname]->SetBinContent(ii+1,ll+1,rLP);
                                else histo2F[hname]->SetBinContent(ietaN[i]-ii,ll+1,rLP);

				if (rLP>0) {
        				mLE += rLP;
        				nmLE++;
        				drLP=rLP*sqrt(pow(1./hdatPx[pname]->Integral(),2)+
                      			pow(hdatPx[pname]->GetMeanError()/hdatPx[pname]->GetMean(),2));
					if(ieta>0)histo2F[hname]->SetBinError(ii+1,ll+1,drLP);
					else histo2F[hname]->SetBinError(ietaN[i]-ii,ll+1,drLP);
		        	
				}
      				else {
					if(ieta>0)histo2F[hname]->SetBinError(ii+1,ll+1,0);
					histo2F[hname]->SetBinError(ietaN[i]-ii,ll+1,0);
				}

			}

			if (nmLE>0) mLE /= nmLE;
    			else mLE=0;
			printf("ieta %2d :  <E>= %8.1f \n",ieta,mLE);
		
			for (int ll=0;ll<72;ll++) {
				 Int_t iphi=ll+1;
				 Int_t jj=ll;
				 if(i>3){
					if(iphi%2==0)continue;
					jj = iphi/2;
					}
                                 std::string pname;
                                 if(ieta>0)pname="E_+"+std::to_string(abs(ieta))+"_"+std::to_string(iphi)+"_"+std::to_string(j);
                                 else pname="E_-"+std::to_string(abs(ieta))+"_"+std::to_string(iphi)+"_"+std::to_string(j);
                                 if(hdatPx.count(pname.c_str())==0)continue;

				 for (nIter=1;nIter<nIterN;nIter++) {

					 std::string cname = "corr_"+det[i]+std::to_string(j);
					 std::string ename = "E_"+det[i]+std::to_string(j);  
				 	 if(ieta>0){
						 if(histo2F[hname]->GetBinContent(ii+1,jj+1)<=0) continue;
					 	 corrL=histo2F[cname]->GetBinContent(ii+1,jj+1);
					 }
					 else {
						 if(histo2F[hname]->GetBinContent(ietaN[i]-ii,jj+1)<=0) continue;
					         corrL=histo2F[cname]->GetBinContent(ietaN[i]-ii,jj+1);
					 }
				  	 hdatP[pname]->Reset();
					 for (int kk=1;kk<=hdatPx[pname]->GetNbinsX();kk++) {
         					 xx[kk-1]=hdatPx[pname]->GetBinCenter(kk);
          				         yy[kk-1]=hdatPx[pname]->GetBinContent(kk);
        				 }
        				 tt = new TSpline5("tt",xx,yy,1000,"",10,20);

        				 for (int kk=1;kk<=hdatP[pname]->GetNbinsX();kk++) {
          					 x=hdatP[pname]->GetBinCenter(kk);
          					 hdatP[pname]->Fill(x*corrL,tt->Eval(x)/10.0);
        				 }
        				 tt->Delete();

					 hdatP[pname]->SetAxisRange(Ethr1,Ethr2);
        				 rLP = hdatP[pname]->Integral()*hdatP[pname]->GetMean();
        				 dcorrL=(rLP-mLE)/mLE; if (fabs(dcorrL)>0.7) dcorrL=0.7*dcorrL/fabs(dcorrL);
					 if (rLP>0) drLP=sqrt(pow(hdatP[pname]->GetMeanError()/hdatP[pname]->GetMean(),2)+1.f/hdatP[pname]->Integral()+pow(dcorrL/(1.0+sqrt((float) nIter)),2));
        				 else drLP=1.e-6;
					 if (fabs(dcorrL)>0.001) {
					 corrL*=1-dcorrL/(1.0+sqrt((float) nIter));
					 	if(ieta>0){
							histo2F[cname]->SetBinContent(ii+1,jj+1,corrL);
						 	histo2F[cname]->SetBinError(ii+1,jj+1,corrL*drLP);
          				 	 	histo2F[ename]->SetBinContent(ii+1,jj+1,rLP);
					 	}
						else {
							histo2F[cname]->SetBinContent(ietaN[i]-ii,jj+1,corrL);
                                                        histo2F[cname]->SetBinError(ietaN[i]-ii,jj+1,corrL*drLP);
                                                        histo2F[ename]->SetBinContent(ietaN[i]-ii,jj+1,rLP);

						}

					 }
					 else {
						printf("%2d : %2d / %2d / %d %7.3f %8.4f %8.4f\n",
                 					nIter,ieta,iphi,j,dcorrL,corrL,corrL*drLP);
						if(ieta>0){
                                                        histo2F[cname]->SetBinError(ii+1,jj+1,corrL*drLP);
                                                        histo2F[ename]->SetBinContent(ii+1,jj+1,rLP);
                                                }
                                                else {
                                                        histo2F[cname]->SetBinError(ietaN[i]-ii,jj+1,corrL*drLP);
                                                        histo2F[ename]->SetBinContent(ietaN[i]-ii,jj+1,rLP);
						}
                                                
						break;
						

					 }
					 if (nIter==nIterN-1)printf("%2d : %2d / %2d / %d %8.4f %8.4f %8.4f\n",nIter,ieta,iphi,j,dcorrL,corrL,corrL*drLP);

				 } 



			}	

		}		

	FILE *ft1;
	std::string txtf="txtfile/"+det[i]+std::to_string(j)+".txt"; 
  	if ((ft1 = fopen(txtf.c_str(),"w"))==NULL){               // Open new file
        	return 0;
        }
        printf("\n\n File '%s' open \n\n",txtf.c_str());


	std::string dname=det[i]+"_depth"+std::to_string(j);
        cdDet[dname]->cd();

	TH1D *hprL[ietaN[i]],*hprL0[ietaN[i]],*hprcL[ietaN[i]];
  	TCanvas *cpr[ietaN[i]],*ccc[ietaN[i]];
  	TLine *lin1 = new TLine(0,1,71,1); lin1->SetLineWidth(1);

	int noff=0;
	for (int ii=0;ii<ietaN[i];ii++) {

		Int_t ieta=ii+sieta[i];
        	if(i==1||i==3||i==5)ieta*=-1;

		sprintf(ctit,"%s_corr_%d_%d",det[i].c_str(),ieta,j); 
		if(ieta>0)hprcL[ii] = histo2F[cname]->ProjectionY(ctit,ii+1,ii+1);
 	   	else hprcL[ii] = histo2F[cname]->ProjectionY(ctit,ietaN[i]-ii,ietaN[i]-ii);
		if(hprcL[ii]->GetEntries()==0)continue;
		hprcL[ii]->SetTitle(ctit);
    		ccc[ii] = new TCanvas(ctit,ctit,800,100,500,500);
    		hprcL[ii]->SetMinimum(0.41);
    		hprcL[ii]->SetMaximum(hprcL[ii]->GetMaximum()*1.1);
    		hprcL[ii]->SetTitleOffset(0.9,"X");
    		hprcL[ii]->Draw("e");
    		lin1->Draw();
		ccc[ii]->Modified();
  		ccc[ii]->Update();
  		ccc[ii]->Write();	

		sprintf(ctit,"%s_E0_%d_%d;i#phi;",det[i].c_str(),ieta,j);
    		if(ieta>0)hprL0[ii] = histo2F[hname]->ProjectionY(ctit,ii+1,ii+1);
   		else hprL0[ii] = histo2F[hname]->ProjectionY(ctit,ietaN[i]-ii,ietaN[i]-ii);
		hprL0[ii]->SetTitle(ctit);
    		sprintf(ctit,"%s_E_%d_%d;i#phi;",det[i].c_str(),ieta,j);
		if(ieta>0)hprL[ii] = histo2F[ename]->ProjectionY(ctit,ii+1,ii+1);
		else hprL[ii] = histo2F[ename]->ProjectionY(ctit,ietaN[i]-ii,ietaN[i]-ii);
    		if (i<2&&abs(ieta)>20) {
      		hprL[ii]->Rebin();
      		hprL0[ii]->Rebin();
    		}
		if (i>3&&abs(ieta)>39){
		hprL[ii]->Rebin(2);
                hprL0[ii]->Rebin(2);
		}	
    		cpr[ii] = new TCanvas(ctit,ctit,800,100,500,500);
    		hprL0[ii]->SetFillColor(3);hprL0[ii]->SetLineColor(3);hprL0[ii]->SetLineWidth(1);
    		hprL0[ii]->SetTitleOffset(0.9,"X");
    		hprL0[ii]->SetMinimum(0);
    		hprL0[ii]->Draw("hist");
    		hprL[ii]->SetLineWidth(2);
    		hprL[ii]->Draw("samehist");
    		cpr[ii]->Modified();
		cpr[ii]->Update();
		cpr[ii]->Write();

		for (int ll=0;ll<72;ll++) {
			Int_t iphi=ll+1;
	        	Int_t jj=ll;
			
			if(i>3){
			if(iphi%2==0)continue;
			jj = iphi/2;
			}
			if (abs(ieta)==29) continue;
			if(i<2){if (abs(ieta)>20 && iphi%2==0) continue;}
			
			double_t corrL,dcorrL;
			if(ieta>0){
			corrL=histo2F[cname]->GetBinContent(ii+1,jj+1);
      			dcorrL=histo2F[cname]->GetBinError(ii+1,jj+1);
			}
			else {
			corrL=histo2F[cname]->GetBinContent(ietaN[i]-ii,jj+1);
                     	dcorrL=histo2F[cname]->GetBinError(ietaN[i]-ii,jj+1);	
			}
			std::string crname="hcorr_"+det[i]+std::to_string(j);
                	histo1F[crname]->Fill(corrL);
			noff++;
			printf("%2d : %2d / %2d / %d %8.4f %8.4f\n",noff,ieta,iphi,j,corrL,dcorrL);
			fprintf(ft1,"%2d   %2d   %d %8.4f %8.4f\n",ieta,iphi,j,corrL,dcorrL);
			}

	    }


	TCanvas *chisto1F = new TCanvas("chisto1F","chisto1F",0,0,650,600);
	std::string crname="hcorr_"+det[i]+std::to_string(j);
	histo1F[crname]->Draw("HIST");
        chisto1F->Update();
	
	fclose(ft1);	

	TCanvas *chisto2F = new TCanvas("chisto2F","chisto2F",0,0,650,600);
  	chisto2F->cd(); chisto2F->SetRightMargin(0.12); chisto2F->SetLogz();
  	//histo2F[ename]->SetAxisRange(histo2F[ename]->GetBinContent(2,1)/2,-1111,"Z");
  	histo2F[ename]->Draw("colz");
  	chisto2F->Update();

  	TCanvas *chisto2F0 = new TCanvas("chisto2F0","chisto2F0",0,0,650,600);
  	chisto2F0->cd(); chisto2F0->SetRightMargin(0.12); chisto2F0->SetLogz();
  	//histo2F[hname]->SetAxisRange(histo2F[hname]->GetBinContent(2,1)/2,-1111,"Z");
  	histo2F[hname]->Draw("colz");
  	chisto2F0->Update();
  	
	TCanvas *chisto2Fc = new TCanvas("chisto2Fc","chisto2Fc",0,0,650,600);
  	chisto2Fc->cd(); chisto2Fc->SetRightMargin(0.12);
  	histo2F[cname]->SetAxisRange(0.6,2,"Z");
  	histo2F[cname]->Draw("colz");
  	chisto2Fc->Update();

	chisto1F->Write();
  	chisto2F->Write();
  	chisto2F0->Write();
  	chisto2Fc->Write();

	delete chisto2F;
	delete chisto2F0;
	delete chisto2Fc;
	delete chisto1F;	

	hdatP.clear();
	hdatPx.clear();

 	}

      }


cout<<"finish  ---------"<<endl;

return 0;

}
