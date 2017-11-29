Iterative method for HCAL Phi Symmetry.

1) Prepare the enviroment

-setup CMSSW

-git clone https://github.com/CMSHCALCalib/PhiSymmetry/Iterative


2) Prepare RecHit energy spectra for HB, HE, HF: src/phiSym.cc, test/d92XphiSymRECO.py

3) produce PhiSym calibration coefficients 
  
   for HB depth 1,2
   
   -------  cetaflatHBMt.C, cetaflatHBPt.C, cetaflatHBM2t.C, cetaflatHBP2t.C
     
   for HE depth 1,2,3 
    
   -------  cetaflatHEP1t.C, cetaflatHEP2t.C, cetaflatHEP3t.C,
   -------  cetaflatHEM1t.C, cetaflatHEM2t.C, cetaflatHEM3t.C)
     
   for HF long and short fibers
   
   -------   cetaflatHFM12.C, cetaflatHFP12.C
