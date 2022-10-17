# FCLS_IDL  
 
; This program is used for linear spectral unmixing with fully constrained least squared method (FCLS) for multi-spectral image.  
;  
## Linear Spectral Mixture Model  
;      Mixed_spectrum=f1×endmember_1 + f2×endmember_2 +...+ fn×endmember_n  
;      constrains: f1+f2+...+fn=1 & f1, f2,... fn >=0  
;      n: number of endmembers    
  
## METHOD:  
;  Fully Constrained Least Squared Method: Fraction is solved with sum-to-one and nonnegative constrained method  
;  This FCLS code is based on a method proposed by Dr. Xi LI in Wuhan University. The method was documented in an email from Dr. Xili to Xuehong Chen
;  The code is designed in an vector-style with high efficenciy  

## HOW TO USE   
### Option 1: Compile and Run the program in ENVI+IDL environment, then set input and output data path through GUI.   

 INPUT:  
;  Multispectral Image: ENVI standard file  
;  Endmember: csv file (first line is endmember name, first row is band number/wavelength, others are endmember data)  
;      For Example  
;     “  
;         wavelength(nm), vegetation,  impervious,  soil  
;             442.6, 0.031106667, 0.18706, 0.101763333  
;             551.8, 0.035743333, 0.191026667, 0.13075  
;             632, 0.022243333, 0.1835,  0.146583333  
;             698.3, 0.04047, 0.19548, 0.167933333  
;             896.8, 0.230096667, 0.190093333, 0.193353333  
;             1020.6,  0.226776667, 0.17783, 0.195176667  
;       ”    
    
 OUTPUT:  
; Fraction data (ENVI standard file): n+1 bands, the last band is residual error (RMSE).  

### Option 2: use the function of 'FCLS_unmixing' in IDL code  
;  SYNTAX  
;   Result=FCLS_unmixing(ImgData, endmember, residual=residual, end_names=end_names)  
;  
;  RETURN VALUE: Result data of 3-dimension matrix(ns×nl×n_endmember)  (n_endmember is number of endmembers, ns and nl are the size of input image)  
;  
;  INPUT:    
;   ImgData: 3-dimension matrix(ns×nl×nb)  (nb is the number of spectral bands)  
;   endmember: 2-dimension matrix (n×nb) of endmember spectra  
;    
； KEYWORDS:  
;   residual: 2-dimension matrix (ns×nl), residual of linear spectral mixture model  
;   end_names: names of endmember  
;  
## AUTHOR:  
; CHEN Xuehong (chenxuehong@bnu.edu.cn) @ Beijing Normal University, 2022-09.  
