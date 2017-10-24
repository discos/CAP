pro calflux, namecal, nu, intSnu, intSnu_err

  ; Last Edited: Oct 6, 2017
  
  ; Adapted from SDI_CALFLUX_RELISE_V1: CREATION DATE 4/5/2016; VERSION 1 BY Elise Egron, Alberto Pellizzoni.
  ; Addition of Ott et al flux density estimate: Simona Righini - such part of the code is now muted, though
  ; P&B *interpolated* flux densities are provided. 
  
  ; EXAMPLE: calflux, '3C295', 7.2
  ; EXAMPLE (with outputs): calflux, '3C295', 7.2, Snu, intSnu, Snu_err, intSnu_err, pbchoice

  ; INPUTS:
  ; namecal = calibrator name (string)
  ; nu = required frequency (GHz)

  ; OUTPUTS:
  ;
  ; intSnu: flux density (Jy) at nu (error in intSnu_err )
  ; ("INTERPOLATION" based on polynomial coeff. from Perley et al. 2013)
  ;
  ; AVAILABLE CALIBRATORS; 3C123, 3C286, 3C295, 3C48, 3C147, NGC7027


  ; ******************************

  ; CALIBRATORS' DATA from Perley et al. 2013

  ; frequency table in GHz (33 entries)

  f_tab=[0.3275,1.015,1.275,1.465,1.865,2.565,3.565,4.535,4.835,4.885,6.135,6.885,7.465,8.435,8.485,8.735,11.06,12.890,14.635,14.715,14.915,14.965,17.422,18.230,18.485,18.585,20.485,22.460,22.835,24.450,25.836,26.485,28.450]


  ; CALIBRATOR DATA: fluxes (Jy), flux errors (Jy), polynomial coeff. (log Jy), polynomial coeff. errors (log Jy), RA/DEC coords. J2000 (deg)

  ; ****** 3C123 ******

  S_3C123=[145.0,66.2,46.6,47.8,38.7,28.9,21.4,16.9,16.0,15.88,12.81,11.20,11.01,9.20,9.10,8.86,6.73,6.035,5.34,5.02,5.132,5.092,4.272,4.181,4.090,3.934,3.586,3.297,3.334,2.867,2.697,2.716,2.436]

  S_3C123_err=[4.3,4.3,3.2,0.5,0.6,0.3,0.8,0.2,0.2,0.1,0.15,0.14,0.2,0.04,0.15,0.05,0.15,0.15,0.05,0.05,0.025,0.028,0.07,0.07,0.055,0.055,0.055,0.022,0.06,0.03,0.06,0.05,0.06]

  A_3C123=[1.8077d,-0.8018d,-0.1157d,0.d]
  A_3C123_err=[0.0036d,0.0081d,0.0047d,0.d]

  coord_3C123=[69.268230d,29.670505d]

  ; NOTES:
  ; missing flux value at 12.890 GHz: estimated (6.73+5.34)/2.=6.035
  ; missing flux value at 18.230 GHz: estimated (4.272+4.090)/2.=4.181
  ; missing error value at 12.890 GHz: estimated 0.15
  ; missing error value at 18.230 GHz: estimated 0.07
  ; coords. 04:37:04.3752, +29:40:13.818


  ; ****** 3C286 ******

  S_3C286=[26.1,18.4,13.8,15.0,13.2,10.9,9.5,7.68,7.33,7.297,6.49,5.75,5.70,5.059,5.045,4.930,4.053,3.662,3.509,3.375,3.399,3.387,2.980,2.860,2.925,2.880,2.731,2.505,2.562,2.387,2.181,2.247,2.079]

  S_3C286_err=[0.8,4.3,2.0,0.2,0.2,0.2,0.1,0.1,0.2,0.046,0.15,0.05,0.10,0.021,0.07,0.024,0.08,0.070,0.040,0.040,0.016,0.015,0.04,0.045,0.045,0.04,0.05,0.0160,0.05,0.03,0.06,0.05,0.05]

  A_3C286=[1.2515d,-0.4605d,-0.1715d,0.0336d]
  A_3C286_err=[0.0048d,0.0163d,0.0208d,0.0082d]

  coord_3C286=[202.784533d,30.509155d]

  ; NOTES:
  ; coords. 13:31:08.2879, +30:30:32.958


  ; ****** 3C295 ******

  S_3C295=[60.8,30.8,21.5,22.2,17.9,12.8,9.62,6.96,6.45,6.37,4.99,4.21,4.13,3.319,3.295,3.173,2.204,1.904,1.694,1.630,1.626,1.617,1.311,1.222,1.256,1.221,1.089,0.952,0.967,0.861,0.770,0.779,0.689]

  S_3C295_err=[1.8,7.3,3.0,0.5,0.3,0.2,0.2,0.09,0.15,0.04,0.05,0.05,0.07,0.014,0.05,0.016,0.05,0.04,0.04,0.03,0.008,0.007,0.025,0.05,0.020,0.015,0.015,0.005,0.015,0.020,0.02,0.0200,0.020]

  A_3C295=[1.4866d,-0.7871d,-0.3440d,0.0749d]
  A_3C295_err=[0.0036d,0.0110d,0.0160d,0.0070d]

  coord_3C295=[212.835279d, 52.202644d]

  ; NOTES:
  ; coords. 14:11:20.467, +52:12:09.52


  ; ****** 3C48 ******

  A_3C48=[1.3324d, -0.7690d, -0.1950d, 0.059d]
  A_3C48_err=[0.005d,0.01d,0.02d,0.02d]    ; NOT REPORTED IN PERLEY (Derived by ee, ap)

  coord_3C48=[24.422081d, 33.159759d]

  ; NOTES:
  ; polynomial coeff. for 2012 (table 11 p.13 Perley et al. 2013); ERRORS NOT AVAILABLE !
  ;
  ; coords. 01:37:41.2994, +33:09:35.132


  ; ****** 3C147 ******

  A_3C147=[1.4616d, -0.7187d, -0.2424d, 0.079d]
  A_3C147_err=[0.01d, 0.02d, 0.02d, 0.03d]    ; NOT REPORTED IN PERLEY (Derived by ee, ap)

  coord_3C147=[85.650575d, 49.852009d]

  ; NOTES:
  ; polynomial coeff. for 2012 (table 11 p.13 Perley et al. 2013); ERRORS NOT AVAILABLE !
  ;
  ; coords. 05:42:36.1379, +49:51:07.233


  ; ****** NGC7027 ******

  f_tab_NGC7027=[1.465, 4.885, 8.435, 14.965, 22.460, 43.340]  ; required frequency table

  S_NGC7027=[1.561d, 5.392d, 5.851d, 5.702d, 5.480d, 5.091d]
  S_NGC7027_err=[2., 5., 5., 7., 13., 42.]/1000d

  S_NGC7027_delta=[3.6, -3.0, -7.1, -7.4, -7.0, -5.0]/1000d

  S_NGC7027=S_NGC7027+S_NGC7027_delta*12.   ; VALID FOR 2012

  f00=5.5215-3.8*15./1000.
  f11=5.9294-7.12*15./1000.

  coord_NGC7027=[316.756638d, 42.236163d]

  ; NOTES:
  ; polynomial coeff. NOT AVAILABLE
  ; Assuming interpolated values for approx. 2012 (table 12 p.17 Perley et al. 2013) !
  ;
  ; coords. 21:07:01.5930, +42:14:10.186


  ;***************************


  S_name=['3C123','3C286','3C295', '3C48', '3C147', 'NGC7027']

  nn=n_elements(S_name)

  nfreq_tab=n_elements(f_tab)
  S_tot=fltarr(nn,2,nfreq_tab)
  S_tot(0,0,*)=S_3C123
  S_tot(1,0,*)=S_3C286
  S_tot(2,0,*)=S_3C295
  S_tot(0,1,*)=S_3C123_err
  S_tot(1,1,*)=S_3C286_err
  S_tot(2,1,*)=S_3C295_err


  if (strtrim(strupcase(namecal),2)) eq 'NGC7027' then begin
    f_tab=f_tab_NGC7027
    nfreq_tab=n_elements(f_tab_NGC7027)
    S_tot=fltarr(nn,2,nfreq_tab)
    S_tot(5,0,*)=S_NGC7027
    S_tot(5,1,*)=S_NGC7027_err
  endif


  A_tot=dblarr(nn,2,4)
  A_tot(0,0,*)=A_3C123
  A_tot(1,0,*)=A_3C286
  A_tot(2,0,*)=A_3C295
  A_tot(3,0,*)=A_3C48
  A_tot(4,0,*)=A_3C147
  A_tot(0,1,*)=A_3C123_err
  A_tot(1,1,*)=A_3C286_err
  A_tot(2,1,*)=A_3C295_err
  A_tot(3,1,*)=A_3C48_err
  A_tot(4,1,*)=A_3C147_err

  Coord_tot=dblarr(nn,2)
  Coord_tot(0,*)=coord_3C123
  Coord_tot(1,*)=coord_3C286
  Coord_tot(2,*)=coord_3C295
  Coord_tot(3,*)=coord_3C48
  Coord_tot(4,*)=coord_3C147
  Coord_tot(5,*)=coord_NGC7027

  wname=where(strtrim(strupcase(namecal),2) eq S_name,ncal)

  if ncal eq 0 then begin
    print
    print,'CALIBRATOR NOT FOUND!'
    print
    stop
  endif

  S=S_tot(wname,0,*)
  S_err=S_tot(wname,1,*)

  A=A_tot(wname,0,*)
  A_err=A_tot(wname,1,*)

  ra=Coord_tot(wname,0)
  dec=Coord_tot(wname,1)

  wmin=where(f_tab le nu)
  numin=max(f_tab(wmin))
  wmin=where(f_tab eq numin)

  wmax=where(f_tab ge nu)
  numax=min(f_tab(wmax))
  wmax=where(f_tab eq numax)

  Snu=(S(wmax)*(nu-numin)+S(wmin)*(numax-nu))/(numax-numin)
  Snu_err=(S_err(wmax)*(nu-numin)+S_err(wmin)*(numax-nu))/(numax-numin)
  Snu=Snu(0)
  Snu_err=Snu_err(0)

  intSnu_=A(0)+A(1)*alog10(nu)+A(2)*(alog10(nu))^2.+A(3)*(alog10(nu))^3.
  intSnu=(10.d0)^intSnu_
  intSnu_err=A_err(0)       ; TBD more precise error propagation ****** !!!!!! *******
  intSnu=intSnu(0)
  intSnu_err=intSnu_err(0)

  if not keyword_set(verb) then verb=1   ; for standalone usage of procedure, otherwise it fails
  
  if verb eq 1 then begin

;    print 
;    print,'**************************************************'
;    print,'CALIBRATOR: ', strtrim(namecal,2)
;    print,'**************************************************'
;    print,'>>> Fluxes by Perley et al. 2013'

;    if intSnu_ ne 0 then begin
;;      print, 'INTERPOLATED FLUX (Jy): ', strtrim(string(intSnu),2), ' +/- ', strtrim(string(intSnu_err),2)
;    endif
;
;    if Snu ne 0 then begin
;;      print, 'EXTRAPOLATED FLUX (Jy): ',strtrim(string(Snu),2),' +/- ',strtrim(string(Snu_err),2)
;    endif
;
;    if Snu ne 0 and intSnu_ ne 0 then begin
;;      print, 'DISCREPANCY (Jy): ', strtrim(string(abs(Snu-intSnu)),2)
;;      print, 'DISCREPANCY %: ', strtrim(string(abs(Snu-intSnu)/Snu*100.),2)
;      pbchoice='INT'
;    endif
;
;    if Snu eq 0 or Snu eq '-NaN' then begin
;;      print, 'EXTRAPOLATED FLUX NOT AVAILABLE!'
;      Snu=intSnu
;      Snu_err=intSnu_err
;      pbchoice='INT'
;    endif
;
    if A(0) eq 0 and A(1) eq 0 and A(2) eq 0 and A(3) eq 0 then begin
     ; print, 'INTERPOLATED FLUX NOT AVAILABLE! Assigning extrapolated one'
      intSnu=Snu
      intSnu_err=Snu_err
      pbchoice='EXT'
    endif


;    OTT ET AL. PART

;    print,'--------------------------------------------------'
;    print,'>>> Fluxes by Ott et al. 1994'

    case strupcase(namecal) of
      '3C48': begin
        a=2.465
        b=-0.004
        c=-0.125
        low=1.408
        high=23.780
      end

      '3C123': begin
        a=2.525
        b=0.246
        c=-0.1638
        low=1.408
        high=23.780
      end

      '3C147': begin
        a=2.806
        b=-0.140
        c=-0.1031
        low=0.500
        high=23.780
      end

      '3C286': begin
        a=0.956
        b=0.584
        c=-0.1644
        low=1.408
        high=43.200
      end

      '3C295': begin
        a=1.470
        b=0.765
        c=-0.2545
        low=1.408
        high=32.000
      end

      'NGC7027': begin
        a=1.322
        b=-0.134
        c=0.0
        low=10.550
        high=43.200
      end

    endcase

    ott_flux=10^(a+b*alog10(nu*1000.0)+c*alog10(nu*1000.0)*alog10(nu*1000.0))
    
    inrange='Y'
    ; warning on improper use of polynomial, i.e. when observed frequency is outside the validity range
    if (nu lt low) or (nu gt high) then begin
;      print, 'WARNING: requested frequency lies 
;      print, '         outside Ott et al. range
      inrange='N'
    endif
    

;    print, 'OTT INTERPOLATED FLUX (Jy): ', strtrim(string(ott_flux,format='(D7.4)'),2)
;    print, 'DISCREPANCY WITH PERLEY (Jy): ', strtrim(string(abs(ott_flux-intSnu),format='(D7.4)'),2)
;    print, 'DISCREPANCY WITH PERLEY %: ', strtrim(string(abs(ott_flux-intSnu)/ott_flux*100.,format='(D5.2)'),2)
;    print,'**************************************************'
;    print

  endif

  ;stop



end



