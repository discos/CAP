pro opacity, elmin=elmin, elmax=elmax, t0=t0, m=m, q=q, allplots=allplots

  ; *****************
  ;       INFO
  ; *****************
  ; 
  ; The program estimates Tau0 values from skydips (DISCOS FITS format). 
  ;
  ; At startup a dialog box appears: select the folder holding
  ; either the skydip FITS files or the subfolders containing them. 
  ;
  ; The fitting procedure exploits the following algorithm: 
  ; CURVEFIT, exponential fitting, with partial derivatives, 2 free parameters
  ;    Using a fixed Tatm, it estimates T0 e Tau0
  ;    The employed function is
  ;      Tsys - T0 = Tatm*(1- exp(-Tau0*x))
  ;    where x = airmass = 1/sin(El)
  ;
  ; By default, the program uses initial guess values for T0:
  ;   Tric = T0 = 9 K @ 4.5-7.7 GHz
  ;   Tric = T0 = 24 K @ 8.0-10.0 GHz
  ;   Tric = T0 = 30 K @ 18.0-26.5 GHz
  ;    (the observed freq is read insite the FITS file itself)
  ; Users can indicate a custom value for T0, by passing the argument t0=ddd.d
  ; in the command line. 
  ;
  ; Tatm is derived from the weather parameters recorded in the FITS content:
  ;   Tamb = average of the measurements stored in the file, converted to K
  ;   Tatm = 0.683*Tamb + 78
  ; this value is kept fixed within the fitting procedure. 
  ; If users want to use a different relation, they might do so by providing
  ; their own parameters (m and q) so that Tatm will be computed as
  ;   Tatm = m*Tamb +q 
  ;
  ; As concerns the zenithal opacity Tau0, the pre-fitting guess values
  ; are estimated by means of a simplified relation, deriving the opacity 
  ; through the Tsys values measured near 30° and 90° of Elevation: 
  ;   init_tau=alog(2/(1+sqrt(1-4*(Tsys30-Tsys90)/Tatm)))
  ; Such Tsys values are extracted from the Tant datastream inside the skydip,
  ; by selecting the nearest available samples.   
  ;
  ; By default, a skydip is considered to span between 20° and 87° of elevation.
  ; The elmin ed elmax optional parameters can be used when needing to specify custom ranges.
  ;
  ; OUTPUT:
  ; The program generates a Tau.txt file containig the Tau0 estimates, plus a skydip_fit.ps 
  ; file showing various plots in multi-plot pages.
  ; Using the /allplots option, each FITS will produce an on-screen plot, automatically 
  ; saved into a separate .jpg file.  
  ;
  ; EXAMPLE:
  ; IDL> opacity, elmin=25, elmax=85      exploits the skydip data in the 25°-85° elevation range only
  ; IDL> opacity, elmin=28                exploits the skydip data in the 28°-87° elevation range
  ; IDL> opacity, t0=40                   imposes that T0 is 40 K
  ; IDL> opacity, /allplots               produces individual plots (on screen and in JPG files)
  ; IDL> opacity, m=0.5, q=100            computes Tatm as Tatm = 0.5*Tamb + 100
  ;
  ; Last edited: Dec 6th, 2017

  ; ************************
  ;        PROGRAMMA
  ; ************************

  if not keyword_set(elmin) then elmin=20.0
  if not keyword_set(elmax) then elmax=87.0
  if not keyword_set(m) then m=0.683
  if not keyword_set(q) then q=78.0

  ; Resuming the Skydip files and extrapolating the Tau0 values
  print, ' '
  print, '**********************************************'
  print, 'Select folder with SKYDIPS files or subfolders'
  print, '**********************************************'
  print, ' '
  fin=dialog_pickfile(/DIRECTORY)
  mypath=fin

  close, /all
  openw, 5, fin+'Tau.txt'
  printf, 5, 'Folder: '+fin
  printf, 5, '********************************************************************************'
  printf, 5, 'MJD           Freq     Ts30L  Ts90L  Ts30R  Ts90R  empTauL empTauR  TauL    TauR'

  ; setup for PS file output
  set_plot,'ps'
  device,file=fin+'skydip_fit.ps',/landscape,/decomposed
  !p.multi=[0,2,2]

  ; Searching for the input files (includes subfolders)
  sky_list=file_search(fin,'2*.fits',count=num,/fully_qualify_path)

  ; Array for the data streams and output streams 
  timeline=dblarr(4,num)
  alltsys=dblarr(4,num)
  colore=lonarr(2)
  tausL=dblarr(num)
  tausR=dblarr(num)


  ; Cycling on the FITS files
  for i=0, num-1 do begin
    ; reading data and other parameters from the different FITS binary tables
    RFinfo=mrdfits(sky_list[i],2,/silent)
    sky_bandwidth=RFinfo[0].bandWidth
    freq=RFinfo[0].frequency+sky_bandwidth/2.0
    timeline(1,i)=freq
    data=mrdfits(sky_list[i],4,/silent)
    Xelall=data.el*180.0/!pi
    time=data.time
    ; average time for this acquisition
    timemom=moment(time)
    timeav=timemom[0]       ; average subscan MJD time
    timeline(0,i) = timeav
    
    elspan=abs(Xelall[-1]-Xelall[0])
    if elspan lt 40.0 then begin
      ; this is not a skydip. Maybe it is a misplaced FITS file. 
      print, 'File ',sky_list[i],' seems not to be a skydip. Skipping it.'
      continue
    endif

    ; Isolating the FITS filename 
    pathsep=path_sep()
    pathelems=strsplit(sky_list[i],pathsep,/extract)
    fitsname=pathelems[-1]

    ; reading weather parameters from FITS content
    weather=data.weather
    tmom = moment(weather[1])
    tamb = tmom[0]+273.15   ; average ground temperature
    tatm = m*tamb + q  ; extrapolating the atmospheric temperature - if no user-defined values were given, default parameters by F.Buffa (SRT) are used

    ; data streams in Tant (K)
    datak=mrdfits(sky_list[i],5,/silent)
    str0=datak.ch0
    str1=datak.ch1

    ; searching for the samples acquired nearest to 30 and 90 degrees of elevation
    around30=where(abs(Xelall-30.0) eq min(abs(Xelall-30.0)))
    around90=where(abs(Xelall-90.0) eq min(abs(Xelall-90.0)))
    Tsys30L=str0[around30]
    Tsys90L=str0[around90]
    Tsys30R=str1[around30]
    Tsys90R=str1[around90]

    ; empirical formula to pre-compute Tau0 values:
    ; we use these as initial guess values for the fitting procedure
    tau0L=alog(2/(1+sqrt(1-4*(Tsys30L-Tsys90L)/tatm)))
    tau0R=alog(2/(1+sqrt(1-4*(Tsys30R-Tsys90R)/tatm)))

    Xelall=reform(temporary(Xelall))
    elstart=where(abs(Xelall-elmin) eq min(abs(Xelall-elmin)))
    elstop=where(abs(Xelall-elmax) eq min(abs(Xelall-elmax)))
    ini=min([elstart[0],elstop[0]])
    fin=max([elstart[0],elstop[0]])
    range=fin-ini+1

    Xel=dblarr(range)
    Xel=Xelall[ini:fin]

    ; further storage arrays (maybe not all are subsequently exploited:
    ; some might be relics from older versions of the code. Pretty dirty coding, I know)
    streams=dblarr(range,2)
    fit=dblarr(range)
    fitb=dblarr(range)
    tau_channels=dblarr(2)
    mod1fit=dblarr(range,2)

    ; analysis, section by section (CH0 e CH1)
    for j=0, 1 do begin
      ch_str=strcompress(string(j), /remove_all)
      case j of
        0: begin
          streams(*,j)=str0[ini:fin]
          lab=fitsname
          datacol='000000'x
          fitcol='00ff00'x
          init_tau=tau0L
        end
        1: begin
          streams(*,j)=str1[ini:fin]
          lab=' '
          datacol='000000'x
          fitcol='0000ff'x
          init_tau=tau0R
        end
      endcase

      ; initial value for T0
      if not keyword_set(t0) then begin
        if (freq ge 4.5E+3) and (freq le 7.7E+3) then t0 = 9.0
        if (freq ge 8.0E+3) and (freq le 10.0E+3) then t0 = 24.0
        if (freq ge 18.E+3) and (freq le 26.5E+3) then t0 = 30.0
      endif

      am = 1./sin(Xel*!dpi/180.)  ; airmasses

      ; Duplication of variables, because the fitting procedure manipulates them 
      ; and updates their content
      tatm_0=tatm
      init_tau0=init_tau
      t0_0=t0
      init_par=[tatm,init_tau,t0]

      ; *** FITTING ***
      tau_mod1, j, am, streams(*,j), init_par, tau_val, fit, lab, thischi
      tau_channels(j)=tau_val
      mod1fit(*,j)=fit


      ; *** PLOTTING ***
      plot, Xelall, streams(*,j), ys=1, $
        title=lab, xtitle='Elevation [deg]', ytitle='Tsys', $
        charsize=chars, color=datacol
      xyouts, 60, max(streams(*,j))*0.9, 'CH'+ch_str
      oplot, Xelall, mod1fit(*,j), color=fitcol

      sp=path_sep()
      path=strsplit(sky_list[i],sp,/extract)
      file=strsplit(path[-1],'.',/extract)
      root=file[0]

      ; optional object graphics plots

      if keyword_set(allplots) then begin

        skyplot=plot(am, streams(*,j))
        skyplot.xtitle='Airmass'
        skyplot.ytitle='Tsys [K]'
        skyplot.title=root+'  (CH'+ch_str+')'
        skyplot.symbol='circle'
        skyplot.sym_size='0.5'
        skyplot.linestyle=' '
        fitplot=plot(am, fit, overplot=1)
        fitplot.color='red'
        fitplot.linestyle=' '
        fitplot.symbol='circle'
        fitplot.sym_size='0.3'

        ypos=max(streams(*,j))-(max(streams(*,j))-min(streams(*,j)))/4.0
        label=text(1.5,ypos,'tau0 = '+string(tau_val, format='(D5.3)'),/DATA, FONT_SIZE=14) ; per 20160211

        skyplot.save,  mypath+root+'_'+ch_str+'.jpg'
        wait, 1
        skyplot.Close

      endif

    endfor

    ; filling in the overall arrays (at present unemployed - leaving them for possible future updates)
    timeline(2,i)=tau_channels[0]
    timeline(3,i)=tau_channels[1]

    tausL[i]=tau0L
    tausR[i]=tau0R

    alltsys[0,i]=Tsys30L
    alltsys[1,i]=Tsys90L
    alltsys[2,i]=Tsys30R
    alltsys[3,i]=Tsys90R

    ; scrittura nell'output testuale
    if thischi lt 1.0 then begin
      printf, 5, timeav, freq, Tsys30L, Tsys90L, Tsys30R, Tsys90R, tau0L, tau0R, tau_channels[0], tau_channels[1], format='(D12.6,2X,D7.1,2X,D5.1,2X,D5.1,2X,D5.1,2X,D5.1,3X,D5.3,3X,D5.3,3X,D5.3,3X,D5.3)'
    endif else begin
      printf, 5, timeav, freq, Tsys30L, Tsys90L, Tsys30R, Tsys90R, tau0L, tau0R, tau_channels[0], tau_channels[1], 'Chi-squared=',thischi, format='(D12.6,2X,D7.1,2X,D5.1,2X,D5.1,2X,D5.1,2X,D5.1,3X,D5.3,3X,D5.3,3X,D5.3,3X,D5.3,1X,A12,D6.3)'
    endelse
  endfor

  device,/close
  close, /all
  print, ' '
  print, '*******************'
  print, '****** DONE! ******'
  print, '*******************'

  return
end


; FITTING MODEL: fixed Tatm, fitting data and extracting T0 and Tau0

PRO exptau_1, X, A, F, pder
  bx = EXP(-A[1]*X)
  F = A[0]*(1. - bx) + A[2]
  ; partial derivatives:
  pder = [[1.-bx], [A[0]*X*bx], [replicate(1.0, N_ELEMENTS(X))]]
END

pro tau_mod1, ch, X, Y, init_pars, tau_ch, unc_fit, filename, chisquared
  weights = 0*Y+1.0
  A=init_pars
  init_tatm=A[0]
  init_t0=A[2]
  ;unc_fit = CURVEFIT(X, Y, weights, A, SIGMA, FUNCTION_NAME='exptau_1', FITA=[1,1,1], tol=1.e-3, itmax=40, /double) ; use this if wanting to also fit for Tatm
  unc_fit = CURVEFIT(X, Y, weights, A, SIGMA, FUNCTION_NAME='exptau_1', FITA=[0,1,1], tol=1.e-3, itmax=40, /double)
  tau=A[1]
  tau_ch=tau
  tau_s = string(tau, format = '(D5.3)')
  chisq,Y,weights,unc_fit,n_elements(Y)-3,chi
  chisquared=chi
  if filename ne ' ' then begin
    print, '------------------------------------------------'
    print, 'FILE ', filename
    print, ' '
  endif
  print, '* FIT RESULTS for CH'+strcompress(string(ch),/remove_all)
  ;  print, '--> Fitted Tatm [K] = ', A[0], ' vs. init_Tatm ', init_tatm, format='(A,D6.2,A,D6.2)' ; use this if wanting to also fit for Tatm
  print, '--> Fixed Tatm [K] = ', init_tatm, format='(A,D6.2)'
  print, '--> Fitted T0 [K] = ', A[2], ' vs. initial T0 = ', init_t0, format='(A,D5.1,A,D5.1)'
  print, '--> Fitted tau0 = ', tau_s, format='(A,D5.3)'
  print, 'Chi-squared = ', chi
  print, ' '
end
