pro calibresiduals, X,Y,datafit,Unit2,peak,err, x_peak, rasd, decsd, tipo, SNR, peak_cnt, err_cnt, Level, n_off, off, residual

  ;Minor section, to compute the rms-noise and the signal-to-noise ratio

  fit = LINFIT(X[0:20], Y[0:20], YFIT=datalinfit)
  subtracted = Y[0:20]-datalinfit
  mom = moment(subtracted)
  noise = sqrt(mom[1])
  residual = Y-datafit
  mom = moment(residual)
  resnoise = sqrt(mom[1])
  PRINTF, Unit2, "noise: ", noise, " counts; noise residuals: ", resnoise, " counts"
  ;  ;;print, "noise: ", noise, " counts; noise residuals: ", resnoise, " counts"
  SNR=peak/noise
  PRINTF, Unit2, "SNR: ", SNR
  peak_cnt = peak
  err_cnt = err
  Level = (Y[-1] + Y[0])/2.

  case tipo of
    0: begin
      meanlat=x_peak/180d*!dpi          ; radians
      off=(x_peak-decsd)                ; degrees
    end
    1: begin
      meanlat=decsd   ; radians
      off=(x_peak-rasd)*cos(meanlat)    ; degrees
    end
  endcase
  n_off=1

  return
end



pro getcnt2jy, flux, peak_cnt, err_cnt, tau0, datael, cnt2Jy, err_cnt2Jy, Unit2, Level
  
  ; Minor section, for the computation of conversion factor (counts --> Jy),
  ; inclusive of opacity correction. 
  
  peak_cnt=peak_cnt*exp(mean(tau0)/sin(datael)) 
  err_cnt=err_cnt*exp(mean(tau0)/sin(datael))
  cnt2Jy = flux/peak_cnt
  err_cnt2Jy = cnt2Jy*err_cnt/peak_cnt
  
  if flux gt 1000 then begin 
    cnt2Jy = -99.000
    err_cnt2Jy = -99.000
  endif  
  
  PRINTF, Unit2, " "
  PRINTF, Unit2, "cnt2Jy = (", cnt2Jy, "+-", err_cnt2Jy, " [Jy/cnt]"
  PRINTF, Unit2, " "
  PRINTF, Unit2, "Level = ", Level, " [cnt]"
  PRINTF, Unit2, " "
  return
end



pro calibfit, scanflag,stacflag,polyflag,section,tipo,allpath,namefile,Out3,flux,datael,tau0, xx,yy,ii,ff,x_mid,Nsamples,sd_sub,gpos,decsd,rasd, fwhm, n_off, off, p, e, c, d, SNR, plo, doplot


; Main procedure, devoted to the fitting operations, both in the linear and cubic cases, for targets.
;
; Authors: Marcello Giroletti, Simona Righini
; Last edited: Oct 5, 2017
;

  ; scanflag = scan buono (=1) o cattivo (=0)
  ; stacflag = single o stacked scan
  ; polyflag = gaussiana più grado del polinomio (primo='linear', terzo='cubic')
  ; section = Ch_0 o Ch_1
  ; tipo = flag di direzione - convertito qua sotto in direflag
  if (tipo eq 1) then direflag = 'RA' else direflag = 'DEC'
  
  ;  common logistics
  sep=path_sep()
  if scanflag eq 0 then begin
  
    n_off = 0
    off = 0.
    p = -99.
    e = -99.
    c = -99.
    d = -99.
    SNR = -99.
    if doplot eq 'y' then dataplot, 0, [0,1], [0,1], namefile, allpath, section, 0, 0, 0, 0, 0, direflag, polyflag, plo, stacflag
  
  endif else begin
    
    decs=decsd*!dpi/180.0
    media = (yy[ii]+yy[ff])/2.
    ; define new coordinates in which the source is centred at x=0 and the baseline is about 0 at the source position
    X=xx[ii:ff]-gpos
    Y=yy[ii:ff]-media

    if (polyflag eq 'linear') then begin
   
      datafit = GAUSSFIT(X, Y, A_fit, ESTIMATES=[(Y[x_mid]),0.0,sd_sub,(X[0]*Y[-1]-X[-1]*Y[0])/(X[0]-X[-1]),(Y[-1]-Y[0])/(X[-1]-X[0])], NTERMS=5, SIGMA=A_sig)
;      print, '                  Ampl(cnt)   Peak_pos(deg)      sigma(deg)          q(cnt)      m(cnt/deg)'
;      print, 'Estimates =',(Y[x_mid]),0.0,sd_sub,(X[0]*Y[-1]-X[-1]*Y[0])/(X[0]-X[-1]),(Y[-1]-Y[0])/(X[-1]-X[0])
;      print, 'Fit       =',A_fit
;      print, 'Sigma     =',A_sig
      if (tipo eq 0) then begin
        fwhm=2*SQRT(2*ALOG(2))*A_fit[2]*60.0  ; in arcmin, if fit is done on coordinates in degrees  ; *(speed/60.)*60. to be used if done in samples
      endif else begin
        fwhm=2*SQRT(2*ALOG(2))*A_fit[2]*60.0*cos(decs)  ; in arcmin, if fit is done on coordinates in degrees  ; *(speed/60.)*60. to be used if done in samples
      endelse
      x_peak=A_fit[1]+gpos ; coordinate corresponding to gaussian peak
      PRINTF, Out3, "     GAUSSIAN + LINEAR  "
      PRINTF, Out3, section+" Gaussian fit -- amplitude:", A_fit[0], " +- ", A_sig[0], " counts"
      PRINTF, Out3, section+" Gaussian fit -- sigma:    ", A_fit[2], " +- ", A_sig[2], " [s]"
      PRINTF, Out3, "    corresponding to FWHM:    ", fwhm, " [arcmin]"

      calibresiduals, X,Y,datafit,Out3,A_fit[0],A_sig[0], x_peak, rasd, decsd, tipo, SNR, peak_cnt, err_cnt, Level, n_off, off, residual
     
      ; computing and accounting for the RMS of the residuals in the overall amplitude error estimate
      nres=n_elements(residual)
      resrange=ceil(nres/5.0)   
      rescut=[residual[0:resrange],residual[-1*resrange,-1]]   ; avoiding the central part of the subscan, where artifacts can be present due to sidelobes
      resstat=moment(rescut)
      res_rms=sqrt(resstat[1])
      err_cnt=sqrt(err_cnt^2+res_rms^2+(0.03*peak_cnt)^2)  ; updated error for the amplitude measurement, including a default 3% uncertainty on calibrator flux-amplitude
      
      getcnt2Jy, flux, peak_cnt, err_cnt, tau0, datael, cnt2Jy, err_cnt2Jy, Out3, Level
      if doplot eq 'y' then dataplot, 1, X, Y, namefile, allpath, section, datafit, cnt2Jy, Level, SNR, residual, direflag, polyflag, plo, stacflag
 
    endif else begin

      ; estimates for fit with A0*exp(-(X-A1)^2/(2*A2))+A3+A4*X+A5*X^2+A6*X^3; ie cubic+gaussian
      A = [Y[x_mid],0.0,sd_sub,(X[0]*Y[-1]-X[-1]*Y[0])/(X[0]-X[-1]),(Y[-1]-Y[0])/(X[-1]-X[0]),0.0,0.0]
      A_GUESS = A
      ; compute the parameters without weights (DUM_W not set)
      datafit = CURVEFIT(X,Y,DUM_W, A, SIGMA, FITA=[1,1,1,1,1,1,1], FUNCTION_NAME='cal_gfunct', /DOUBLE, STATUS=suc_fit)
;      print, '                          A0=Ampl(cnt) A1=Peak_pos(deg)      A2=sigma(deg)        A3(cnt)     A4(cnt/deg)    A5(cnt/deg2)    A6(cnt/deg3)'
;      print, 'Function estimates:  ', A_GUESS
;      print, 'Function parameters: ', A
;      print, 'Sigma parameters:    ', SIGMA
      if (tipo eq 0) then begin
        fwhm=2*SQRT(2*ALOG(2))*A[2]*60.0  ; in arcmin, if fit is done on coordinates in degrees  ; *(speed/60.)*60. to be used if done in samples
      endif else begin
        fwhm=2*SQRT(2*ALOG(2))*A[2]*60.0*cos(decs)  ; in arcmin, if fit is done on coordinates in degrees  ; *(speed/60.)*60. to be used if done in samples
      endelse
      x_peak=A[1]+gpos ; coordinate corresponding to gaussian peak
      PRINTF, Out3, "     GAUSSIAN + CUBIC  "
      PRINTF, Out3, section+" Gaussian fit -- amplitude:", A[0], " +- ", SIGMA[2], " counts"
      PRINTF, Out3, section+" Gaussian fit -- sigma:    ", A[0], " +- ", SIGMA[0], " [s]"
      PRINTF, Out3, "    corresponding to FWHM:    ", fwhm, " [arcmin]"

      calibresiduals, X,Y,datafit,Out3,A[0],A[2], x_peak, rasd, decsd, tipo, SNR, peak_cnt, err_cnt, Level, n_off, off, residual
      
      ; computing and accounting for the RMS of the residuals in the overall amplitude error estimate
      nres=n_elements(residual)
      resrange=ceil(nres/5.0)
      rescut_b=[residual[0:resrange],residual[-1*resrange,-1]]   ; avoiding the central part of the subscan, where artifacts can be present due to sidelobes
      resstat=moment(rescut_b)
      res_rms=sqrt(resstat[1])
      err_cnt=sqrt(err_cnt^2+res_rms^2+(0.03*peak_cnt)^2)  ; updated error for the amplitude measurement
      
      getcnt2Jy, flux, peak_cnt, err_cnt, tau0, datael, cnt2Jy, err_cnt2Jy, Out3, Level
      if doplot eq 'y' then dataplot, 1, X, Y, namefile, allpath, section, datafit, cnt2Jy, Level, SNR, residual, direflag, polyflag, plo, stacflag
  
    endelse

    el_d=datael*180.0/!dpi
    p=peak_cnt
    e=err_cnt
    c=cnt2Jy
    d=err_cnt2Jy
 
  endelse
 
  return

end


