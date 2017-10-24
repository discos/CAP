pro runcalib, pickpath=pickpath, skypath=skypath, sub=sub, linear=linear, cubic=cubic, plot=plot, extlist=extlist, skipgc=skipgc

  ; This procedure is devoted to the reduction of cross scans acquired on flux calibrators and
  ; achieve the conversion factors cnt->Jy.
  ; It plots and analyses the data obtained by stacking all the subscans contained in any scan
  ; (i.e. subfolder) stored in the working folder.
  ; If users want to graphically choose the data folder, they must use:
  ;
  ;  IDL> runcalib, /pickpath    in order to select the path to the data parent folder (otherwise the chosen path is ./CALIBRATORS)
  ;  IDL> runcalib, /skypath     in order to select the path to the folder containing skydips (otherwise the chosen path is ./SKYDIPS)
  ;
  ; Measurements for all the single subscans can be obtained by explicitly choosing this option:
  ;
  ;  IDL> runcalib, /sub
  ;
  ; The gaussian fitting performed on the (sub)scans can be run using three different baseline-fitting options,
  ; in particular:
  ;
  ;  IDL> runcalib, /linear    performs, and outputs, only the linear fitting of baselines
  ;  IDL> runcalib, /cubic     performs, and outputs, only the cubic fitting of baselines
  ;  IDL> runcalib, /both      performs, and outputs, both the cubic and linear fitting of baselines
  ;
  ; By default, the choice is 'both'.
  ;
  ; At present, graphic output is exceedingly slow, in particular for single subscans (if the /sub option is selected).
  ; So it is disabled by default.
  ; To enable it, explicitly set the option:
  ;
  ;  IDL> runcalib, /plot
  ;
  ; If non-standard calibrators are to be considered, users must specify the /extlist option,
  ; which allows them to pick a proper external list, where source names and nominal flux
  ; densities (Jy) are provided in the simple format (A,D).
  ;
  ; The option
  ;
  ;  IDL> runcalib, /skipgc
  ;
  ; overrides the use of the standard gain curves; it poses Gain=1 for all the elevations,
  ; thus practically avoiding the elevation-dependent gain compensation.
  ; 
  ; COMPATIBILITY WITH CONVERTED FILES
  ; FITS files achieved with older versions of DISCOS (dating back to years 2008-2015)
  ; need to be converted into the proper FITS format (ask for procedure named updatefits.pro). 
  ; The program automatically handles the TP-like files resulting from the
  ; conversion of SARDARA acquisitions, which have raw counts levels in the
  ; orders of magnitude of 10E+06..10E+07, producing properly-formatted output tables.
  ;
  ; Authors: Marcello Giroletti, Simona Righini
  ; Last edited: Oct 18, 2017
  ;


  common logistics, workpath, calpath, sep, skydippath, myextlist, fitchoice, dosingle, doplot

  close, /all  ; to close possible pending text files

  sep=path_sep()

  if keyword_set (sub) then dosingle='y' else dosingle='n'

  if keyword_set(pickpath) then begin
    workpath=dialog_pickfile(/DIRECTORY, TITLE='Please pick PARENT FOLDER containing CALIBRATION SCANS')
  endif else begin
    CD, CURRENT=curr
    workpath=curr+sep+'CALIBRATORS'+sep
  endelse

  if keyword_set(skypath) then begin
    skydippath=dialog_pickfile(/DIRECTORY, TITLE='Please pick folder containing file Tau.txt')
  endif else begin
    CD, CURRENT=curr
    workpath=curr+sep+'SKYDIPS'+sep
  endelse

  if keyword_set(linear) then begin
    fitchoice='linear'
  endif else begin
    if keyword_set(cubic) then fitchoice='cubic' else fitchoice='both'
  endelse

  if keyword_set(plot) then doplot='y' else doplot='n'

  if keyword_set(extlist) then begin
    myextlist=dialog_pickfile(TITLE='Please pick EXTERNAL LIST of calibrators')
  endif else begin
    myextlist=' '
  endelse

  cal_stack, freq=freq
  return

end


PRO cal_gfunct, X, A, F, pder
  EX2 = A[0]*EXP(-(X-A[1])^2/(2*A[2]^2))
  F = EX2 + A[3] + A[4]*X +A[5]*X^2 +A[6]*X^3
  IF N_PARAMS() GE 4 THEN $
    pder = [[EX2/A[0]], [EX2*(X-A[1])/(A[2]^2)], [EX2*((X-A[1])^2)/(A[2]^3)], [replicate(1.0, N_ELEMENTS(X))], [X], [X^2], [X^3]]
END


function cal_wmean, val, dval, nan = nan, error = error
  compile_opt idl2
  on_error, 2

  ;- check inputs
  if n_params() ne 2 then begin
    return, !values.f_nan
  endif
  sz = n_elements(val)
  if n_elements(dval) ne sz then $
    message, 'val and dval must contain the same number of elements'


  ;- the calculation
  ;- requires at least one finite value
  if keyword_set(nan) then begin
    hit = where(finite(val) and finite(dval), ct)
    if ct eq 0 then begin
      error = !values.f_nan
      return, !values.f_nan
    endif
  endif

  dw = total(1D / dval^2, nan = keyword_set(nan))
  sum = total(1D * val / dval^2, nan = keyword_set(nan)) / dw
  error = sqrt(1D / dw)

  return, sum
end


pro cal_stack, path=path, out=out, plot=plot, beam=beam, speed=speed, dt=dt, source=source, flux=flux, freq=freq
  common logistics

  ;  !EXCEPT=2   ; enabling the indications of location/line numbers in exception messages
  !EXCEPT=0   ; silencing the printout of exception messages

  ; reads first file from first scan to set basic info in output files
  list = FILE_SEARCH(workpath+'20*', COUNT=number, /TEST_DIRECTORY)

  timings=dblarr(number+1)
  gap=dblarr(number+1)
  timings[0]=systime(1)
  gap[0]=0

  sublist = FILE_SEARCH(list[0]+sep+'20*.fits', COUNT=subnumber)
  ; XXX ToDO: ripetere la lettura per un file di ogni sub-folder e verificare che siano omogenei
  data=MRDFITS(sublist[0],4,/SILENT)
  RFinfo=MRDFITS(sublist[0],2,/SILENT)
  Sectinfo=MRDFITS(sublist[0],1,/SILENT)

  ; read frequency and bandwidth from RF table
  bw=RFinfo[0].bandWidth

  freq=RFinfo[0].frequency+bw/2.0
  dt=1.0/(Sectinfo[0].sampleRate*1e6)

  ; which telescope was used? Let's get to a 3-char code to be used afterwards
  firstmainh=mrdfits(sublist[0],0,firsthead0,/silent)
  findh0key, firsthead0, '*ANTENNA =*', keylab, sitename, infolab, siteflag
  splitsitename=strsplit(sitename,"'",/extract)
  cleansite=strupcase(strcompress(splitsitename[1],/remove_all))
  site=strmid(cleansite,0,3) ; this way, possible cases will be MED, SRT and NOT

  case site of
    'MED': begin
      beam=38.7/(freq/1000.0)  ; beam in arcmin
      ; stop
      ; frequency dependent normalized gain polynomials
      ; http://www.med.ira.inaf.it/ManualeMedicina/English/5.%20Efficiency.htm
      ; taken on 2015 September 17
      if (freq gt 4000 and freq lt 6000) then begin
        ; this is C band
        A_0=[7.9212981E-1,6.2403816E-3,-4.6834953E-5,0]
        A_1=[7.9212981E-1,6.2403816E-3,-4.6834953E-5,0]
      endif
      if (freq gt 8000 and freq lt 10000) then begin
        ; this is X band
        A_0=[0.730487,0.00773332,-5.54743e-05,0]   ; updated 2015, but it is not opacity-corrected
        A_1=[0.730487,0.00773332,-5.54743e-05,0]   ; updated 2015, but it is not opacity-corrected
       ; A_0=[6.1059261E-01,1.0623634E-02,-7.2457279E-05,0]  ; origin???
       ; A_1=[6.1059261E-01,1.0623634E-02,-7.2457279E-05,0]  ; origin???
      endif
      if (freq gt 18000 and freq le 18500) then begin
        ; this is K-low band
        A_0=[0.8373589,0.005140111,-4.058370E-5,0] ; as of Sept.30, 2015
        A_1=[0.8373589,0.005140111,-4.058370E-5,0] ; as of Sept.30, 2015
      endif
      if (freq gt 18500 and freq le 26500) then begin
        ; this is K-high band
        A_0=[0.7929185,0.005900533,-4.203179E-5,0] ; as of Sept.30, 2015
        A_1=[0.7929185,0.005900533,-4.203179E-5,0] ; as of Sept.30, 2015
      endif
    end
    'SRT': begin
      ; frequency dependent normalized gain polynomials
      ; as provided by Andrea Orlati (commissioning)
      if (freq gt 4000 and freq lt 8000) then begin
        ; this is C band
        beam=18.937/(freq/1000.0)  ; beam in arcmin
        A_0=[0.939813,0.00139,-0.000009,0] ; XXX To be updated
        A_1=[0.939813,0.00139,-0.000009,0] ; XXX To be updated
      endif
      if (freq ge 18000 and freq lt 27000) then begin
        ; this is K band
        beam=18.264/(freq/1000.0)  ; beam in arcmin
        if freq lt 21000 then begin
          A_0=[0.508145,0.014839,-0.000111,0] ; from official SRT page
          A_1=[0.508145,0.014839,-0.000111,0] ; from official SRT page
        endif
        if freq ge 21000 and freq le 23000 then begin
          A_0=[0.606733,0.01132,-0.000081,0] ; from official SRT page
          A_1=[0.606733,0.01132,-0.000081,0] ; from official SRT page
        endif
        if freq gt 23000 then begin
          A_0=[0.533025,0.015561,-0.000129,0] ; from official SRT page
          A_1=[0.533025,0.015561,-0.000129,0] ; from official SRT page
        endif
      endif
    end
    'NOT': begin
      ;beam=38.7/(freq/1000.0)  ; beam in arcmin - CLONED FROM MEDICINA XXX
      if (freq gt 4000 and freq lt 6000) then begin
        ; this is C band
        beam=7.5  ; fixed value!
       ; A_0=[0.98463659,0.00017053038,1.9348601e-09,0] ; old
       ; A_1=[0.98463659,0.00017053038,1.9348601e-09,0] ; old
        A_0=[0.95227062926,0.00198726460815,-2.0685473288e-05,0] ; Cassaro, Oct 16th 2017
        A_1=[0.95227062926,0.00198726460815,-2.0685473288e-05,0] ; Cassaro, Oct 16th 2017
      endif
      if (freq gt 8000 and freq lt 10000) then begin
        ; this is X band
        print, 'X Curve for Noto is not available, yet. Returning.'
        return
      endif
      if (freq gt 18000 and freq le 26500) then begin
        ; this is K band
        print, 'K Curve for Noto is not available, yet. Returning.'
        return
      endif
    end
    else: begin
      print, 'UNKNOWN SITE: cannot apply gain curves. Returning.'
      return
    end
  endcase
  beamd=beam/60.
  sd=beamd/(2*SQRT(2*ALOG(2)))


  ; open the main output files
  if fitchoice eq 'linear' or fitchoice eq 'both' then begin
    openw, Unit0, workpath+'cal_stack_lin_final.txt', /GET_LUN
    openw, Unit1, workpath+'cal_stack_lin_sep.txt', /GET_LUN
    printf, Unit0, FORMAT = '(A7,1X,A3)', "Site =",site
  endif
  if fitchoice eq 'cubic' or fitchoice eq 'both' then begin
    openw, Unit10, workpath+'cal_stack_cub_final.txt', /GET_LUN
    openw, Unit11, workpath+'cal_stack_cub_sep.txt', /GET_LUN
  endif

  if strmatch(sublist[0],'*TPlike*.*',/FOLD_CASE) then begin
    TPlike=1   ; this is a SARDARA file converted/integrated in TotalPower-like format
    if fitchoice eq 'linear' or fitchoice eq 'both' then begin
      printf, Unit0, '        Name    HHMMSS       # MJD    Freq   Elev         Amp_0        Amp_1      e_Amp_0      e_Amp_1        c2J_0        c2J_1      e_c2J_0      e_c2J_1'
      printf, Unit0, '                                       MHz    deg         count        count        count        count       Jy/cnt       Jy/cnt       Jy/cnt       Jy/cnt'
      printf, Unit1, '        Name    HHMMSS       # MJD    Freq   Elev         Amp_0        Amp_1      e_Amp_0      e_Amp_1        c2J_0        c2J_1      e_c2J_0      e_c2J_1    SNR_0    SNR_1   Offset    FW/HP   Type'
      printf, Unit1, '                                       MHz    deg         count        count        count        count       Jy/cnt       Jy/cnt       Jy/cnt       Jy/cnt                        deg                '
    endif
    if fitchoice eq 'cubic' or fitchoice eq 'both' then begin
      printf, Unit10, '        Name    HHMMSS       # MJD    Freq   Elev         Amp_0        Amp_1      e_Amp_0      e_Amp_1        c2J_0        c2J_1      e_c2J_0      e_c2J_1'
      printf, Unit10, '                                       MHz    deg         count        count        count        count       Jy/cnt       Jy/cnt       Jy/cnt       Jy/cnt'
      printf, Unit11, '        Name    HHMMSS       # MJD    Freq   Elev         Amp_0        Amp_1      e_Amp_0      e_Amp_1        c2J_0        c2J_1      e_c2J_0      e_c2J_1    SNR_0    SNR_1   Offset    FW/HP   Type'
      printf, Unit11, '                                       MHz    deg         count        count        count        count       Jy/cnt       Jy/cnt       Jy/cnt       Jy/cnt                        deg                '
    endif
  endif else begin
    TPlike=2   ; this is an original TotalPower FITS
    if fitchoice eq 'linear' or fitchoice eq 'both' then begin
      printf, Unit0, '        Name    HHMMSS       # MJD    Freq     Elev    Amp_0    Amp_1  e_Amp_0  e_Amp_1    c2J_0    c2J_1  e_c2J_0  e_c2J_1'
      printf, Unit0, '                                       MHz      deg    count    count    count    count   Jy/cnt   Jy/cnt   Jy/cnt   Jy/cnt'
      printf, Unit1, '        Name    HHMMSS       # MJD    Freq     Elev    Amp_0    Amp_1  e_Amp_0  e_Amp_1    c2J_0    c2J_1  e_c2J_0  e_c2J_1    SNR_0    SNR_1   Offset    FW/HP   Type'
      printf, Unit1, '                                       MHz      deg    count    count    count    count   Jy/cnt   Jy/cnt   Jy/cnt   Jy/cnt                        deg                '
    endif
    if fitchoice eq 'cubic' or fitchoice eq 'both' then begin
      printf, Unit10, '        Name    HHMMSS       # MJD    Freq     Elev    Amp_0    Amp_1  e_Amp_0  e_Amp_1    c2J_0    c2J_1  e_c2J_0  e_c2J_1'
      printf, Unit10, '                                       MHz      deg    count    count    count    count   Jy/cnt   Jy/cnt   Jy/cnt   Jy/cnt'
      printf, Unit11, '        Name    HHMMSS       # MJD    Freq     Elev    Amp_0    Amp_1  e_Amp_0  e_Amp_1    c2J_0    c2J_1  e_c2J_0  e_c2J_1    SNR_0    SNR_1   Offset    FW/HP   Type'
      printf, Unit11, '                                       MHz      deg    count    count    count    count   Jy/cnt   Jy/cnt   Jy/cnt   Jy/cnt                        deg                '
    endif
  endelse



  openw, Unit2, workpath+'cal_prints.txt', /GET_LUN


  ; reads skydip information
  ans=file_test(skydippath+'Tau.txt')
  if ans eq 1 then begin
    readcol, skydippath+'Tau.txt', tau_mjd, tau_freq, t30l, t90l, t30r, t90r, tau_empL, tau_empR, tau0L, tau0R, format='(d,d,d,d,d,d,d,d,d,d)', /SILENT
    for t=0,n_elements(tau_mjd)-1 do begin
      ; if for any reason tau_emp is considered more reliable than tau0, uncomment the following two lines
      ; tau0L[t]=tau_empL[t]
      ; tau0R[t]=tau_empR[t]
    endfor
    caltau = 1
  endif else begin
    print, ' '
    print, '+++++++++++++++++++++++++++++++++++++++++++++'
    print, 'Beware: no Tau.txt available, assuming tau0=0
    print, '+++++++++++++++++++++++++++++++++++++++++++++'
    wait, 0.5
    tau0L=[0,0]
    tau0R=[0,0]
    caltau = 0
  endelse

  ; definition of arrays
  p0=fltArr(2)
  p1=fltArr(2)
  e0=fltArr(2)
  e1=fltArr(2)
  c0=fltArr(2)                  ; count to jansky array
  c1=fltArr(2)
  d0=fltArr(2)                  ; error on count to jansky array
  d1=fltArr(2)
  SNR0=fltArr(2)
  SNR1=fltArr(2)
  offset=dblArr(2)
  p0c=fltArr(2)
  p1c=fltArr(2)
  e0c=fltArr(2)
  e1c=fltArr(2)
  c0c=fltArr(2)                  ; count to jansky array
  c1c=fltArr(2)
  d0c=fltArr(2)                  ; error on count to jansky array
  d1c=fltArr(2)
  SNR0c=fltArr(2)
  SNR1c=fltArr(2)
  offsetc=dblArr(2)
  fw_ratio0=fltarr(2)
  fw_ratio1=fltarr(2)
  fw_ratio0c=fltarr(2)
  fw_ratio1c=fltarr(2)


  ; scan-based analysis
  for i =0,number-1 do begin
    print, ' '
    print, '========================================='
    print, 'Now analysing scan ', list[i]
    print, '========================================='
    sommaL=0
    sommaR=0
    sommaLA = 0
    sommaRA = 0
    numLA = 0
    numRA = 0
    iiA = 0
    ffA = 0
    sommaLB = 0
    sommaRB = 0
    numLB = 0
    numRB = 0
    iiB = 0
    ffB = 0

    splitfoldname1=strsplit(list[i],sep,/extract)
    foldname=splitfoldname1[-1]
    splitfoldname2=strsplit(foldname,'-',/extract)
    hhmmss = splitfoldname2[1]
    sublist = FILE_SEARCH(list[i]+sep+'2*.fits', COUNT=subnumber)
    subscan=file_basename(sublist)

    ; reads flagging info
    err=0
    checkfile = list[i]+sep+'checkfile_'+foldname+'.txt'
    openr, Unit9, checkfile, ERROR=err, /GET_LUN
    flagcode=strarr(subnumber)
    flagcode[*]='+1cx'
    if (err lt 0) then begin
      print, "No valid flagging information for this scan"
      print, "Temporarily assuming all scans are good"
    endif else begin
      readcol, checkfile, F='a,a,a', subscan_t, flagcode_t, who, /SILENT
      if (n_elements(subscan_t) lt subnumber) then begin
        print, "Some subscans miss flagging information"
        print, "Temporarily assuming those subscans are good"
      endif
      for n_flag = 0, n_elements(subscan_t)-1 do begin
        for j_flag = 0, subnumber-1 do begin
          if (subscan_t[n_flag] eq subscan[j_flag]) then begin
            flagcode[j_flag]=flagcode_t[n_flag]
          endif
        endfor
      endfor
      close, Unit9
      free_lun, Unit9
    endelse


    if  dosingle eq 'y' then begin
      ; creating output files for single-subscan analysis
      if fitchoice eq 'linear' or fitchoice eq 'both' then begin
        openw, Unit3, list[i]+sep+'cal_single_linear.txt', /GET_LUN
        printf, Unit3,'        Name         # MJD    El(r)    El(d)    cnt_0    cnt_1    e_c_0    e_c_1    offL(d)  offR(d)    FW/HP scantype'
      endif
      if fitchoice eq 'cubic' or fitchoice eq 'both' then begin
        openw, Unit4, list[i]+sep+'cal_single_cubic.txt', /GET_LUN
        printf, Unit4,'        Name         # MJD    El(r)    El(d)    cnt_0    cnt_1    e_c_0    e_c_1    offL(d)  offR(d)    FW/HP scantype'
      endif
      openw, Unit5, list[i]+sep+'cal_prints.txt', /GET_LUN
    endif

    ; Dynamic arrays (for amplitudes and errors, for each section)

    Lampl=[ ]
    Rampl=[ ]
    Lerrs=[ ]
    Rerrs=[ ]
    offset_sub=[ ]
    type=[ ]
    mean_lat=[ ]

    ; subscan-based measurements
    for j=0,subnumber-1 do begin
      ; READ DATA TABLE (binary data table of the MegaTable FITS file)
      data=MRDFITS(sublist[j],4,/SILENT)
      RFinfo=MRDFITS(sublist[j],2,/SILENT)
      Sectinfo=MRDFITS(sublist[j],1,/SILENT)
      ; read frequency and bandwith from RF table
      bw=RFinfo[0].bandWidth
      freq=RFinfo[0].frequency+bw/2.0
      dt=1.0/(Sectinfo[0].sampleRate*1e6)

      ndat = (size(data.time))

      ; takes source name and subscantype from header
      mainh=mrdfits(sublist[j],0,head0,/silent)

      ; assessing whether this is an old converted FITS
      findh0key, head0, '*Converted FITS*', keylab, converted, infolab, convertedflag
      findh0key, head0, '*Obtained from ESCS 0.1 TP FITS*', keylab, version, infolab, oldestflag

      if convertedflag eq 0 then begin ; this IS NOT a converted file: ESCS 0.3 keyword names hold
        findh0key, head0, '*SOURCE*', keylab, target, infolab, tarflag
        findh0key, head0, '*SubScanType =*', keylab, orscantype, infolab, tarflag
        findh0key, head0, '*Azimuth Offset =*', keylab, hAZoff, infolab, typeflag ; offset utente e di sistema
        findh0key, head0, '*Elevation Offset =*', keylab, hELoff, infolab, typeflag ; offset utente e di sistema
        findh0key, head0, '*RightAscension Offset =*', keylab, hRAoff, infolab, typeflag ; offset utente e di sistema
        findh0key, head0, '*Declination Offset =*', keylab, hDECoff, infolab, typeflag ; offset utente e di sistema
        findh0key, head0, '*RightAscension =*', keylab, ras, infolab, typeflag   ; RA nominale del target, radianti
        findh0key, head0, '*Declination =*', keylab, decs, infolab, typeflag      ; DEC nominale del target, radianti
      endif else begin ; this IS a converted file: custom 8-char keyword names hold according to original FITS date, and offsets are not present
        findh0key, head0, '*SOURCE*', keylab, target, infolab, tarflag
        findh0key, head0, '*SUBSCANT=*', keylab, orscantype, infolab, tarflag
        hAZoff=0.0d
        hELoff=0.0d
        hRAoff=0.0d
        hDECoff=0.0d
        if oldestflag eq 0 then begin
          findh0key, head0, '*RightAscension =*', keylab, ras, infolab, typeflag   ; RA nominale del target, radianti
          findh0key, head0, '*Declination =*', keylab, decs, infolab, typeflag      ; DEC nominale del target, radianti
        endif else begin ; these are the oldest FITS, where several more keywords have converted names
          findh0key, head0, '*RIGHTASC=*', keylab, ras, infolab, typeflag   ; RA nominale del target, radianti
          findh0key, head0, '*DECLINAT=*', keylab, decs, infolab, typeflag      ; DEC nominale del target, radianti
        endelse
      endelse


      ; assigning labels and coefficients to calibrators
      target=strcompress(target, /remove_all)
      cleantarget=strsplit(target,"'",/extract)
      target=cleantarget[0]
      sourcename=target
      method='none'
      case strupcase(target) of
        '3C48': begin
          calflux, strupcase(target), freq/1000.0, intSnu, intSnu_err
          method='P&B'
        end
        '3C123': begin
          ; here follows a patch made necessary by a schedule mistake: 3C48 was observed instead of 3C123, using 3C123 as target name!
          if abs(ras*180.0/!dpi - 24.4221) lt 0.1 and abs(decs*180.0/!dpi - 33.1598) lt 0.1 then begin
            target='3C48'
            sourcename=target
          endif
          calflux, strupcase(target), freq/1000.0, intSnu, intSnu_err
          method='P&B'
        end
        '3C147': begin
          calflux, strupcase(target), freq/1000.0, intSnu, intSnu_err
          method='P&B'
        end
        '3C161': begin
          a=1.250
          b=+0.726
          c=-0.2286
          low=1408
          high=10550
          flux=10^(a+b*alog10(freq)+c*alog10(freq)*alog10(freq))
          method='Ott'
        end
        '3C286': begin
          calflux, strupcase(target), freq/1000.0, intSnu, intSnu_err
          method='P&B'
        end
        '3C295': begin
          calflux, strupcase(target), freq/1000.0, intSnu, intSnu_err
          method='P&B'
        end
        'NGC7027': begin
          calflux, strupcase(target), freq/1000.0, intSnu, intSnu_err
          method='P&B'
        end
        'DR21': begin
          a=1.810
          b=-0.122
          c=0.0
          low=7000
          high=31000
          flux=10^(a+b*alog10(freq)+c*alog10(freq)*alog10(freq))
          method='Ott'
        end
        'MYSRC': begin
          ; assuming flux density = 1 Jy
          a=0.00
          b=0.00
          c=0.00
          flux=10^(a+b*alog10(freq)+c*alog10(freq)*alog10(freq))
          method='CustomSRC'
        end
        'CAL01': begin
          ; fake source for tests, having flux density=0.1 Jy
          a=-1.00
          b=0.00
          c=0.00
          flux=10^(a+b*alog10(freq)+c*alog10(freq)*alog10(freq))
          method='CustomSRC'
        end
        'CAL05': begin
          ; fake source for tests, having flux density=0.5 Jy
          a=-0.301030
          b=0.00
          c=0.00
          flux=10^(a+b*alog10(freq)+c*alog10(freq)*alog10(freq))
          method='CustomSRC'
        end
        'CAL1': begin
          ; fake source for tests, having flux density=1 Jy
          a=0.00
          b=0.00
          c=0.00
          flux=10^(a+b*alog10(freq)+c*alog10(freq)*alog10(freq))
          method='CustomSRC'
        end
        'CAL2': begin
          ; fake source for tests, having flux density=2 Jy
          a=0.30103
          b=0.00
          c=0.00
          flux=10^(a+b*alog10(freq)+c*alog10(freq)*alog10(freq))
          method='CustomSRC'
        end
        'CAL3': begin
          ; fake source for tests, having flux density=3 Jy
          a=0.4771212
          b=0.00
          c=0.00
          flux=10^(a+b*alog10(freq)+c*alog10(freq)*alog10(freq))
          method='CustomSRC'
        end
        else: begin
          ; reading external list
          if myextlist ne ' ' then begin
            readcol, myextlist, extcalname, extcalflux, format='(A,D)', /SILENT
            foundit=where(strmatch(extcalname,target,/FOLD_CASE) eq 1)
            if foundit eq -1 then begin
              print, ' '
              print, target+' is an unknown calibrator, cannot retrieve a flux density to use:'
              print, 'it will return dummy values in the measurements'
              print, ' '
              a=4.00
              b=0.00
              c=0.00
              flux=10^(a+b*alog10(freq)+c*alog10(freq)*alog10(freq))
              method='CustomSRC'
            endif else begin
              a=alog10(extcalflux[foundit])   ; fixed value, as the external list does not contain polynomials
              b=0.00
              c=0.00
              flux=10^(a+b*alog10(freq)+c*alog10(freq)*alog10(freq))
              method='CustomSRC'
            endelse
          endif else begin
            print, ' '
            print, target+' is an unknown calibrator, cannot retrieve a flux density to use:'
            print, 'it will return dummy values in the measurements'
            print, ' '
            a=0.00
            b=0.00
            flux=10^(a+b*alog10(freq)+c*alog10(freq)*alog10(freq))
            method='CustomSRC'
          endelse
        endelse
      endcase

      ; aligning the P&B estimates to the general nomenclature
      if method eq 'P&B' then flux=intSnu

      if j eq 0 then begin
        ; expliciting the flux density computation result, only for the first subscan
        print, ' '
        print, '***********'
        print, sourcename, 'Flux density', flux, ' Jy at', freq/1000.0, 'GHz', format='(A,1X,A,1X,D6.3,1X,A,1X,D5.2,1X,A3)'
        print, 'Estimate by', method, format='(A,1X,A)'
        print, '***********'
      endif

      ; nominal target coordinates, converted in degrees
      rasd=ras*180.0d/!dpi
      decsd=decs*180.0d/!dpi

      elapsed=(double(data.time)-double(data[0].time))*24.0*3600.0  ; elapsed time from first sample
      duration=dt*ndat[1]/60.0 ; duration of the whole subscan (approx)
      ; RA and DEC span of the subscan
      deltaDec=abs(data[0].decj2000-data[-1].decj2000)/!dpi*180.0
      deltaRA=abs(data[0].raj2000-data[-1].raj2000)*cos(mean(data.decj2000))/!dpi*180.0

      ; determine speed (deg/min) along scan direction
      ; if type is AZ or EL, treats is as DEC or RA depending on which quantity spans a larger range
      splitscantype=strsplit(orscantype,"'",/extract)
      cleanscantype=strcompress(splitscantype[1], /remove_all) ; cleaning the string so that it become more manageable
      case cleanscantype of
        'DEC': begin
          if (data[0].decj2000 lt data[-1].decj2000) then scantype = 1 else scantype = 2
          speed=deltaDec/duration
        end
        'RA': begin
          if (data[0].raj2000 lt data[-1].raj2000) then scantype = 3 else scantype = 4
          speed=deltaRA/duration
        end
        'AZ': begin
          if (deltaDec gt deltaRA) then begin
            if (data[0].decj2000 lt data[-1].decj2000) then scantype = 1 else scantype = 2
            speed=deltaDec/duration
          endif
          if (deltaDec lt deltaRA) then begin
            if (data[0].raj2000 lt data[-1].raj2000) then scantype = 3 else scantype = 4
            speed=deltaRA/duration
          endif
        end
        'EL': begin
          if (deltaDec gt deltaRA) then begin
            if (data[0].decj2000 lt data[-1].decj2000) then scantype = 1 else scantype = 2
            speed=deltaDec/duration
          endif
          if (deltaDec lt deltaRA) then begin
            if (data[0].raj2000 lt data[-1].raj2000) then scantype = 3 else scantype = 4
            speed=deltaRA/duration
          endif
        end
        'UNKNOWN': begin
          if (deltaDec gt deltaRA) then begin
            if (data[0].decj2000 lt data[-1].decj2000) then scantype = 1 else scantype = 2
            speed=deltaDec/duration
          endif
          if (deltaDec lt deltaRA) then begin
            if (data[0].raj2000 lt data[-1].raj2000) then scantype = 3 else scantype = 4
            speed=deltaRA/duration
          endif
        end
      endcase

      ; check for spikes in X-band wide bandwidth data  XXX Perché farlo solo per la banda X?
      if (bw gt 500 and freq lt 10000) then begin
        for k=3, ndat[1]-4 do begin
          temp0=[data[k-3].ch0,data[k-2].ch0,data[k-1].ch0,data[k+1].ch0,data[k+2].ch0,data[k+3].ch0]
          temp1=[data[k-3].ch1,data[k-2].ch1,data[k-1].ch1,data[k+1].ch1,data[k+2].ch1,data[k+3].ch1]
          avech0=mean(temp0)
          avech1=mean(temp1)
          stdch0=stddev(temp0)
          stdch1=stddev(temp1)
          if ((data[k].ch0-avech0) gt 7*stdch0) then begin
            data[k].ch0=avech0
          endif
          if ((data[k].ch1-avech1) gt 7*stdch1) then begin
            data[k].ch1=avech1
          endif
        endfor
      endif

      ; reverse arrays for scans in decreasing RA or Dec
      if (scantype eq 2 or scantype eq 4) then begin
        ;        data=reverse(data)
        data.raj2000=reverse(data.raj2000)
        data.decj2000=reverse(data.decj2000)
        data.ch0=reverse(data.ch0) ; original line
        ; data.ch0=reverse(data.ch2) ; ABC
        data.ch1=reverse(data.ch1)
      endif

      ; computing the sample offset with respect to the source position
      xdec=(data.decj2000-decs)/!dpi*180.0
      xra=(data.raj2000-ras)*cos(mean(data.decj2000))/!dpi*180.0

      scanname = sublist[j]
      case flagcode[j] of
        '-1c0': begin
          k0=1
          k1=0
        end
        '-1c1': begin
          k0=0
          k1=1
        end
        '+0cx': begin
          k0=1
          k1=1
        end
        '+1cx': begin
          k0=1
          k1=1
        end
        '-1cx': begin
          k0=0
          k1=0
        end
      endcase

      ; shift scans of type 'Dec' so they align with first one; sets start and end value of array with data from every scan
      if (scantype eq 1 or scantype eq 2) then begin
        if ((numLA+numRA) eq 0) then begin
          dataA=data
          ndatA=ndat
          xdecA=xdec
          iiA=0
          if (k0+k1 gt 0) then ffA=ndatA[1]-1
          decoffA = min(xdecA,xdec0A,/ABSOLUTE) ; xdec0A is the number of sample at the source declination for scan A
        endif else begin
          decoff = min(xdec,xdec0,/ABSOLUTE)
          ; shift data by (xdec0-xdec0B) samples
          ndelta=xdec0-xdec0A
          if (ndelta gt 0) then begin
            iiA=iiA
            if (k0+k1 gt 0) then ffA=min([ffA,ndat[1]-1-ndelta])
            for k=iiA,ffA do begin
              if (k0+k1 gt 0) then data[k]=data[k+ndelta]
            endfor
          endif
          if (ndelta lt 0) then begin
            if (k0+k1 gt 0) then iiA=max([iiA,0-ndelta])
            ; here we have a problem: we would like to start from ffA and go down to iiA, since ffA can be filled by shifting "data" right. but data only has dimension ndat, not ndat-ndelta (ndelta is <0)
            ;            ffA=min([ffA,ndat[1]-1-ndelta])
            if (k0+k1 gt 0) then ffA=min([ffA,ndat[1]-1])
            for k=ffA,iiA,-1 do begin
              if (k0+k1 gt 0) then data[k]=data[k+ndelta]
            endfor
          endif
          if (ndelta eq 0) then begin
            if (k0+k1 gt 0) then ffA=min([ffA,ndat[1]-1])
          endif
        endelse
      endif

      ; shifts scans of type 'RA' so they align with first one; sets start and end value of array with data from every scan
      if (scantype eq 3 or scantype eq 4) then begin
        if ((numLB+numRB) eq 0) then begin
          dataB=data
          ndatB=ndat
          xraB=xra
          iiB=0
          if (k0+k1 gt 0) then ffB=ndatB[1]-1
          raoffB = min(xraB,xra0B,/ABSOLUTE) ; xra0B is the number of sample at the source right ascension for scan B
        endif else begin
          raoff = min(xra,xra0,/ABSOLUTE)
          ; shift data by (xra0-xra0B) samples
          ndelta=xra0-xra0B
          if (ndelta gt 0) then begin
            iiB=iiB
            if (k0+k1 gt 0) then ffB=min([ffB,ndat[1]-1-ndelta])
            for k=iiB,ffB do begin
              if (k0+k1 gt 0) then data[k]=data[k+ndelta]
            endfor
          endif
          if (ndelta lt 0) then begin
            if (k0+k1 gt 0) then iiB=max([iiB,0-ndelta])
            ; here we have a problem: we would like to start from ffB and go down to iiB, since ffB can be filled by shifting "data" right. but data only has dimension ndat, not ndat-ndelta (ndelta is <0)
            ;            ffB=min([ffB,ndat[1]-1-ndelta])
            if (k0+k1 gt 0) then ffB=min([ffB,ndat[1]-1])
            for k=ffB,iiB,-1 do begin
              if (k0+k1 gt 0) then data[k]=data[k+ndelta]
            endfor
          endif
          if (ndelta eq 0) then begin
            if (k0+k1 gt 0) then ffB=min([ffB,ndat[1]-1])
          endif
        endelse
      endif

      time_s = (data.time-data[round(ndat[1]/2)].time)*24.d0*3600.d0

      midscan=FIX(ndat[1]/2)
      beamspan= beam/speed ; duration of acquisition for one beamwidth
      Nsamples=ROUND(beamspan/dt) ; number of samples in one beamwidth

      el_d=data[midscan].el*180.0/!dpi

      g_0=A_0[0]+A_0[1]*el_d+A_0[2]*el_d^2+A_0[3]*el_d^3
      g_1=A_1[0]+A_1[1]*el_d+A_1[2]*el_d^2+A_1[3]*el_d^3

      if keyword_set(skipgc) then begin
        ; Overriding the gain computation and applying a flat gain curve (constant value = 1.0)
        g_0=1.0
        g_1=1.0
      endif

      ; actually sums counts only for scans flagged as good and sets inputs for single scan fit
      if ( scantype eq 1 or scantype eq 2) then begin
        ; if (k0 ne 0) then sommaLA = sommaLA + k0*data.ch2/g_0  ; ABC
        if (k0 ne 0) then sommaLA = sommaLA + k0*data.ch0/g_0 ; original line
        if (k1 ne 0) then sommaRA = sommaRA + k1*data.ch1/g_1
        numLA = numLA + k0*1
        numRA = numRA + k1*1
        ascissa=(data.decj2000-hDECoff)*180d/!dpi
        gpos=decsd
      endif else begin
        ; if (k0 ne 0) then sommaLB = sommaLB + k0*data.ch2/g_0 ; ABC
        if (k0 ne 0) then sommaLB = sommaLB + k0*data.ch0/g_0 ; original line
        if (k1 ne 0) then sommaRB = sommaRB + k1*data.ch1/g_1
        numLB = numLB + k0*1
        numRB = numRB + k1*1
        ascissa=(data.raj2000-hRAoff/cos(mean(data.decj2000)))*180d/!dpi
        gpos=rasd
      endelse


      if dosingle eq 'y' then begin

        if (fitchoice eq 'linear') or (fitchoice eq 'both') then begin
          ; fit single subscans with linear + gaussian
          Lfwhm=0
          Rfwhm=0
          ; yy0=data.ch2/g_0 ; ABC
          yy0=data.ch0/g_0  ; original line
          yy1=data.ch1/g_1
          if (scantype eq 1 or scantype eq 2) then tipo=0 else tipo=1
          calibfit, k0,'single','linear','Ch_0',tipo,list[i],subscan[j],Unit5,flux,data[0].el,tau0L,ascissa,yy0,0,ndat[1]-1,midscan,Nsamples,sd*(1.+tipo*(1./cos(decs)-1.)),gpos,decsd,rasd, Lfwhm, nL, offL, peak_cnt_0, err_cnt_0, dum3, dum4, dum5, psingle, doplot
          calibfit, k1,'single','linear','Ch_1',tipo,list[i],subscan[j],Unit5,flux,data[0].el,tau0R,ascissa,yy1,0,ndat[1]-1,midscan,Nsamples,sd*(1.+tipo*(1./cos(decs)-1.)),gpos,decsd,rasd, Rfwhm, nR, offR, peak_cnt_1, err_cnt_1, dum3, dum4, dum5, psingle, doplot
          if doplot eq 'y' then begin
            splitscan=strsplit(subscan[j],'.',/extract)
            scanroot=splitscan[0]
            if (j eq subnumber-1) then begin
              ;    psingle.save, list[i]+'_lin_single.pdf', /LANDSCAPE, RESOLUTION=300, /APPEND, /CLOSE
            endif else begin
              ;   psingle.save, list[i]+sep+scanroot+'_single.pdf', /LANDSCAPE, RESOLUTION=300, /APPEND
            endelse
          endif
          if TPlike eq 2 then begin
            printf, Unit3, FORMAT = '(a15,D12.5,9(" ",D08.4),i)', sourcename, data[midscan].time,$
              data[midscan].el, el_d, peak_cnt_0, peak_cnt_1, err_cnt_0, err_cnt_1, offL, offR, Lfwhm/beam, scantype
          endif else begin
            printf, Unit3, FORMAT = '(a15,1X,D12.5,1X,D08.4,1X,D08.4,1X,E12.6,1X,E12.6,1X,E12.6,1X,E12.6,1X,D08.4,1X,D08.4,1X,D08.4,1X,i)', sourcename, data[midscan].time,$
              data[midscan].el, el_d, peak_cnt_0, peak_cnt_1, err_cnt_0, err_cnt_1, offL, offR, Lfwhm/beam, scantype
          endelse
        endif

        if (fitchoice eq 'cubic') or (fitchoice eq 'both') then begin
          ; fit single subscans with cubic + gaussian
          Lfwhm=0
          Rfwhm=0
          ; yy0=data.ch2/g_0 ;ABC
          yy0=data.ch0/g_0  ; original line
          yy1=data.ch1/g_1
          if (scantype eq 1 or scantype eq 2) then tipo=0 else tipo=1
          calibfit, k0,'single','cubic','Ch_0',tipo,list[i],subscan[j],Unit5,flux,data[0].el,tau0L,ascissa,yy0,0,ndat[1]-1,midscan,Nsamples,sd*(1.+tipo*(1./cos(decs)-1.)),gpos,decsd,rasd, Lfwhm, nL, offL, peak_cnt_0, err_cnt_0, dum3, dum4, dum5, psingle, doplot
          calibfit, k1,'single','cubic','Ch_1',tipo,list[i],subscan[j],Unit5,flux,data[0].el,tau0R,ascissa,yy1,0,ndat[1]-1,midscan,Nsamples,sd*(1.+tipo*(1./cos(decs)-1.)),gpos,decsd,rasd, Rfwhm, nR, offR, peak_cnt_1, err_cnt_1, dum3, dum4, dum5, psingle, doplot
          ; NOTA: modificata la gestione del caso /sub; ora il salvataggio del PDF avviene direttamente dentro alla procedura dataplot,
          ; che salva un file per ogni chiamata.

          if TPlike eq 2 then begin
            printf, Unit4, FORMAT = '(a15,D12.5,9(" ",D08.4),i)', sourcename, data[midscan].time,$
              data[midscan].el, el_d, peak_cnt_0, peak_cnt_1, err_cnt_0, err_cnt_1, offL, offR, Lfwhm/beam, scantype
          endif else begin
            printf, Unit4, FORMAT = '(a15,1X,D12.5,1X,D08.4,1X,D08.4,1X,E12.6,1X,E12.6,1X,E12.6,1X,E12.6,1X,D08.4,1X,D08.4,1X,D08.4,1X,i)', sourcename, data[midscan].time,$
              data[midscan].el, el_d, peak_cnt_0, peak_cnt_1, err_cnt_0, err_cnt_1, offL, offR, Lfwhm/beam, scantype
          endelse
        endif

        if (nL+nR gt 0) then offset_sub = [offset_sub,(offL+offR)/(nL+nR)]

      endif

    endfor

    if dosingle eq 'y' then begin
      ; closing the output files and freeing memory
      if fitchoice eq 'linear' or fitchoice eq 'both' then begin
        close, Unit3
        FREE_LUN, Unit3
      endif
      if fitchoice eq 'cubic' or fitchoice eq 'both' then begin
        close, Unit4
        FREE_LUN, Unit4
      endif
      close, Unit5
      free_lun, Unit5, /force
    endif

    ; stacked analysis begins

    for tipo = 0,1 do begin
      ascissa=time_s ; default if all subscans are void

      ; set some inputs
      case tipo of
        0:  begin
          if (numLA+numRA gt 0) then begin
            if (numLA gt 0) then sommaL = sommaLA/numLA else sommaL = sommaLA
            if (numRA gt 0) then sommaR = sommaRA/numRA else sommaR = sommaRA
            numL = numLA
            numR = numRA
            ii=iiA
            ff=ffA
            ; tipo=0 means scan declination, ie scantype Dec, ie constant RA
            ascissa=(dataA.decj2000-hDECoff)*180d/!dpi
            xdiffA = min(ascissa-decsd,x_mid,/ABSOLUTE) ; x_mid contains the sample number nearest to the source nominal position
            gpos=decsd
          endif else begin
            numL=0
            numR=0
          endelse
        end
        1:  begin
          if (numLB+numRB gt 0) then begin
            if (numLB gt 0) then sommaL = sommaLB/numLB else sommaL = sommaLB
            if (numRB gt 0) then sommaR = sommaRB/numRB else sommaR = sommaRB
            numL = numLB
            numR = numRB
            ii=iiA
            ff=ffB
            ; tipo=1 means scan right asc., ie scantype RA, ie constant Dec
            ascissa=(dataB.raj2000-hRAoff/cos(mean(dataB.decj2000)))*180d/!dpi
            xdiffB = min(ascissa-rasd,x_mid,/ABSOLUTE) ; x_mid contains the sample number nearest to the source nominal position
            gpos=rasd
          endif else begin
            numL=0
            numR=0
          endelse
        end
      endcase

      PRINTF, Unit2, " "
      PRINTF, Unit2, "****** "+scanname+" ******"
      PRINTF, Unit2, FORMAT = '("* MJD = ", D014.8)', data[midscan].time
      PRINTF, Unit2, "* nsamples = ", ndat[1]
      PRINTF, Unit2, "* Beamspan = ", beamspan, " [s]"
      PRINTF, Unit2, "* Beam samples = ", Nsamples

      if fitchoice eq 'linear' or fitchoice eq 'both' then begin

        ; LINEAR FIT OF STACKED SCAN
        Lfwhm=0
        Rfwhm=0

        calibfit, numL,'stacked','linear','Ch_0',tipo,list[i],list[i],Unit2,flux,data[0].el,tau0L,ascissa,sommaL,ii,ff,x_mid,Nsamples,sd*(1.+tipo*(1./cos(decs)-1.)),gpos,decsd,rasd, Lfwhm, n_offL, offL, dum1, dum2, dum3, dum4, dum5, pstack, doplot
        p0[tipo]=dum1
        e0[tipo]=dum2
        c0[tipo]=dum3
        d0[tipo]=dum4
        SNR0[tipo]=dum5
        fw_ratio0[tipo]=Lfwhm/beam

        calibfit, numR,'stacked','linear','Ch_1',tipo,list[i],list[i],Unit2,flux,data[0].el,tau0R,ascissa,sommaR,ii,ff,x_mid,Nsamples,sd*(1.+tipo*(1./cos(decs)-1.)),gpos,decsd,rasd, Rfwhm, n_offR, offR, dum1, dum2, dum3, dum4, dum5, pstack, doplot
        p1[tipo]=dum1
        e1[tipo]=dum2
        c1[tipo]=dum3
        d1[tipo]=dum4
        SNR1[tipo]=dum5
        fw_ratio1[tipo]=Rfwhm/beam

        if doplot eq 'y' then begin
          if (tipo eq 1 and fitchoice eq 'linear') then begin
            pstack.save, list[i]+'_stacked.pdf', /LANDSCAPE, RESOLUTION=300, /APPEND, /CLOSE
          endif else begin
            pstack.save, list[i]+'_stacked.pdf', /LANDSCAPE, RESOLUTION=300, /APPEND
          endelse
        endif

        if (n_offL+n_offR gt 0) then offset[tipo]=(offL+offR)/(n_offL+n_offR) else offset[tipo] = -99.

        el_d=data[0].el*180.0/!dpi

        if tipo eq 0 then sclab='DEC'
        if tipo eq 1 then sclab='RA'

        if TPlike eq 2 then begin
          printf, Unit1, FORMAT = '(a15,a7,D12.5,d8.1,13(" ",D8.4),A7)', sourcename, hhmmss, data[0].time, freq, $
            el_d, p0[tipo],p1[tipo],e0[tipo],e1[tipo],c0[tipo],c1[tipo],d0[tipo],d1[tipo],SNR0[tipo], SNR1[tipo], offset[tipo], fw_ratio0[tipo], sclab
        endif else begin
          printf, Unit1, FORMAT = '(a15,a7,D12.5,d8.1,d8.4,8(" ",E12.6),4(" ",D8.4),A7)', sourcename, hhmmss, data[0].time, freq, $
            el_d, p0[tipo],p1[tipo],e0[tipo],e1[tipo],c0[tipo],c1[tipo],d0[tipo],d1[tipo],SNR0[tipo], SNR1[tipo], offset[tipo], fw_ratio0[tipo], sclab
        endelse

      endif

      if fitchoice eq 'cubic' or fitchoice eq 'both' then begin

        ; CUBIC FIT OF STACKED SCAN
        Lfwhm=0
        Rfwhm=0

        calibfit, numL,'stacked','cubic','Ch_0',tipo,list[i],list[i],Unit2,flux,data[0].el,tau0L,ascissa,sommaL,ii,ff,x_mid,Nsamples,sd*(1.+tipo*(1./cos(decs)-1.)),gpos,decsd,rasd, Lfwhm, n_offL, offL, dum1, dum2, dum3, dum4, dum5, pstack, doplot
        p0c[tipo]=dum1
        e0c[tipo]=dum2
        c0c[tipo]=dum3
        d0c[tipo]=dum4
        SNR0c[tipo]=dum5
        fw_ratio0c[tipo]=Lfwhm/beam

        calibfit, numR,'stacked','cubic','Ch_1',tipo,list[i],list[i],Unit2,flux,data[0].el,tau0R,ascissa,sommaR,ii,ff,x_mid,Nsamples,sd*(1.+tipo*(1./cos(decs)-1.)),gpos,decsd,rasd, Rfwhm, n_offR, offR, dum1, dum2, dum3, dum4, dum5, pstack, doplot
        p1c[tipo]=dum1
        e1c[tipo]=dum2
        c1c[tipo]=dum3
        d1c[tipo]=dum4
        SNR1c[tipo]=dum5
        fw_ratio1c[tipo]=Rfwhm/beam

        if doplot eq 'y' then begin
          if (tipo eq 1) then begin
            pstack.save, list[i]+'_stacked.pdf', /LANDSCAPE, RESOLUTION=300, /APPEND, /CLOSE
          endif else begin
            pstack.save, list[i]+'_stacked.pdf', /LANDSCAPE, RESOLUTION=300, /APPEND
          endelse
        endif

        if (n_offL+n_offR gt 0) then offsetc[tipo]=(offL+offR)/(n_offL+n_offR) else offsetc[tipo] = -99.

        el_d=data[0].el*180.0/!dpi

        if tipo eq 0 then sclab='DEC'
        if tipo eq 1 then sclab='RA'

        if TPlike eq 2 then begin
          printf, Unit11, FORMAT = '(a15,a7,D12.5,d8.1,13(" ",D8.4),A7)', sourcename, hhmmss, data[0].time, freq, $
            el_d, p0c[tipo],p1c[tipo],e0c[tipo],e1c[tipo],c0c[tipo],c1c[tipo],d0c[tipo],d1c[tipo],SNR0c[tipo], SNR1c[tipo], offsetc[tipo], Lfwhm/beam, sclab
        endif else begin
          printf, Unit11, FORMAT = '(a15,a7,D12.5,d8.1,d8.4,8(" ",E12.6),4(" ",D8.4),A7)', sourcename, hhmmss, data[0].time, freq, $
            el_d, p0c[tipo],p1c[tipo],e0c[tipo],e1c[tipo],c0c[tipo],c1c[tipo],d0c[tipo],d1c[tipo],SNR0c[tipo], SNR1c[tipo], offsetc[tipo], Lfwhm/beam, sclab
        endelse
      endif
    endfor

    ; If the pointing offsets are available in both directions, the amplitude measurement is compensated for them.
    ; Otherwise, the measurement is rejected.
    ; In this position also the FWHM check against HPBW - i.e. beam - is performed. A tolerance of FWHM/beam is here defined, the
    ; measurements outside the consequent boundaries are rejected as well.

    fw_tol=0.2  ; we accept measurements where FWHM is within 20% of the nominal HPBW

    if (fitchoice eq 'linear') or (fitchoice eq 'both') then begin

      if (offset[0] eq -99. or offset[1] eq -99. or (abs(fw_ratio0[0]-1.0) ge fw_tol and abs(fw_ratio0[1]-1.0) ge fw_tol) or (abs(fw_ratio1[0]-1.0) ge fw_tol and abs(fw_ratio1[1]-1.0) ge fw_tol)) then begin
        peak_cnt_0=-99.0
        peak_cnt_1=-99.0
        err_cnt_0=-99.0
        err_cnt_1=-99.0
        cnt2Jy_0=-99.0
        cnt2Jy_1=-99.0
        err_cnt2Jy_0=-99.0
        err_cnt2Jy_1=-99.0
      endif else begin
        ; riscalo le ampiezze e gli errori della polarizzazione 0 per l'offset incrociato
        p0[0]=p0[0]/exp(-(offset[1]*1.66/beamd)^2.)
        p0[1]=p0[1]/exp(-(offset[0]*1.66/beamd)^2.)
        e0[0]=e0[0]/exp(-(offset[1]*1.66/beamd)^2.)
        e0[1]=e0[1]/exp(-(offset[0]*1.66/beamd)^2.)
        ; riscalo le ampiezze e gli errori della polarizzazione 1 per l'offset incrociato
        p1[0]=p1[0]/exp(-(offset[1]*1.66/beamd)^2.)
        p1[1]=p1[1]/exp(-(offset[0]*1.66/beamd)^2.)
        e1[0]=e1[0]/exp(-(offset[1]*1.66/beamd)^2.)
        e1[1]=e1[1]/exp(-(offset[0]*1.66/beamd)^2.)

        peak_cnt_0 = cal_wmean(p0, e0, /nan)
        peak_cnt_1 = cal_wmean(p1, e1, /nan)
        w0=1/(e0/p0)
        w1=1/(e1/p1)
        err_cnt_0 = sqrt((w0[0]*e0[0])^2+(w0[1]*e0[1])^2)/(w0[0]+w0[1])
        err_cnt_1 = sqrt((w1[0]*e1[0])^2+(w1[1]*e1[1])^2)/(w1[0]+w1[1])
        ;  err_cnt_0 = stddev(p0, /nan)
        ;  err_cnt_1 = stddev(p1, /nan)

        ; riscalo i ct2Jy e gli errori della polarizzazione 0 per l'offset incrociato
        c0[0]=c0[0]*exp(-(offset[1]*1.66/beamd)^2.)
        c0[1]=c0[1]*exp(-(offset[0]*1.66/beamd)^2.)
        d0[0]=d0[0]*exp(-(offset[1]*1.66/beamd)^2.)
        d0[1]=d0[1]*exp(-(offset[0]*1.66/beamd)^2.)
        ; riscalo i ct2Jy e gli errori della polarizzazione 1 per l'offset incrociato
        c1[0]=c1[0]*exp(-(offset[1]*1.66/beamd)^2.)
        c1[1]=c1[1]*exp(-(offset[0]*1.66/beamd)^2.)
        d1[0]=d1[0]*exp(-(offset[1]*1.66/beamd)^2.)
        d1[1]=d1[1]*exp(-(offset[0]*1.66/beamd)^2.)

        cnt2Jy_0 = cal_wmean(c0, d0, /nan)
        cnt2Jy_1 = cal_wmean(c1, d1, /nan)
        w0=1/(d0/c0)
        w1=1/(d1/c1)
        err_cnt2Jy_0 = sqrt((w0[0]*d0[0])^2+(w0[1]*d0[1])^2)/(w0[0]+w0[1])
        err_cnt2Jy_1 = sqrt((w1[0]*d1[0])^2+(w1[1]*d1[1])^2)/(w1[0]+w1[1])

        ;         last check: when calibrators are not known, the resulting dummy cnt2Jy values are to be reset to -99.00
        ;         as the above offset compensation affects even them (and steers them a little bit from the dummy value)
        if (c0[0] lt 0) or (c0[1] lt 0) then begin
          cnt2Jy_0 = -99.0
          err_cnt2Jy_0 = -99.0
        endif
        if (c1[0] lt 0) or (c1[1] lt 0) then begin
          cnt2Jy_1 = -99.0
          err_cnt2Jy_1 = -99.0
        endif
      endelse


      if TPlike eq 2 then begin
        printf, Unit0, FORMAT = '(a15,a7,D12.5,d8.1,9(" ",D8.4))', sourcename, hhmmss, data[0].time, freq, $
          el_d, peak_cnt_0, peak_cnt_1,$
          err_cnt_0, err_cnt_1, cnt2Jy_0, cnt2Jy_1, $
          err_cnt2Jy_0, err_cnt2Jy_1;, SNR0, SNR1
      endif else begin
        printf, Unit0, FORMAT = '(a15,a7,D12.5,d8.1,D8.4,8(" ",E12.6))', sourcename, hhmmss, data[0].time, freq, $
          el_d, peak_cnt_0, peak_cnt_1,$
          err_cnt_0, err_cnt_1, cnt2Jy_0, cnt2Jy_1, $
          err_cnt2Jy_0, err_cnt2Jy_1;, SNR0, SNR1
      endelse

      ; note that SNR is not reported as it was not recalculated for the combined types
    endif


    if (fitchoice eq 'cubic') or (fitchoice eq 'both') then begin

      if (offsetc[0] eq -99. or offsetc[1] eq -99. or (abs(fw_ratio0c[0]-1.0) ge fw_tol and abs(fw_ratio0c[1]-1.0) ge fw_tol) or (abs(fw_ratio1c[0]-1.0) ge fw_tol and abs(fw_ratio1c[1]-1.0) ge fw_tol)) then begin
        peak_cnt_0=-99.0
        peak_cnt_1=-99.0
        err_cnt_0=-99.0
        err_cnt_1=-99.0
        cnt2Jy_0=-99.0
        cnt2Jy_1=-99.0
        err_cnt2Jy_0=-99.0
        err_cnt2Jy_1=-99.0
      endif else begin
        ; riscalo le ampiezze e gli errori della polarizzazione 0 per l'offset incrociato
        p0c[0]=p0c[0]/exp(-(offsetc[1]*1.66/beamd)^2.)
        p0c[1]=p0c[1]/exp(-(offsetc[0]*1.66/beamd)^2.)
        e0c[0]=e0c[0]/exp(-(offsetc[1]*1.66/beamd)^2.)
        e0c[1]=e0c[1]/exp(-(offsetc[0]*1.66/beamd)^2.)
        ; riscalo le ampiezze e gli errori della polarizzazione 1 per l'offset incrociato
        p1c[0]=p1c[0]/exp(-(offsetc[1]*1.66/beamd)^2.)
        p1c[1]=p1c[1]/exp(-(offsetc[0]*1.66/beamd)^2.)
        e1c[0]=e1c[0]/exp(-(offsetc[1]*1.66/beamd)^2.)
        e1c[1]=e1c[1]/exp(-(offsetc[0]*1.66/beamd)^2.)

        peak_cnt_0 = cal_wmean(p0c, e0c, /nan)
        peak_cnt_1 = cal_wmean(p1c, e1c, /nan)
        w0=1/(e0c/p0c)
        w1=1/(e1c/p1c)
        err_cnt_0 = sqrt((w0[0]*e0[0])^2+(w0[1]*e0[1])^2)/(w0[0]+w0[1])
        err_cnt_1 = sqrt((w1[0]*e1[0])^2+(w1[1]*e1[1])^2)/(w1[0]+w1[1])
        ;  err_cnt_0 = stddev(p0, /nan)
        ;  err_cnt_1 = stddev(p1, /nan)

        ; riscalo i ct2Jy e gli errori della polarizzazione 0 per l'offset incrociato
        c0c[0]=c0c[0]*exp(-(offsetc[1]*1.66/beamd)^2.)
        c0c[1]=c0c[1]*exp(-(offsetc[0]*1.66/beamd)^2.)
        d0c[0]=d0c[0]*exp(-(offsetc[1]*1.66/beamd)^2.)
        d0c[1]=d0c[1]*exp(-(offsetc[0]*1.66/beamd)^2.)
        ; riscalo i ct2Jy e gli errori della polarizzazione 1 per l'offset incrociato
        c1c[0]=c1c[0]*exp(-(offsetc[1]*1.66/beamd)^2.)
        c1c[1]=c1c[1]*exp(-(offsetc[0]*1.66/beamd)^2.)
        d1c[0]=d1c[0]*exp(-(offsetc[1]*1.66/beamd)^2.)
        d1c[1]=d1c[1]*exp(-(offsetc[0]*1.66/beamd)^2.)

        cnt2Jy_0 = cal_wmean(c0c, d0c, /nan)
        cnt2Jy_1 = cal_wmean(c1c, d1c, /nan)
        w0=1/(d0/c0)
        w1=1/(d1/c1)
        err_cnt2Jy_0 = sqrt((w0[0]*d0[0])^2+(w0[1]*d0[1])^2)/(w0[0]+w0[1])
        err_cnt2Jy_1 = sqrt((w1[0]*d1[0])^2+(w1[1]*d1[1])^2)/(w1[0]+w1[1])

        ; last check: when calibrators are not known, the resulting dummy cnt2Jy values are to be reset to -99.00
        ; as the above offset compensation affects even them (and steers them a little bit from the dummy value)
        if (c0[0] lt 0) or (c0[1] lt 0) then begin
          cnt2Jy_0 = -99.0
          err_cnt2Jy_0 = -99.0
        endif
        if (c1[0] lt 0) or (c1[1] lt 0) then begin
          cnt2Jy_1 = -99.0
          err_cnt2Jy_1 = -99.0
        endif

      endelse

      if TPlike eq 2 then begin
        printf, Unit10, FORMAT = '(a15,a7,D12.5,d8.1,9(" ",D8.4))', sourcename, hhmmss, data[0].time, freq, $
          el_d, peak_cnt_0, peak_cnt_1,$
          err_cnt_0, err_cnt_1, cnt2Jy_0, cnt2Jy_1, $
          err_cnt2Jy_0, err_cnt2Jy_1;, SNR0, SNR1
      endif else begin
        printf, Unit10, FORMAT = '(a15,a7,D12.5,d8.1,D8.4,8(" ",E12.6))', sourcename, hhmmss, data[0].time, freq, $
          el_d, peak_cnt_0, peak_cnt_1,$
          err_cnt_0, err_cnt_1, cnt2Jy_0, cnt2Jy_1, $
          err_cnt2Jy_0, err_cnt2Jy_1;, SNR0, SNR1
      endelse
      ; note that SNR is not reported as it was not recalculated for the combined types
    endif

  endfor

  ; closing all txt output files
  close, /ALL

  ; freeing logical units
  if fitchoice eq 'linear' or fitchoice eq 'both' then begin
    free_lun, Unit0, /force
    free_lun, Unit1, /force
  endif
  if fitchoice eq 'cubic' or fitchoice eq 'both' then begin
    free_lun, Unit10, /force
    free_lun, Unit11, /force
  endif

  free_lun, Unit2, /force


  print, ' '
  print, ' ********************************'
  print, ' ******* RUNCALIB is DONE *******'
  print, ' *** Next step is CJ2TIMELINE ***'
  print, ' ********************************'
  print, ' '

  return
end

