pro c2jtimeline, pickpath=pickpath, range=range, linf=linf, detplot=detplot

  ; This procedure performs an averaging or a linear fitting of the cnt2Jy values, considering the
  ; lists produced by runcalib.pro and grouping data according to a time range defined by the user (expressed in hours).
  ; By default, this gap is set to 24 hours.
  ;
  ; Users must specify if they want to extract linearly fitted cnt2Jy values (linfit value inside each gap),
  ; by setting the /linf option:
  ;
  ; IDL> c2jtimeline, /linf
  ;
  ; By default, runcalib's output files are searched for inside the working folder.
  ; Users might graphically select a different folder by using the /pickpath option:
  ;
  ; IDL> c2jtimeline, /pickpath
  ;
  ; c2jtimeline's detailed results are plotted into .ps and jpg files, one devoted to the EXLIN results, one to the EXCUB results.
  ; On screen, by default only the whole dataset/timeline is shown. To also display on screen the detailed plots,
  ; each showing the cnt2Jy trends inside a single time interval, use:
  ;
  ; IDL> c2jtimeline, /detplot
  ;
  ; Notice: the program automatically handles the measurements achieved on TP-like files resulting from the
  ; conversion of SARDARA acquisitions, which have raw counts levels in the
  ; orders of magnitude of 10E+06, 10E+07, producing properly-formatted output tables.
  ;
  ; Authors: Marcello Giroletti, Simona Righini
  ; Last edited: Jan 20, 2022
  ;

  sep=path_sep()

  if keyword_set(pickpath) then begin
    workpath=dialog_pickfile(/DIRECTORY)
  endif else begin
    cd, current=curr
    workpath=curr+sep+'CALIBRATORS'+sep
  endelse


  ;  if not keyword_set(range) then range=0.08333333       ;  5 minutes
  ;  if not keyword_set(range) then range=0.16666666       ; 10 minutes
  ;  if not keyword_set(range) then range=0.25             ; 15 minutes
  ;  if not keyword_set(range) then range=0.50             ; 30 minutes
  ;  if not keyword_set(range) then range=1.00             ; 1 hour
  if not keyword_set(range) then range=24.0                ; 24 hours


  if keyword_set(linf) then mode='linf' else mode='aver'

  infiles=strarr(2)
  outfiles=strarr(2)
  outplot=strarr(2)
  infiles[0]='cal_stack_lin_final.txt'
  ; infiles[0]='cal_stack_cub_final.txt' ; temporary substitution for 2010-2012 ra-only data
  infiles[1]='cal_stack_cub_final.txt'
  outfiles[0]='cnt2jy_exlin.txt'
  outfiles[1]='cnt2jy_excub.txt'
  outplot[0]='cnt2Jy_lin.ps'
  outplot[1]='cnt2Jy_cub.ps'

  sym=['circle','square','triangle','star','diamond','asterisk','X','plus','circle','square','triangle','star','diamond','asterisk','X','plus']
  colors=['blue','red','green','dark orange','plum','gold','saddle brown','hot pink','green yellow','medium orchid','midnight blue','firebrick','medium sea green','dark slate gray','teal','dark magenta']


  for f=0, 1 do begin

    checkinput=file_search(workpath+infiles[f],count=fin,/test_read)
    if fin eq 0 then begin
      print, ' '
      print, '***************************'
      print, 'File ',infiles[f],' not found in folder ', workpath
      print, '***************************'
      break
    endif

    ; read and fits the c2Jy from the calibrator scans as a constant or as a function of time

    readcol,  workpath+infiles[f], format='(a12,a7,D12.5,d8.1,d8.4,d8.4,d8.4,d8.4,d8.4,d8.4,d8.4,d8.4,d8.4)', name, hms, time, freq, el_d, counts_0, counts_1, err_cnt_0, err_cnt_1, cnt2Jy_0_el, cnt2Jy_1_el, err_cnt2Jy_0_el, err_cnt2Jy_1_el, /SILENT, count=ngood
    if ngood eq 0 then begin
      readcol,  workpath+infiles[f], format='(a12,a7,D12.5,d8.1,d8.4,E12.6,E12.6,E12.6,E12.6,E12.6,E12.6,E12.6,E12.6)', name, hms, time, freq, el_d, counts_0, counts_1, err_cnt_0, err_cnt_1, cnt2Jy_0_el, cnt2Jy_1_el, err_cnt2Jy_0_el, err_cnt2Jy_1_el, /SILENT, count=ngood
    endif

    day=strcompress(string(long(time[0])),/remove_all)
    t_h=24*(time-floor(time[0]))
    timespan=(time[-1]-time[0])*24.0d
    elapsed=(time-time[0])*24.0

    num=n_elements(name)
    intnum=ceil(timespan/range) ; number of time intervals with "gap" length

    if timespan eq 0.0 then intnum=1   ; to handle cases where only one calibration measurement is available


    ; finding the number of unique calibrator names
    nameu = name[UNIQ(name, SORT(name))]
    numcals=n_elements(nameu)

    ; setup for the output text and graphic file
    openw, Unit, workpath+outfiles[f], /GET_LUN
    printf, Unit, FORMAT = '(a,a)','Selected interpolation mode = ',mode
    printf, Unit, FORMAT = '(a,a,a)','Time range = ',strcompress(string(range,format='(D6.3)'),/remove_all),' hours'
    printf, Unit, FORMAT = '(a,a)', 'Number of calibration solution intervals = ',strcompress(string(intnum,format='(I3)'),/remove_all)

    if counts_0[0] le 3000 then begin
      ; printf, Unit, FORMAT = '(a2,a9,2(a11),9(a11))',"n","Freq","t_i","t_f","m_0","q_0","e_m_0","e_q_0","m_1","q_1","e_m_1","e_q_1"
      printf, Unit, '   n     Freq        t_i        t_f          m_0          q_0        e_m_0        e_q_0          m_1          q_1        e_m_1        e_q_1'
    endif else begin
      printf, Unit, ' n     Freq        t_i        t_f          m_0          q_0        e_m_0        e_q_0          m_1          q_1        e_m_1        e_q_1'
    endelse
    loadct, 13, /silent, ncolors=numcals

    set_plot,'ps'
    device, /COLOR, bits=8
    device, file=workpath+outplot[f], /COLOR, /landscape


    if f eq 0 then origin='EXLIN' else origin='EXCUB'
    
    real0=where(cnt2Jy_0_el gt 0)
    real1=where(cnt2Jy_1_el gt 0)
    nreal0=n_elements(real0)
    nreal1=n_elements(real1)
    c2j0stats=moment(cnt2Jy_0_el[real0])
    c2j1stats=moment(cnt2Jy_1_el[real1])
    c2j0av=c2j0stats[0]
    c2j1av=c2j1stats[0]
    c2j0stdev=sqrt(c2j0stats[1])
    c2j1stdev=sqrt(c2j1stats[1])
    
    wholeplot=errorplot(el_d[real0], cnt2Jy_0_el[real0], el_d[real0]*0, err_cnt2Jy_0_el[real0], name='CH0')
    wholeplot.title='All cnt2Jy measurements - '+origin
    wholeplot.xtitle='Elevation (deg)'
    wholeplot.ytitle='Jy/cnt'
    wholeplot.linestyle=' '
    wholeplot.color='midnight blue'
    wholeplot.symbol='square'
    wholeplot.sym_filled=1
    wholeplot.sym_size=1.0
    wholeplot.xrange=[0,90]
    wholeplot2=errorplot(el_d[real1], cnt2Jy_1_el[real1], el_d[real1]*0, err_cnt2Jy_1_el[real1],/overplot, name='CH1')
    wholeplot2.linestyle=' '
    wholeplot2.color='deep sky blue'
    wholeplot2.symbol='triangle'
    wholeplot2.sym_filled=1
    wholeplot2.sym_size=1.0
    wholeplot3=plot([min(el_d[real0]),max(el_d[real0])], [c2j0av,c2j0av],/overplot, name='Aver0')
    wholeplot3.linestyle='-'
    wholeplot3.color='midnight blue'
    wholeplot4=plot([min(el_d[real1]),max(el_d[real1])], [c2j1av,c2j1av],/overplot, name='Aver1')
    wholeplot4.linestyle='-'
    wholeplot4.color='deep sky blue'
    wholeplot5=plot([min(el_d[real0]),max(el_d[real0])], [c2j0av+c2j0stdev,c2j0av+c2j0stdev],/overplot, name='StDev0')
    wholeplot5.linestyle='dash_dot'
    wholeplot5.color='midnight blue'
    wholeplot6=plot([min(el_d[real0]),max(el_d[real0])], [c2j0av-c2j0stdev,c2j0av-c2j0stdev],/overplot)
    wholeplot6.linestyle='dash_dot'
    wholeplot6.color='midnight blue'
    wholeplot7=plot([min(el_d[real1]),max(el_d[real1])], [c2j1av+c2j1stdev,c2j1av+c2j1stdev],/overplot, name='StDev1')
    wholeplot7.linestyle='dash_dot'
    wholeplot7.color='deep sky blue'
    wholeplot8=plot([min(el_d[real1]),max(el_d[real1])], [c2j1av-c2j1stdev,c2j1av-c2j1stdev],/overplot)
    wholeplot8.linestyle='dash_dot'
    wholeplot8.color='deep sky blue'
    

    leg = legend(TARGET=[wholeplot, wholeplot2], POSITION=[0.9,0.9], /AUTO_TEXT_COLOR)
    
    wholeplot.save, workpath+'timeline_'+origin+'.jpg'


    ; overall plot of full dataset

    !p.multi=[1,1,0]
    plot, elapsed, cnt2Jy_0_el, ys=1, xstyle=1, psym=2, ytitle="cnt2Jy", xtitle="Elapsed time [hh.h] since "+day
    ; plot, elapsed, cnt2Jy_0_el, yrange=[0,max(cnt2Jy_0_el)+1], xstyle=1, psym=2, ytitle="cnt2Jy", xtitle="Elapsed time [hh.h] since "+day
    oplot, elapsed, cnt2Jy_1_el, psym=6


    ; preparing multi-plot page for details
    !p.multi=[0,2,2]

    for j=0, intnum-1 do begin

      thisgap=where(elapsed ge elapsed[0]+range*j and elapsed le elapsed[0]+range*(j+1))  ; these are the measurements contained in the time interval at hand

      if thisgap[0] ne -1 then begin

        ti=time[thisgap[0]]
        tf=time[thisgap[-1]]

        index0 = counts_0[thisgap]
        index1 = counts_1[thisgap]
        t_0=time[thisgap]
        t_1=time[thisgap]
        el_0=el_d[thisgap]
        el_1=el_d[thisgap]
        c2j_0=cnt2Jy_0_el[thisgap]
        c2j_1=cnt2Jy_1_el[thisgap]
        err_c2j_0=err_cnt2Jy_0_el[thisgap]
        err_c2j_1=err_cnt2Jy_1_el[thisgap]
        C_0=dblarr(2)  ; for actual interpolation, to be written in output file
        C_1=dblarr(2)  ; for actual interpolation, to be written in output file
        C_sig0=dblarr(2)  ; for actual interpolation, to be written in output file
        C_sig1=dblarr(2)  ; for actual interpolation, to be written in output file
        Ch_0=dblarr(2)    ; for plot only, in elapsed time
        Ch_1=dblarr(2)    ; for plot only, in elapsed time
        Ch_sig0=dblarr(2) ; for plot only, in elapsed time
        Ch_sig1=dblarr(2) ; for plot only, in elapsed time


        ; selecting only the meaningful values (discarding the "-99.00" cases)

        valid0=WHERE(index0 GT 0, valnum0)
        valid1=WHERE(index1 GT 0, valnum1)

        if valnum0 ge 1 then begin
          t_0val=dblarr(valnum0)
          el_0val=dblarr(valnum0)
          c2j_0va=dblarr(valnum0)
          err_c2j_0val=dblarr(valnum0)
          th_0val=dblarr(valnum0)
        endif

        if valnum1 ge 1 then begin
          t_1val=dblarr(valnum1)
          el_1val=dblarr(valnum1)
          c2j_1va=dblarr(valnum1)
          err_c2j_1val=dblarr(valnum1)
          th_1val=dblarr(valnum1)
        endif

        e=findgen(1000) ; for elevation axis representation
        e=e/1000.0*90

        if valid0[0] ne -1 then begin

          ; there are usable measurements
          t_0val=t_0[valid0]
          el_0val=el_0[valid0]
          c2j_0val=c2j_0[valid0]
          err_c2j_0val=err_c2j_0[valid0]
          th_0val=elapsed[valid0]

          if valnum0 gt 1 then begin
            ; there are at least two usable measurements
            w_0=1/(err_c2j_0val/c2j_0val)^2
            meanerr,c2j_0val,err_c2j_0val,x_0mean,sigma_0m,sigma_0d
            C_0 = LINFIT( t_0val, c2J_0val, CHISQ=chisq, /DOUBLE, MEASURE_ERRORS=err_c2J_0val, PROB=lik_0, SIGMA=C_sig0 )
            Ch_0 = LINFIT(th_0val, c2J_0val, CHISQ=chisq_h, /DOUBLE, MEASURE_ERRORS=err_c2J_0val_h, PROB=lik_0_h, SIGMA=Ch_sig0)
            el0plot=el_0val
            c2J0plot=c2j_0val
            errc2j0plot=err_c2j_0val
            th0plot=th_0val
          endif else begin
            ;only one measurement is available: linear fitting cannot be performed!
            x_0mean=c2j_0val[0]
            sigma_0d=err_c2j_0val[0]
            sigma_0m=err_c2j_0val[0]
            C_0[0]=c2j_0val[0]
            C_0[1]=0.0
            C_sig0[0]=err_c2j_0val[0]
            C_sig0[1]=0.0
            Ch_0[0]=c2j_0val[0]
            Ch_0[1]=0.0
            Ch_sig0[0]=err_c2j_0val[0]
            Ch_sig0[1]=0.0
            el0plot=[el_0[valid0],el_0[valid0]]
            c2J0plot=[c2j_0val,c2j_0val]
            th0plot=[th_0val,th_0val]
            errc2j0plot=[err_c2j_0val,err_c2j_0val]
          endelse

          t=th_0val
          y_0=Ch_0[0]+Ch_0[1]*t
          y_0min=Ch_0[0]+(Ch_0[1])*th_0val+Ch_sig0[0]
          y_0max=Ch_0[0]+(Ch_0[1])*th_0val-Ch_sig0[0]

          ; printouts for test usage only
          ;          print, '***********************'
          ;          print, 'Interval ',j+1
          ;          print, 'thisgap ', thisgap
          ;          print, 'ti, tf ', ti, tf
          ;          print, 'c2J_0 ',c2j_0
          ;          print, 'c2J_0val ',c2j_0val
          ;          print, 'valid0 ', valid0
          ;          print, 'th_0val ',th_0val
          ;          print, '***********************'
          ;          print, ' '
          ;          stop

          realc2j0=where(c2J0plot gt 0)
          ystep=(max(c2J0plot[realc2j0]+errc2J0plot[realc2j0])-min(c2J0plot[realc2j0]-errc2J0plot[realc2j0]))/10.0
          labstep=ystep
          if ystep eq 0 then begin
            ystep=0.1
            labstep=ystep/numcals
          endif

          istherec2J0=n_elements(c2J0plot)
          if istherec2J0 ne 0 then begin
            ploterror, el0plot[realc2j0], c2J0plot[realc2j0], el_0val[realc2j0]*0, errc2j0plot[realc2j0], /NODATA, xrange=[0,90], ys=1, ytitle="cnt2Jy_0", xtitle="Elevation (deg)", charsize=1.0
            ;plot, el0plot, c2J0plot, /NODATA, xrange=[0,90], yrange=[min(c2J0plot-ystep),max(c2J0plot+ystep)], ytitle="cnt2Jy_0", xtitle="Elevation (deg)"
          endif
          for i=0,n_elements(nameu)-1 do begin
            thiscal=where(name[valid0] eq nameu[i], namenum)
            if namenum ge 1 then begin
              if istherec2J0 ne 0 then begin
                oplot, el0plot[thiscal], c2J0plot[thiscal], color=i,  psym=2
                oploterror, el0plot[thiscal], c2J0plot[thiscal], el_0val[thiscal]*0, errc2j0plot[thiscal], color=i,  psym=2
                xyouts, 0, max(c2J0plot)-i*labstep, nameu[i], color=i, charsize=0.8, charthick=2
              endif
            endif
          endfor


          xstep=(max(th0plot)-min(th0plot))/10.0 
          if xstep eq 0 then xstep=1

          if keyword_set(detplot) then begin
            if istherec2J0 ne 0 then begin
              displot=errorplot(th0plot[realc2j0], c2J0plot[realc2j0], th0plot[realc2j0]*0, errc2j0plot[realc2j0], /NODATA, layout=[1,2,1])
              displot.title=origin+' cnt2Jy_0 vs Time, INT'+strcompress(string(j),/remove_all)
              displot.xtitle='Elapsed time (hh.h)'
              displot.ytitle='Jy/cnt'
           ;   displot.xrange=[min(th0plot)-xstep,max(th0plot)+xstep]
           ;   displot.yrange=[min(c2J0plot[realc2j0])-ystep,max(th0plot[realc2j0])+ystep]
            endif
          endif

          if istherec2J0 ne 0 then begin
            ploterror, th0plot[realc2j0], c2J0plot[realc2j0], th0plot[realc2j0]*0, errc2j0plot[realc2j0], /NODATA, xrange=[min(th0plot)-xstep,max(th0plot)+xstep], ys=1, ytitle="cnt2Jy_0", xtitle="Elapsed hours since "+day, XTICKFORMAT='(f5.2)', charsize=1.0
            ; ploterror, th0plot, c2J0plot, th0plot*0, errc2j0plot, /NODATA, xs=1, yrange=[min(c2J0plot-ystep),max(c2J0plot+ystep)], ytitle="cnt2Jy_0", xtitle="Elapsed hours since "+day, XTICKFORMAT='(f5.2)'
            ;
          endif
          for i=0,n_elements(nameu)-1 do begin
            thiscal=where(name[valid0] eq nameu[i], namenum)
            if namenum ge 1 then begin

              if istherec2J0 ne 0 then begin
                oplot, th0plot[thiscal], c2J0plot[thiscal], color=i, psym=2
                oploterror, th0plot[thiscal], c2J0plot[thiscal], th0plot[thiscal]*0, errc2j0plot[thiscal], color=i, psym=2
              endif

              if keyword_set(detplot) then begin
                if istherec2J0 ne 0 then begin
                  ; object graphics
                  displot2=errorplot(th0plot[thiscal], c2J0plot[thiscal], th0plot[thiscal]*0, errc2j0plot[thiscal], NAME=nameu[i], /overplot)
                  displot2.symbol=sym[i]
                  displot2.linestyle=' '
                  displot2.color=colors[i]
                  displot2.sym_size=1.5
                  displot2.ERRORBAR_COLOR=colors[i]
                  displot2.ERRORBAR_CAPSIZE=0.2
                  displot2.sym_filled=1
                  displot3=plot(th_0val, y_0,/overplot)
                  displot3.color='dark slate gray'
                  displot3.linestyle='solid'
                  displot4=plot(th_0val, y_0min,/overplot)
                  displot4.color='dark slate gray'
                  displot4.linestyle='dashed'
                  displot5=plot(th_0val, y_0max,/overplot)
                  displot5.color='dark slate gray'
                  displot5.linestyle='dashed'
                  lab=text(max(th0plot)-xstep,max(c2J0plot+errc2j0plot)-(i+1)*ystep, nameu[i], color=colors[i],/data, font_size=10)
                endif
              endif

            endif
          endfor


          oplot, th_0val, y_0, linestyle=0
          oplot, th_0val, y_0min, linestyle=1
          oplot, th_0val, y_0max, linestyle=1

        endif else begin
          ; there are no usable measurements for CH0!
          print, infiles[f]+' Time range '+string(j, format='(I2)')+' - Invalid measurements for CH0'
        endelse


        if valid1[0] ne -1 then begin
          ; there are usable measurements
          t_1val=t_1[valid1]
          el_1val=el_1[valid1]
          c2j_1val=c2j_1[valid1]
          err_c2j_1val=err_c2j_1[valid1]
          th_1val=elapsed[valid1]

          if valnum1 gt 1 then begin
            ; there are at least two usable measurements
            w_1=1/(err_c2j_1val/c2j_1val)^2
            meanerr,c2j_1val,err_c2j_1val,x_1mean,sigma_1m,sigma_1d
            C_1 = LINFIT(t_1val, c2J_1val, CHISQ=chisq, /DOUBLE, MEASURE_ERRORS=err_c2J_1val, PROB=lik_1, SIGMA=C_sig1 )
            Ch_1 = LINFIT(th_1val, c2J_1val, CHISQ=chisq_h, /DOUBLE, MEASURE_ERRORS=err_c2J_1val_h, PROB=lik_1_h, SIGMA=Ch_sig1)
            el1plot=el_1val
            c2J1plot=c2j_1val
            errc2j1plot=err_c2j_1val
            th1plot=th_1val
          endif else begin
            ;only one measurement is available: linear fitting cannot be performed!
            x_1mean=c2j_1val[0]
            sigma_1d=err_c2j_1val[0]
            sigma_1m=err_c2j_1val[0]
            C_1[0]=c2j_1val[0]
            C_1[1]=0.0
            C_sig1[0]=err_c2j_1val[0]
            C_sig1[1]=0.0
            Ch_1[0]=c2j_1val[0]
            Ch_1[1]=0.0
            Ch_sig1[0]=err_c2j_1val[0]
            Ch_sig1[1]=0.0
            el1plot=[el_1[valid1],el_1[valid1]]
            c2J1plot=[c2j_1val,c2j_1val]
            th1plot=[th_1val,th_1val]
            errc2j1plot=[err_c2j_1val,err_c2j_1val]
          endelse

          t=th_1val
          y_1=Ch_1[0]+Ch_1[1]*th_1val
          y_1min=Ch_1[0]+(Ch_1[1])*th_1val+Ch_sig1[0]
          y_1max=Ch_1[0]+(Ch_1[1])*th_1val-Ch_sig1[0]

          realc2j1=where(c2J1plot gt 0)
          ystep=(max(c2J1plot[realc2j1]+errc2J1plot[realc2j1])-min(c2J1plot[realc2j1]-errc2J1plot[realc2j1]))/10.0
          labstep=ystep
          if ystep eq 0 then begin
            ystep=0.1
            labstep=ystep/numcals
          endif

          istherec2J1=n_elements(c2J1plot)
          if istherec2J1 ne 0 then begin
            ploterror, el1plot[realc2j1], c2J1plot[realc2j1], el1plot[realc2j1]*0, errc2J1plot[realc2j1], /NODATA, xrange=[0,90], ys=1, ytitle="cnt2Jy_1", xtitle="Elevation (deg)", charsize=1.0
            ;   plot, el1plot, c2J1plot, /NODATA, xrange=[0,90], yrange=[min(c2J1plot-ystep),max(c2J1plot+ystep)], ytitle="cnt2Jy_1", xtitle="Elevation (deg)"
          endif
          for i=0,n_elements(nameu)-1 do begin
            thiscal=where(name[valid1] eq nameu[i], namenum)
            if namenum ge 1 then begin
              if istherec2J1 ne 0 then begin
                oplot, el1plot[thiscal], c2J1plot[thiscal], color=i,  psym=6
                oploterror, el1plot[thiscal], c2J1plot[thiscal], el1plot[thiscal]*0, errc2J1plot[thiscal], color=i,  psym=6
                xyouts, 0, max(c2J1plot)-i*labstep, nameu[i], color=i, charsize=0.8, charthick=2
              endif
            endif
          endfor


          xstep=(max(th1plot)-min(th1plot))/10.0
          if xstep eq 0 then xstep=1

          if keyword_set(detplot) then begin
            if istherec2J1 ne 0 then begin
              displot6=errorplot(th1plot[realc2j1], c2J1plot[realc2j1], th1plot[realc2j1]*0, errc2j1plot[realc2j1], /NODATA, layout=[1,2,2], /current)
              displot6.title=origin+' cnt2Jy_1 vs Time, INT'+strcompress(string(j),/remove_all)
              displot6.xtitle='Elapsed time (hh.h)'
              displot6.ytitle='Jy/cnt'
            ;  displot6.xrange=[min(th1plot)-xstep,max(th1plot)+xstep]
            ;  displot6.yrange=[min(c2J1plot[realc2j1])-ystep,max(th1plot[realc2j1])+ystep]
            endif
          endif

          if istherec2J1 ne 0 then begin
            ploterror, th1plot[realc2j1], c2J1plot[realc2j1], th1plot[realc2j1]*0, errc2J1plot[realc2j1], /NODATA, xrange=[min(th1plot)-xstep,max(th1plot)+xstep], ys=1, ytitle="cnt2Jy_1", xtitle="Elapsed hours since "+day, XTICKFORMAT='(f5.2)', charsize=1.0
            ; ploterror, th1plot, c2J1plot, th1plot*0, errc2J1plot, /NODATA, xs=1, yrange=[min(c2J1plot-ystep),max(c2J1plot+ystep)], ytitle="cnt2Jy_1", xtitle="Elapsed hours since "+day, XTICKFORMAT='(f5.2)'
          endif
          for i=0,n_elements(nameu)-1 do begin
            thiscal=where(name[valid1] eq nameu[i], namenum)
            if namenum ge 1 then begin
              if istherec2J1 ne 0 then begin
                oplot, th1plot[thiscal], c2J1plot[thiscal], color=i,  psym=6
                oploterror, th1plot[thiscal], c2J1plot[thiscal], th1plot[thiscal]*0, errc2J1plot[thiscal], color=i, psym=6
              endif

              if keyword_set(detplot) then begin
                if istherec2J1 ne 0 then begin
                  ; object graphics
                  displot7=errorplot(th1plot[thiscal], c2J1plot[thiscal], th1plot[thiscal]*0, errc2j1plot[thiscal], NAME=nameu[i], /overplot)
                  displot7.symbol=sym[i]
                  displot7.linestyle=' '
                  displot7.color=colors[i]
                  displot7.sym_size=1.5
                  displot7.ERRORBAR_COLOR=colors[i]
                  displot7.ERRORBAR_CAPSIZE=0.2
                  displot7.sym_filled=1
                  displot8=plot(th_1val, y_1,/overplot)
                  displot8.color='dark slate gray'
                  displot8.linestyle='solid'
                  displot9=plot(th_1val, y_1min,/overplot)
                  displot9.color='dark slate gray'
                  displot9.linestyle='dashed'
                  displot10=plot(th_1val, y_1max,/overplot)
                  displot10.color='dark slate gray'
                  displot10.linestyle='dashed'
                endif
              endif

            endif

          endfor
          oplot, th_1val, y_1, linestyle=0
          oplot, th_1val, y_1min, linestyle=1
          oplot, th_1val, y_1max, linestyle=1

        endif else begin
          ; there are no usable measurements for CH1!
          print, infiles[f]+' Time range '+string(j, format='(I2)')+' - Invalid measurements for CH1'
        endelse


        case mode of
          'aver': begin
            printf, Unit, FORMAT = '(i4," ",d8.1,2(d11.4),8(g13.7))',j+1,freq[0],ti,tf,0,x_0mean,0,sigma_0m,0,x_1mean,0,sigma_1m
            print, x_0mean, sigma_0m, sigma_0d
            print, x_1mean, sigma_1m, sigma_1d
          end
          'linf': begin
            printf, Unit, FORMAT = '(i4," ",d8.1,2(d11.4),8(g13.7))',j+1,freq[0],ti,tf,C_0[1],C_0[0],C_sig0[1],C_sig0[0],C_1[1],C_1[0],C_sig1[1],C_sig1[0]
          end
        endcase

      endif else begin
        print, 'No measurements in time interval number ', strcompress(string(j+1),/remove_all)
        ; printf, Unit, 'Interval number ',strcompress(string(j+1),/remove_all),' does not contain measurements'
      endelse

    endfor

  endfor

  device, /close
  close, /ALL

  print, ' ****************************** '
  print, ' ****  C2JTIMELINE is DONE **** '
  print, ' *** Next step is RUNTARGET *** '
  print, ' ****************************** '

  return
end



;pro fit,ascissa,somma,ch,


pro meanerr,x,sigmax,xmean,sigmam,sigmad
  ;+
  ; NAME:
  ;meanerr
  ; PURPOSE: (one line)
  ;Calculate the mean and estimated errors for a set of data points
  ; DESCRIPTION:
  ;This routine is adapted from Program 5-1, XFIT, from "Data Reduction
  ;and Error Analysis for the Physical Sciences", p. 76, by Philip R.
  ;Bevington, McGraw Hill.  This routine computes the weighted mean using
  ;Instrumental weights (w=1/sigma^2).
  ; CATEGORY:
  ;Statistics
  ; CALLING SEQUENCE:
  ;meanerr,x,sigmax,xmean,sigmam,sigmad
  ; INPUTS:
  ;x      - Array of data points
  ;sigmax - array of standard deviations for data points
  ; OPTIONAL INPUT PARAMETERS:
  ;None.
  ; KEYWORD PARAMETERS:
  ;None.
  ; OUTPUTS:
  ;xmean  - weighted mean
  ;sigmam - standard deviation of mean
  ;sigmad - standard deviation of data
  ; COMMON BLOCKS:
  ;None.
  ; SIDE EFFECTS:
  ;None.
  ; RESTRICTIONS:
  ;None.
  ; PROCEDURE:
  ; MODIFICATION HISTORY:
  ;Written by Marc W. Buie, Lowell Observatory, 1992 Feb 20
  ;-

  if n_elements(x) eq 1 then begin
    xmean  = x[0]
    sigmam = sigmax[0]
    sigmad = sigmax[0]
  endif else begin
    weight = 1.0/sigmax^2
    sum    = total(weight)
    if sum eq 0.0 then print,'MEANERR: sum is zero.'
    sumx   = total(weight*x)
    xmean  = sumx/sum
    sigmam = sqrt(1.0/sum)
    sigmad = sqrt(total((x-xmean)^2)/(n_elements(x)-1))
  endelse

end