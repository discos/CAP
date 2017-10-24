pro dataplot, scanflag, X, Y, allnamefile, mypath, section, datafit, cnt2Jy, Level, SNR, residual, direflag, tipofit, pdis, stac
            
  ; This procedure is in charge of the plotting operations
  ;
  ; Authors: Marcello Giroletti, Simona Righini
  ; Last edited: Oct 5, 2017
  ;

  sep=path_sep()
  if strmatch(allnamefile, '*'+sep+'*') then begin
    splitted=strsplit(allnamefile,sep,/extract)
    namefile=splitted[-1]
  endif else begin
    namefile=allnamefile
  endelse
  plot_title = namefile+" "+section+" "+direflag+' '+tipofit

  if (scanflag eq 0) then begin
    if (section eq 'Ch_0') then begin
      pdis=plot(X, Y, buffer=1, layout=[2,2,1], dimensions=[1024,768], title = plot_title, /NODATA)
      pdis=plot(X, Y, buffer=1, layout=[2,2,2], dimensions=[1024,768], title = "RESIDUALS", /NODATA, /CURRENT)
    endif else begin
      pdis=plot(X, Y, buffer=1, layout=[2,2,3], dimensions=[1024,768], title = plot_title, /NODATA, /CURRENT)
      pdis=plot(X, Y, buffer=1, layout=[2,2,4], dimensions=[1024,768], title = "RESIDUALS", /NODATA, /CURRENT)
    endelse
  endif else begin
    if (section eq 'Ch_0') then begin
      nw=1
      ; pdis=plot(X, Y, layout=[2,2,nw], dimensions=[1024,768],$
      pdis=plot(X, Y, buffer=1, layout=[2,2,nw], dimensions=[1024,768],$
        title = plot_title, xtitle = "x_J2000 [deg]", xrange = [X[0],X[-1]], ytitle = "signal [cnt]")
    endif else begin
      nw=3
      ; pdis=plot(X, Y, layout=[2,2,nw], /current, $
      pdis=plot(X, Y, buffer=1, layout=[2,2,nw], /current, $
        title = plot_title, xtitle = "x_J2000 [deg]", xrange = [X[0],X[-1]], ytitle = "signal [cnt]")
    endelse
    pdis.symbol=2
    pdis.sym_size=0.2

    ; overplot of the gaussian fit
    pdis=plot(X, datafit, /overplot)
    ; plot residuals in side panel
    ; pdis=plot(X, residual, layout=[2,2,nw+1], /current, title="RESIDUALS")
    pdis=plot(X, residual,buffer=1, layout=[2,2,nw+1], /current, title="RESIDUALS")
    pdis.xtitle = "x_J2000 [deg]"
    pdis.ytitle = "signal [cnt]"
  endelse
  if stac eq 'single' then begin
    if tipofit eq 'linear' then pdis.save, mypath+sep+namefile+'_lin_single.pdf', /LANDSCAPE, RESOLUTION=300
    if tipofit eq 'cubic' then pdis.save, mypath+sep+namefile+'_cub_single.pdf', /LANDSCAPE, RESOLUTION=300
  endif
  return
end
