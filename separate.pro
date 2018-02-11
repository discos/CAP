
; This procedure separates the acquisitions into sub-folders, depending on the frequency band and the
; calibrator/skydip/target nature of the scan.
;
; The usable RF bands are hard-coded.
; The original version of this procedure considers the following bands and frequency boundaries (MHz):
;
; pbounds=[300,1000]
; lbounds=[1000,2000]
; sbounds=[2000,4000]
; cbounds=[4000,8000]
; xbounds=[8000,10000]
; kbounds=[18000,26500]
;
; If observing with more than one sub-band in the K-band RF range, specify /multiK as in
;
; IDL> separa, /multiK    (in this case, kbounds=[18000,20000,23000,26500])
;
; Default flux density calibrators are considered to be 3C48, 3C123, 3C286, 3C295, NGC7027.
; When wanting to specify a custom list of calibrators, provide an external text file
; containing as columns the source names and flux densities (the same to be used when running "runcalib"),
; and then use
;
; IDL> separa, /extlist
;
; Latest update: Feb 11th, 2018


pro separate, multiK=multiK, extlist=extlist

  print, ' '
  sp=path_sep()
  pin=dialog_pickfile(/DIRECTORY, title='Select main folder, containing the YYYYMMDD data folders')
  print, 'Selected folder is ', pin

  if not keyword_set(multiK) then Kband='single' else Kband='multi'

  homecals=['3C48','3C123','3C286','3C295','NGC7027']   ; default flux density calibrators

  if keyword_set(extlist) then begin
    myextlist=dialog_pickfile(TITLE='Please pick EXTERNAL LIST of calibrators')
    readcol, myextlist, extcalname, extcalflux, format='(A,D)', /SILENT
    calname=strupcase(extcalname)
  endif else begin
    myextlist=' '
    calname=strupcase(homecals)
  endelse

  ; In this module the frequency bands are defined
  ; Users might want to add/edit some
  ; ***********************************************

  pbounds=[300,1000]
  lbounds=[1000,2000]
  sbounds=[2000,4000]
  cbounds=[4000,8000]
  xbounds=[8000,10000]

  ; NOTICE: in case of editing, edit the boundaries implementation as well! Search for string "edit here" in this code.

  if Kband eq 'single' then begin
    bands=['P','L','S','C','X','K']
    kbounds=[18000,26500]
  endif else begin
    bands=['P','L','S','C','X','K-low','K-mid','K-hi']
    kbounds=[18000,20000,23000,26500] ; definition of sub-bands boundaries for K-band
  endelse

  ; ***********************************************

  bandnum=n_elements(bands)
  mybounds=make_array(1,bandnum,/integer)


  ; Checking whether the folder for DAT files (produced by the CalibrationTool) exists.
  ; If not, it is created.
  ans0=file_test(pin+'DATS',/DIRECTORY)
  if ans0 eq 0 then begin
    print, 'Creating subfolder for .dat files...'
    file_mkdir, pin+'DATS'
  endif

  print, 'Creating subfolders for band-separated files...'
  for b=0, bandnum-1 do begin

    ; Checking whether the band-dependant and type-dependant folders alreasy exist.
    ; If not, they are created.

    ans1=file_test(pin+bands[b],/DIRECTORY)

    if ans1 eq 0 then begin
      file_mkdir, pin+bands[b]
    endif

    ans2=file_test(pin+bands[b]+sp+'CALIBRATORS',/DIRECTORY)
    ans3=file_test(pin+bands[b]+sp+'TARGETS',/DIRECTORY)
    ans4=file_test(pin+bands[b]+sp+'SKYDIPS',/DIRECTORY)

    if ans2 eq 0 then begin
      ; print, 'Creating folder for '+bands[b]+'-band CALIBRATORS files'
      file_mkdir, pin+bands[b]+sp+'CALIBRATORS'
    endif
    if ans3 eq 0 then begin
      ; print, 'Creating folder for '+bands[b]+'-band TARGETS files'
      file_mkdir, pin+bands[b]+sp+'TARGETS'
    endif
    if ans4 eq 0 then begin
      ; print, 'Creating folder for '+bands[b]+'-band SKYDIP files'
      file_mkdir, pin+bands[b]+sp+'SKYDIPS'
    endif

  endfor


  list = FILE_SEARCH(pin+'20*', COUNT=number, /TEST_DIRECTORY)

  for i=0,number-1 do begin  ; loop over 'number' days
    ; print, list[i]
    datlist = file_search(list[i]+sp+'*.dat', COUNT=datnum)
    if datnum ne 0 then begin
      for j=0,datnum-1 do begin
        file_move, datlist[j], pin+'DATS'+sp, /OVERWRITE
      endfor
    endif

    sublist = FILE_SEARCH(list[i]+sp+'20*', /TEST_DIRECTORY, COUNT=subnum)
    for j=0,subnum-1 do begin  ; loop over 'subnum' scans
      subsublist = FILE_SEARCH(sublist[j]+sp+'*.fits', COUNT=subsubnum)
      ; obtaining info from the first subscan (i.e. the first FITS in the folder)
      RFinfo=mrdfits(subsublist[0],2,/silent)
      sky_bandwidth=RFinfo[0].bandWidth
      sky_freq=RFinfo[0].frequency+sky_bandwidth/2.0
      mainh=mrdfits(subsublist[0],0,head0,/silent)
      findh0key, head0, '*SOURCE*', keylab, target, infolab, tarflag

      ; elaborating the target name string
      target=strcompress(target, /remove_all)
      cleantarget=strsplit(target,"'",/extract)
      target=cleantarget[0]

      ; directing the file to the proper location
      type='TARGETS'+sp   ; default choice

      checkcal=where(target eq calname, isitcal)
      if isitcal ne 0 then begin
        type = 'CALIBRATORS'+sp ; overriding the default choice: the source is a calibrator
        if subsubnum eq 1 then begin ; to handle cases when SKYDIPS have been mis-named
          data=mrdfits(subsublist[0],4,/SILENT)
          span=abs(data[0].el-data[-1].el)
          if span gt 45/180*!dpi then type = 'SKYDIPS'+sp
        endif
      endif

      if (strmatch(sublist[j],'*sky*',/FOLD_CASE)) then type = 'SKYDIPS'+sp ; overriding the default choice: the file is a SKYDIP

      ; actually moving the folders
      if (i eq 0) and (j eq 0) then print, 'Moving files...'
      for b=0, bandnum-1 do begin
        case bands[b] of
          ; edit here in case you need to add more bands
          'P': mybounds=pbounds
          'L': mybounds=lbounds
          'S': mybounds=sbounds
          'C': mybounds=cbounds
          'X': mybounds=xbounds
          'K': mybounds=kbounds
          'K-low': mybounds=[kbounds[0],kbounds[1]]
          'K-mid': mybounds=[kbounds[1],kbounds[2]]
          'K-hi': mybounds=[kbounds[2],kbounds[3]]
        endcase
        intervals=n_elements(mybounds)-1
        for int=0, intervals-1 do begin
          if (sky_freq ge mybounds[int] and sky_freq lt mybounds[int+1]) then begin
            outband = pin+bands[b]+sp
            file_move, sublist[j], outband+type, /OVERWRITE
            break
          endif
        endfor
      endfor

    endfor ; ending cycle on scans

  endfor; ending cycle on days

  ; removing empty/unused folders (the non-empty ones will be by default spared!)
  print, 'Removing empty subfolders...'
  for b=0, bandnum-1 do begin
    file_delete, pin+bands[b]+sp+'CALIBRATORS',/quiet
    file_delete, pin+bands[b]+sp+'TARGETS',/quiet
    file_delete, pin+bands[b]+sp+'SKYDIPS',/quiet
    file_delete, pin+bands[b]+sp,/quiet
  endfor
  file_delete, pin+'DATS',/quiet

  tobedel=file_search(pin+'20*', /TEST_DIRECTORY, count=delnum)
  for dn=0, delnum-1 do begin
    file_delete, tobedel[dn],/quiet
  endfor

  print, ' '
  print, '**********'
  print, '   DONE    '
  print, '**********'
  print, ' '
  
  return

end


