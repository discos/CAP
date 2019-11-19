pro findh0key, header, stringa, key, value, info, rflag
  
  ; Procedure to read keyword values from FITS headers
  ;
  ; Author: Simona Righini
  ; Last edited: Oct 6, 2017 

  
  rflag=1
  search=strmatch(header, stringa)
  index=where(search eq 1)
  sp=path_sep()
  if index ne -1 and stringa eq '*Converted FITS*' then return 
  if index ne -1 and stringa eq '*Obtained from ESCS 0.1 TP FITS*' then return
  if index ne -1 then begin
    content=strsplit(header[index],'=',/extract)
    key=content[0]
    value_and_info=strsplit(string(content[1]),'/',/extract)
    if stringa ne '*ScheduleName =*' then begin
      value=value_and_info[0]
      info=value_and_info[1]
    endif else begin
      steps=n_elements(value_and_info)
      info=value_and_info[-1]
      value=value_and_info[0]
      for s=1, steps-2 do begin
        value=value+sp+value_and_info[s]
      endfor
    endelse
  endif else begin
    rflag=0
  endelse
end