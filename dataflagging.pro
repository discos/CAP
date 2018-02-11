; Program to visual inspect and flag continuum FITS acquisitions (ESCS 0.2-0.3)
; Created on Feb. 6th 2013 by A.Mattana and S.Righini
; Debug history:
;      7 feb 2013 - added "app" variable in "start", to create/append checkfile
;     12 feb 2013 - fixed problem on Logical Units, by adding free_lun
;     28 feb 2013 - insertion of "BACK" and "SKIP" buttons
; New releases:
;     27 apr 2016 - the working folder is now the "parent"one, containing all the subfolders 
;                   where the FITS files are stored. Users can now flag data easily passing 
;                   from one subfolder to the following/precedent one

; USAGE:
; 1) Compile and run dataflagging.pro; 
; 2) (optional) select the user name's initials from the upper drop-down menu;
; 3) click on BROWSE button and select the parent folder containing the FITS-full folders; 
; 4) click on GO!;
; 5) if a checkfile (i.e. the list of flagged files) already exists for a FITS folder, 
;    a dialog box will appear. By answering "YES", users indicate they want to keep those 
;    flags and go on appending new ones for the so-far-unflagged data (if any). 
;    By answering "NO", they decide to start over, thus rewriting the checkfile; 
; 6) A new window will appear, where the content of each FITS file will be displayed. 
;    The two plots show: 
;      left --> Section CH0
;      right --> Section CH1  
;    Users are asked to express their opinion on the sections to be KEPT (i.e. sent
;    to the following phases of the data reduction) by clicking on NONE, LEFT, RIGHT, BOTH. 
;    The button "BOTH (V. good)" can be ignored, or used to pinpoint the particularly 
;    good acquisitions (for future uses, such as the implementation of AI tools);
; 7) the "Skip" button skips the FITS presently displayed, leaving it unflagged (beware: CAP will by defaut consider it good!);
; 8) using the "Back" button, the previous FITS is re-displayed and re-flagged (the old flag is overwritten).
;    It is not possible to invoke "Back" after a "Skip" (for this reason the button is disabled);
; 9) to interrupt the procedure, yet saving all the previously-performed flagging, click on "STOP and EXIT".


pro WID_BASE_0, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_

  common wid, WID_BASE_0, WID_BOX, wid_left_id, wid_right_id $
    ,WID_DRAW_LEFT, Button_BROWSE $
    ,TEXT_OBS_PATH, Button_START, WID_TAB $
    ,Button_LEFT,Button_RIGHT,Button_PERFECT,Button_ALL,Button_NONE,WID_BUTTON_STOP $
    ,WID_BUTTON_SKIP, WID_BUTTON_BACK, WID_BUTTON_PREVDIR, WID_BUTTON_NEXTDIR

  common cfg, folder, dirlist, dirindex, userlist, user, outname, sp, list, $
    flagged, tobeflagged, flagnext, flags, app, index

  sp=path_sep()

  Resolve_Routine, 'WID_BASE_0_event';,/COMPILE_FULL_FILE  ; Load event callback routines

  WID_BASE_0 = Widget_Base( GROUP_LEADER=wGroup, UNAME='WID_BASE_0'  $
    ,XOFFSET=5 ,YOFFSET=5 ,SCR_XSIZE=1200 ,SCR_YSIZE=750  $
    ,TITLE="FLAGGING: WHAT'S TO KEEP?" ,SPACE=3 ,XPAD=3 ,YPAD=3)

  WID_BOX = Widget_Base( GROUP_LEADER=boxes, UNAME='WID_BOX'  $
    ,XOFFSET=5 ,YOFFSET=5 ,SCR_XSIZE=400 ,SCR_YSIZE=300  $
    ,TITLE="ENTER CHECKFILE CUSTOM NAME" ,SPACE=3 ,XPAD=3 ,YPAD=3)

  WID_TAB = Widget_Tab(WID_BASE_0, UNAME='WID_TAB'  $
    ,XOFFSET=5 ,YOFFSET=5 ,SCR_XSIZE=1200 ,SCR_YSIZE=750  $
    )

  WID_TAB_CFG = Widget_Base(WID_TAB, UNAME='WID_TAB_CFG'  $
    ,XOFFSET=5 ,YOFFSET=5  $
    ,TITLE="    CONFIG   ")

  WID_TAB_PLOT = Widget_Base(WID_TAB, UNAME='WID_TAB_PLOT'  $
    ,XOFFSET=5 ,YOFFSET=5  $
    ,TITLE="     PLOT    ")

  WID_DRAW_LEFT = Widget_Draw(WID_TAB_PLOT, UNAME='WID_DRAW_LEFT'  $
    ,XOFFSET=83 ,YOFFSET=25 ,SCR_XSIZE=450 ,SCR_YSIZE=350)

  WID_DRAW_RIGHT = Widget_Draw(WID_TAB_PLOT, UNAME='WID_DRAW_RIGHT'  $
    ,XOFFSET=566 ,YOFFSET=25 ,SCR_XSIZE=450 ,SCR_YSIZE=350)


  ;Materiale nel TAB CFG
  userlist=['NA','U1','U2','U3','U4','U5','U6','U7','U8','U9']
  flags = ['-1c0','-1c1','+0cx','+1cx','-1cx']
  user='NA'

  WID_USER_LABEL = Widget_Label(WID_TAB_CFG, UNAME='WID_USER_LABEL'  $
    ,XOFFSET=300 ,YOFFSET=50 ,SCR_XSIZE=250 ,SCR_YSIZE=50  $
    ,/ALIGN_RIGHT ,VALUE='Select User ID (or leave default value)')

  DROPLIST_USERS = Widget_Droplist(WID_TAB_CFG,  $
    UNAME='DROPLIST_USERS' ,XOFFSET=300 ,YOFFSET=100  $
    ,SCR_XSIZE=250 ,SCR_YSIZE=40, value=userlist)

  WID_FOLDER_LABEL = Widget_Label(WID_TAB_CFG, UNAME='WID_FOLDER_LABEL'  $
    ,XOFFSET=300 ,YOFFSET=200 ,SCR_XSIZE=400 ,SCR_YSIZE=40  $
    ,/ALIGN_LEFT ,VALUE='Select main folder (containing the subfolders with FITS files)')

  Button_BROWSE = Widget_Button(WID_TAB_CFG, UNAME='Button_BROWSE'  $
    ,XOFFSET=200 ,YOFFSET=240 ,SCR_XSIZE=100 ,SCR_YSIZE=50  $
    ,/ALIGN_CENTER ,VALUE='BROWSE')

  TEXT_OBS_PATH = Widget_Text(WID_TAB_CFG, UNAME='TEXT_OBS_PATH'  $
    ,XOFFSET=300 ,YOFFSET=240 ,SCR_XSIZE=600 ,SCR_YSIZE=50 $
    ,XSIZE=40,YSIZE=1)

  Button_START = Widget_Button(WID_TAB_CFG, UNAME='Button_START'  $
    ,XOFFSET=450 ,YOFFSET=400, SCR_XSIZE=150 ,SCR_YSIZE=150  $
    ,/ALIGN_CENTER ,VALUE='START!')

  WID_HEADER_LABEL = Widget_Label(WID_TAB_PLOT, UNAME='WID_HEADER_LABEL'  $
    ,XOFFSET=375 ,YOFFSET=390 ,SCR_XSIZE=350 ,SCR_YSIZE=40  $
    ,/ALIGN_CENTER ,VALUE='Indicate the section(s) to be kept for analysis')

  Button_LEFT = Widget_Button(WID_TAB_PLOT, UNAME='Button_LEFT'  $
    ,XOFFSET=200 ,YOFFSET=450 ,SCR_XSIZE=100 ,SCR_YSIZE=100  $
    ,/ALIGN_CENTER ,VALUE='LEFT')

  Button_RIGHT = Widget_Button(WID_TAB_PLOT, UNAME='Button_RIGHT'  $
    ,XOFFSET=350 ,YOFFSET=450 ,SCR_XSIZE=100 ,SCR_YSIZE=100  $
    ,/ALIGN_CENTER ,VALUE='RIGHT')

  Button_PERFECT = Widget_Button(WID_TAB_PLOT, UNAME='Button_PERFECT'  $
    ,XOFFSET=500 ,YOFFSET=450 ,SCR_XSIZE=100 ,SCR_YSIZE=100  $
    ,/ALIGN_CENTER ,VALUE='BOTH (V.GOOD)')

  Button_ALL = Widget_Button(WID_TAB_PLOT, UNAME='Button_ALL'  $
    ,XOFFSET=650 ,YOFFSET=450 ,SCR_XSIZE=100 ,SCR_YSIZE=100  $
    ,/ALIGN_CENTER ,VALUE='BOTH')

  Button_NONE = Widget_Button(WID_TAB_PLOT, UNAME='Button_NONE'  $
    ,XOFFSET=800 ,YOFFSET=450 ,SCR_XSIZE=100 ,SCR_YSIZE=100  $
    ,/ALIGN_CENTER ,VALUE='NONE')

  WID_BUTTON_PREVDIR = Widget_Button(WID_TAB_PLOT,  $
    UNAME='WID_BUTTON_PREVDIR' ,XOFFSET=180 ,YOFFSET=600  $
    ,SCR_XSIZE=120 ,SCR_YSIZE=50 ,/ALIGN_CENTER ,VALUE='<< PREV DIR')

  WID_BUTTON_BACK = Widget_Button(WID_TAB_PLOT,  $
    UNAME='WID_BUTTON_BACK' ,XOFFSET=350 ,YOFFSET=600  $
    ,SCR_XSIZE=100 ,SCR_YSIZE=50 ,/ALIGN_CENTER ,VALUE='< BACK')

  WID_BUTTON_STOP = Widget_Button(WID_TAB_PLOT,  $
    UNAME='WID_BUTTON_STOP' ,XOFFSET=500 ,YOFFSET=600  $
    ,SCR_XSIZE=100 ,SCR_YSIZE=50 ,/ALIGN_CENTER ,VALUE='STOP and EXIT')

  WID_BUTTON_SKIP = Widget_Button(WID_TAB_PLOT,  $
    UNAME='WID_BUTTON_SKIP' ,XOFFSET=650 ,YOFFSET=600  $
    ,SCR_XSIZE=100 ,SCR_YSIZE=50 ,/ALIGN_CENTER ,VALUE='SKIP >')

  WID_BUTTON_NEXTDIR = Widget_Button(WID_TAB_PLOT,  $
    UNAME='WID_BUTTON_NEXTDIR' ,XOFFSET=800 ,YOFFSET=600  $
    ,SCR_XSIZE=120 ,SCR_YSIZE=50 ,/ALIGN_CENTER ,VALUE='NEXT DIR >>')


  Widget_Control, /REALIZE, WID_BASE_0
  WIDGET_CONTROL, WID_DRAW_RIGHT, GET_VALUE=wid_right_id
  WIDGET_CONTROL, WID_DRAW_LEFT, GET_VALUE=wid_left_id
  widget_control, WID_BUTTON_PREVDIR, SENSITIVE=0
  widget_control, WID_BUTTON_BACK, SENSITIVE=0
  widget_control, WID_BUTTON_STOP, SENSITIVE=0
  widget_control, WID_BUTTON_SKIP, SENSITIVE=0
  widget_control, WID_BUTTON_NEXTDIR, SENSITIVE=0
  widget_control, Button_LEFT, SENSITIVE=0
  widget_control, Button_RIGHT, SENSITIVE=0
  widget_control, Button_PERFECT, SENSITIVE=0
  widget_control, Button_ALL, SENSITIVE=0
  widget_control, Button_NONE, SENSITIVE=0



  XManager, 'WID_BASE_0', WID_BASE_0, /NO_BLOCK

end


pro WID_BASE_0_event, Event

  common wid

  common cfg

  Device, decomposed=0

  ;help,event,/str
  wTarget = (widget_info(Event.id,/NAME) eq 'TREE' ?  $
    widget_info(Event.id, /tree_root) : event.id)

  wWidget =  Event.top

  case wTarget of

    ;Selezione utente
    Widget_Info(wWidget, FIND_BY_UNAME='DROPLIST_USERS'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_DROPLIST' )then begin $
        user=userlist(event.index)
    endif
  end


  ;Selezione FITS folder
  Widget_Info(wWidget, FIND_BY_UNAME='Button_BROWSE'): begin
    if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then begin $
      folder=dialog_pickfile(/DIRECTORY)
    Widget_Control, TEXT_OBS_PATH, set_value=folder
    dirlist=file_search(folder,'2*',count=dirnum,/TEST_DIRECTORY,/FULLY_QUALIFY_PATH)
    outname=strarr(dirnum)
    for thisfol=0, dirnum-1 do begin
      dirs=strsplit(dirlist[thisfol],sp,/extract)
      outname[thisfol]='checkfile_'+strcompress(dirs[-1],/remove_all)+'.txt'
    endfor
  endif
end

;Inizio del flagging
Widget_Info(wWidget, FIND_BY_UNAME='Button_START'): begin
  if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then begin $
    if keyword_set(folder) then begin
    list=file_search(dirlist[0]+sp+'2*.fits', count=numb)
    flagged=intarr(numb)
    flagged[*]=0
    if numb gt 0 then begin
      dirindex=0
      start
    endif else begin
      choice=dialog_message("No FITS in first folder!",/INFO,/CENTER)
    endelse
  endif else begin
    choice=dialog_message("No folder was selected!",/INFO,/CENTER)
  endelse
endif
end

;Accetta solo il canale left
Widget_Info(wWidget, FIND_BY_UNAME='Button_LEFT'): begin
  if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then begin
    widget_control, WID_BUTTON_BACK, SENSITIVE=1
    flagga, list[tobeflagged[index]], 0, app, tobeflagged[index]
  endif
end

;Accetta solo il canale right
Widget_Info(wWidget, FIND_BY_UNAME='Button_RIGHT'): begin
  if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then begin
    widget_control, WID_BUTTON_BACK, SENSITIVE=1
    flagga, list[tobeflagged[index]], 1, app, tobeflagged[index]
  endif
end

;Accetta entrambi e valutali "perfetti"
Widget_Info(wWidget, FIND_BY_UNAME='Button_PERFECT'): begin
  if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then begin
    widget_control, WID_BUTTON_BACK, SENSITIVE=1
    flagga, list[tobeflagged[index]], 2, app, tobeflagged[index]
  endif
end

;Accetta entrambi nonostante imperfezioni
Widget_Info(wWidget, FIND_BY_UNAME='Button_ALL'): begin
  if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then begin
    widget_control, WID_BUTTON_BACK, SENSITIVE=1
    flagga, list[tobeflagged[index]], 3, app, tobeflagged[index]
  endif
end

;Rifiuta entrambi i canali
Widget_Info(wWidget, FIND_BY_UNAME='Button_NONE'): begin
  if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then begin
    widget_control, WID_BUTTON_BACK, SENSITIVE=1
    flagga, list[tobeflagged[index]], 4, app, tobeflagged[index]
  endif
end

;Back to previous subscan
Widget_Info(wWidget, FIND_BY_UNAME='WID_BUTTON_BACK'): begin
  if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then begin
    if tobeflagged[index] eq 0 then begin
      choice=dialog_message("Not applicable, this is the first subscan!",/INFORMATION,/CENTER)
      return
    endif
    readcol, dirlist[dirindex]+sp+outname[dirindex], fname, flag, usid, format='(A,A,A)'
    nitems=n_elements(fname)
    openw, Unit, dirlist[dirindex]+sp+outname[dirindex], /get_lun
    for i=0, nitems-2 do begin
      printf, Unit, fname[i], flag[i], usid[i], format='(A,1X,A,1X,A)'
    endfor
    close, Unit
    free_lun, Unit
    index=index-1
    plot_ch
  endif
end

;Skip this subscan, go to the next one
Widget_Info(wWidget, FIND_BY_UNAME='WID_BUTTON_SKIP'): begin
  if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then begin $
    if index eq n_elements(tobeflagged)-1 then begin
    choice=dialog_message("Not applicable: this is the last subscan!",/INFORMATION,/CENTER)
    return
  endif
  index=index+1
  widget_control, WID_BUTTON_BACK, SENSITIVE=0
  plot_ch
endif
end

;Stop and exit
Widget_Info(wWidget, FIND_BY_UNAME='WID_BUTTON_STOP'): begin
  if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
    Widget_Control, WID_BASE_0, /DESTROY
end

;Back to previous folder
Widget_Info(wWidget, FIND_BY_UNAME='WID_BUTTON_PREVDIR'): begin
  if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then begin
    if dirindex eq 0 then begin
      choice=dialog_message("Not applicable, this was the first folder!",/INFORMATION,/CENTER)
      return
    endif
    dirindex=dirindex-1
    list=file_search(dirlist[dirindex]+sp+'2*.fits', count=numb)
    flagged=intarr(numb)
    flagged[*]=0
    if numb gt 0 then begin
      start
    endif else begin
      choice=dialog_message("No FITS in first folder!",/INFO,/CENTER)
    endelse
  endif else begin
    choice=dialog_message("No folder was selected!",/INFO,/CENTER)
  endelse
end


;On to next folder
Widget_Info(wWidget, FIND_BY_UNAME='WID_BUTTON_NEXTDIR'): begin
  if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then begin
    if dirindex eq n_elements(dirlist)-1 then begin
      choice=dialog_message("Not applicable, this was the last folder!",/INFORMATION,/CENTER)
      return
    endif
    dirindex=dirindex+1
    list=file_search(dirlist[dirindex]+sp+'2*.fits', count=numb)
    flagged=intarr(numb)
    flagged[*]=0
    if numb gt 0 then begin
      start
    endif else begin
      choice=dialog_message("No FITS in first folder!",/INFO,/CENTER)
    endelse
  endif else begin
    choice=dialog_message("No folder was selected!",/INFO,/CENTER)
  endelse
end

else:
endcase

end




;
; Empty stub procedure used for autoloading.
;
pro dataflagging, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_
  WID_BASE_0, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_
end

;##############################################################
pro Show_Message, comment

  common wid

  widget_control, WID_LABEL_COMMENT, Set_Value=comment

end
;##############################################################

pro plot_ch
  common cfg
  common wid
  ;how many files to be flagged?
  numb=n_elements(list)
  ; load FITS file
  data=MRDFITS(list[tobeflagged[index]],4,/silent)
  ; extract the filename from full path
  segments=strsplit(list[tobeflagged[index]],sp,/extract)
  segcount=n_elements(segments)
  filename=segments[segcount-1]
  ; count the data samples and generate x-axis
  nsample=n_elements(data.time)
  xaxis=indgen(nsample)
  ; ch0 plot
  wset, wid_left_id
  plot, xaxis[10:nsample-10], data[10:nsample-10].ch0, $
    ys=1, background=255, color=0, title=filename
  ; ch1 plot
  wset, wid_right_id
  plot, xaxis[10:nsample-10], data[10:nsample-10].ch1, $
    ys=1, background=255, color=0, title='Dir'+strcompress(string(dirindex+1))+' out of'+strcompress(n_elements(dirlist))+'  -  FITS'+strcompress(string(tobeflagged[index]+1))+' out of'+strcompress(string(numb))
end


pro start
  common wid
  common cfg

  widget_control, Button_LEFT, SENSITIVE=1
  widget_control, Button_RIGHT, SENSITIVE=1
  widget_control, Button_PERFECT, SENSITIVE=1
  widget_control, Button_ALL, SENSITIVE=1
  widget_control, Button_NONE, SENSITIVE=1
  widget_control, WID_BUTTON_PREVDIR, SENSITIVE=1
  widget_control, WID_BUTTON_BACK, SENSITIVE=1
  widget_control, WID_BUTTON_STOP, SENSITIVE=1
  widget_control, WID_BUTTON_SKIP, SENSITIVE=1
  widget_control, WID_BUTTON_NEXTDIR, SENSITIVE=1

  if file_test(dirlist[dirindex]+sp+outname[dirindex]) then begin
    choice=dialog_message("Checkfile already exists: resume and append? 'NO' means you'll start again and overwrite!!!",/QUESTION,/CENTER)
    if choice eq 'No' then begin
      app='No'
    endif else begin
      readcol, dirlist[dirindex]+sp+outname[dirindex], fname, flag, usid, format='(A,A,A)'
      checklines=n_elements(fname)
      match_index=0
      ; evidenzio i file già flaggati
      for i=0, n_elements(list)-1 do begin
        for j=0, checklines-1 do begin
          if strmatch(list[i],'*'+fname[j]) then flagged[i]=1
        endfor
      endfor
      app='Yes'
      if (total(flagged) ge n_elements(list))  then begin
        choice=dialog_message("All FITS were already flagged. Change folder.",/INFORMATION,/CENTER)
        return
      endif
    endelse
    tobeflagged=where(flagged eq 0, howmany)
  endif else begin
    tobeflagged=where(flagged eq 0, howmany)
    openw, Unit, dirlist[dirindex]+sp+outname[dirindex], /GET_LUN
    close, Unit
    free_lun, Unit
    app='Yes'
  endelse
  Widget_Control, WID_TAB, SET_TAB_CURRENT=1
  index=0
  flagnext=tobeflagged[index]
  plot_ch
end

pro flagga, current, flag_id, appflag, ind
  common cfg
  common wid
  slash=strpos(current,sp,/reverse_search)
  current=strmid(current,slash+1,strlen(current)-slash-1)
  if (appflag eq "No") and (ind eq 0) then begin
    openw, Unit1, dirlist[dirindex]+sp+outname[dirindex], /get_lun
  endif else begin
    openw, Unit1, dirlist[dirindex]+sp+outname[dirindex], /get_lun, /append
  endelse
  printf, Unit1, current, flags[flag_id], user, format='(A,1X,A,1X,A)'
  ; print, 'Writing in ', Unit1
  close, Unit1
  free_lun, Unit1
  if index lt n_elements(tobeflagged)-1 then begin
    index=index+1
    flagnext=tobeflagged[index]
  endif else begin
    finale=dialog_message('No more files to flag!',/INFORMATION,/CENTER)
    widget_control, Button_LEFT, SENSITIVE=0
    widget_control, Button_RIGHT, SENSITIVE=0
    widget_control, Button_PERFECT, SENSITIVE=0
    widget_control, Button_ALL, SENSITIVE=0
    widget_control, Button_NONE, SENSITIVE=0
    widget_control, WID_BUTTON_SKIP, SENSITIVE=0
    widget_control, WID_BUTTON_BACK, SENSITIVE=0
  endelse
  if flagnext lt n_elements(list) then begin plot_ch
endif else begin
  finale=dialog_message('No more files to flag!',/INFORMATION,/CENTER)
  widget_control, Button_LEFT, SENSITIVE=0
  widget_control, Button_RIGHT, SENSITIVE=0
  widget_control, Button_PERFECT, SENSITIVE=0
  widget_control, Button_ALL, SENSITIVE=0
  widget_control, Button_NONE, SENSITIVE=0
  widget_control, WID_BUTTON_SKIP, SENSITIVE=0
  widget_control, WID_BUTTON_BACK, SENSITIVE=0
endelse
end



