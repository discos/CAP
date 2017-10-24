pro srcaver 

; Program to average the LCP and RCP measurements produced by CAP, performing a source-based weighted mean. 
; Just select runtarget's output file from the graphic interface appearing after launching:
; 
;  IDL> srcaver 
;  
;  Notice: results are written in a file named "Combined_measurements.txt" 
;          In case of multiple runs, the generated files are automatically numbered with a sequential suffix. 
;          The produced file includes the full path to the original file given in input.  
; 
; Last edited: Oct 16, 2017

fin=dialog_pickfile()

close, /all

readcol, fin, src, hhmmss, mjd, freq, elev, A0, A1, eA0, eA1, f0, f1, ef0, ef1, format='(A,D,D,D,D,D,D,D,D,D,D,D,D)'

src_list=src[UNIQ(src, sort(src))]

num=n_elements(src_list)

aver_f0=dblarr(num)
aver_f1=dblarr(num)
aver_ef0=dblarr(num)
aver_ef1=dblarr(num)
combined_f=dblarr(num)
combined_ef=dblarr(num)

outf=file_search('Combined_measurements*.txt', count=numf)
if numf eq 0 then begin
  outfile='Combined_measurements.txt'
endif else begin
  outfile='Combined_measurements_'+strcompress(string(numf+1),/remove_all)+'.txt'
endelse

openw, 1, outfile
printf, 1, 'Elaboration of CAP results contained in: '
printf, 1, fin
printf, 1, ' '
printf, 1, '         Name  Freq     F0   eF0     F1   eF1      F    eF'
printf, 1, '                MHz     Jy    Jy     Jy    Jy     Jy    Jy'
 

for i=0, num-1 do begin 
  
  thissrc=where(src eq src_list[i])
  w0=1./ef0[thissrc]^2
  w1=1./ef1[thissrc]^2
  aver_f0[i]=total(f0[thissrc]/ef0[thissrc]^2)/total(w0)
  aver_f1[i]=total(f1[thissrc]/ef1[thissrc]^2)/total(w1)
  aver_ef0[i]=sqrt(1/total(w0))
  aver_ef1[i]=sqrt(1/total(w1))
 
  combined_f[i]=(aver_f0[i]/aver_ef0[i]^2+aver_f1[i]/aver_ef1[i]^2)/(1./aver_ef0[i]^2+1./aver_ef1[i]^2)
  combined_ef[i]=sqrt(1/(1./aver_ef0[i]^2+1./aver_ef1[i]^2))

  printf, 1, src_list[i], freq[0], aver_f0[i], aver_ef0[i], aver_f1[i], aver_ef1[i], combined_f[i], combined_ef[i], format='(A13,1X,I5,1X,D6.3,1X,D5.3,1X,D6.3,1X,D5.3,1X,D6.3,1X,D5.3)'
 
endfor

close, 1

print, '****** DONE ******' 

return 
end