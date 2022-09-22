; -------------------------------------------------------------------------------
;                           FCLS Unmixing Tool
;
;DESCRIPTION:
; This program is used for linear spectral unmixing with fully constrained least squared method (FCLS) for multi-spectral image.
;  Compile and Run the program, then set input and output data through GUI.
;  Or use the function of 'FCLS_unmixing' in IDL code
;      Result=FCLS_unmixing(ImgData,endmember, residual=residual, end_names=end_names)
;
; Linear Spectral Mixture Model
;      Mixed_spectrum=f1*endmember_1+f2*endmember_2+...+fn*endmember_n
;      constrains: f1+f2+...+fn=1 & f1, f2,... fn >=0
;      n: number of endmembers
;
;INPUT:
; Multispectral Image: ENVI standard file
; Endmember: csv file (first line is endmember name, first row is band number/wavelength, others are endmember data)
;      For Example
;      [
;         wavelength(nm), vegetation,  impervious,  soil
;             442.6, 0.031106667, 0.18706, 0.101763333
;             551.8, 0.035743333, 0.191026667, 0.13075
;             632, 0.022243333, 0.1835,  0.146583333
;             698.3, 0.04047, 0.19548, 0.167933333
;             896.8, 0.230096667, 0.190093333, 0.193353333
;             1020.6,  0.226776667, 0.17783, 0.195176667
;       ]
;METHODS:
;  Fully Constrained Least Squared Method: Fraction is solved with sum-to-one and nonnegative constrained method
;  FCLS method is 
;  
;
;OUTPUT:
; Fraction data: n+1 bands, the last band is residual error (RMSE).
;
;AUTHOR:
; CHEN Xuehong (chenxuehong@bnu.edu.cn) @ Beijing Normal University, 2022-09.


Function FCLS_unmixing,ImgData, endmember, residual=residual, end_names=end_names

  ;---------------------------------------------------------------
  ;DESCRIPTION:
  ;   This function is used for fully constrained least square linear spectral unmixing
  ;
  ;SYNTAX
  ;   Result=FCLS_unmixing(ImgData,endmember, residual=residual, end_names=end_names)
  ;
  ;RETURN VALUE: Result data of 3-dimension matrix(ns*nl*n_endmember)  (n_endmember is number of endmembers, ns and nl are the size of input image)
  ;
  ;INPUT:
  ;   ImgData: 3-dimension matrix(ns*nl*nb)  (nb is the number of spectral bands)
  ;   endmember: 2-dimension matrix (n*nb) of endmember spectra

  ;KEYWORDS:
  ;   residual: 2-dimension matrix (ns*nl), residual of linear spectral mixture model
  ;   end_names: names of endmember
  ;---------------------------------------------------------------

  t=systime(1)
  sizeX=size(endmember)
  sizeY=size(ImgData)

  n_endmember=sizeX[1]
  nb_endmember=sizeX[2]
  ns=sizeY[1]
  nl=sizeY[2]
  nb=sizeY[3]

  if KEYWORD_SET(endnames) eq 0 then endnames=string(indgen(n_endmember))

  F_result=fltarr(ns,nl,n_endmember); array for fraction result
  residual=fltarr(ns,nl) ;array for residual error

  ;FCLS algorithm
  N_set=2^n_endmember-1
  label=bytarr(n_endmember)
  residual=fltarr(ns,nl)
  residual_temp=fltarr(ns,nl)+10^6
  for kk=1,N_set do begin
    for k=0,n_endmember-1 do label[k]=kk/(2^k) mod 2
    ;print,label
    end_addr=where(label eq 1,end_num)
    E_set=endmember[end_addr,*]
    if end_num eq 1 then begin
      F_LS_SUM_TO_ONE=fltarr(ns,nl)+1
      minvalue=F_LS_SUM_TO_ONE
    endif else begin
      F_LS_SUM_TO_ONE=LS_SUM_TO_ONE(E_set,Imgdata)
      minvalue=min(F_LS_SUM_TO_ONE,dimension=3)
    endelse
    for i=0,ns-1 do residual[i,*]=sqrt(total((reform(Imgdata[i,*,*])-E_set##reform(F_LS_SUM_TO_ONE[i,*,*]))^2,2))
    addr_mask=(residual le residual_temp) and (minvalue ge 0)
    ek=0
    for k=0,n_endmember-1 do begin
      if label[k] eq 1 then begin
        F_result[*,*,k]=F_result[*,*,k]*(1-addr_mask)+addr_mask*F_LS_SUM_TO_ONE[*,*,ek]
        ek=ek+1
      endif else begin
        F_result[*,*,k]=F_result[*,*,k]*(1-addr_mask)
      endelse
    endfor
    residual_temp=residual_temp*(1-addr_mask)+addr_mask*(residual<residual_temp)
  endfor
  residual=residual_temp
  Return, F_result
END

Function LS_SUM_TO_ONE,Endmember,Mixed_Data
  ;---------------------------------------------------------------
  ;DESCRIPTION:
  ;   This function is least square method with sum_to_one constrain
  ;
  ;SYNTAX:
  ;   Result=LS_SUM_TO_ONE(Endmember,Mixed_data)
  ;
  ;RETURN VALUE: Result data of 3-dimension matrix (ns*nl*n)
  ;
  ;INPUT:
  ;   Endmember: 2-dimension matrix (n*nb) of endmember spectra
  ;   Mixed_Data: 3-dimension matrix (ns*nl*nb)
  ;---------------------------------------------------------------
  X=Endmember
  AA=invert(transpose(X)##X)##transpose(X)
  Size=size(Mixed_data)
  ns=size[1]
  nl=size[2]
  nb=size[3]
  sizeX=size(X)
  n_endmember=sizeX[1]
  Z=fltarr(n_endmember)+1
  BB=invert(transpose(X)##X)##transpose(Z)##invert(Z##invert(transpose(X)##X)##transpose(Z))
  F_LS_SUM_TO_ONE=fltarr(ns,nl,n_endmember)
  residual=fltarr(ns,nl)
  for i=0,ns-1 do begin
    mixed_data_line=reform(mixed_data[i,*,*])
    F_LS=AA##mixed_data_line
    ZFU=z##F_LS
    F_LS_SUM_TO_ONE[i,*,*]= F_LS-BB##(ZFU-1)
  endfor
  return,F_LS_SUM_TO_ONE
END


PRO FCLS_Unmixing_tool
  ;Creat top level base

  tlb=widget_base(column=1, title='FCLS Spectral Unmixing', tlb_frame_attr=1)
  mainer=widget_base(tlb,column=1,frame=1)
  Ifbase=widget_base(mainer,row=1,/base_align_center)
  label=widget_label(Ifbase,value='Image Filename:')
  Ipath=widget_text(Ifbase,/editable,xsize=40)
  butt=widget_button(Ifbase,value='Browse...', event_pro='Input_Ifile')

  ;Creat Endmember file widgets
  Efbase=widget_base(mainer,row=1,/base_align_center)
  label=widget_label(Efbase,value='Endmember csv file:')
  Epath=widget_text(Efbase,/editable,xsize=40)
  butt=widget_button(Efbase,value='Browse...',event_pro='Input_Endfile')

  ;Creat output widgets
  OUTPUTbase=widget_base(tlb,/row)
  OUTPUTlabel=widget_label(OUTPUTbase, value='Output Filename:')
  OUTPUTpath=widget_text(OUTPUTbase,/editable,xsize=40)
  OUTPUTbutt=widget_button(OUTPUTbase,value='Browse...', event_pro='Output_event')

  ;Creat OK and cancel buttons
  buttsize=75
  bbase=widget_base(tlb,row=1,/align_center)
  butt=widget_button(bbase,value='OK', xsize=buttsize,event_pro='Running')
  butt=widget_button(bbase,value='Cancel', xsize=buttsize)

  info={Ipath:Ipath,Epath:Epath,output:OUTPUTpath}

  widget_control,tlb,set_uvalue=info
  ;Realize widgets
  wxy= get_screen_size()
  xy = widget_info(tlb,/geo)
  offsetxy = (wxy-[xy.SCR_XSIZE,xy.SCR_YSIZE])/2
  widget_control, tlb,xoffset=offsetxy[0],yoffset=offsetxy[1] ;put the topbase into the mid of window
  widget_control,tlb,/realize
  XManager, 'Input', tlb, /NO_BLOCK
END

PRO output_event,event
  output_file = envi_pickfile(title = 'Output:')
  if (output_file eq '') then return
  widget_control,event.top,get_uvalue=info
  widget_control,info.Output,SET_VALUE=output_file
END

PRO  Input_Ifile,event
  Ifile = envi_pickfile(title='Open the image file', filter='*.*')
  if (Ifile eq '') then return
  widget_control,event.top,get_uvalue=info
  widget_control,info.Ipath,SET_VALUE=Ifile
End

PRO  Input_Endfile,event
  Efile = envi_pickfile(title='Open the Endmember file',filter='*.csv')
  if (Efile eq '') then return
  widget_control,event.top,get_uvalue=info
  widget_control,info.Epath,SET_VALUE=Efile
END

PRO Running,event
  widget_control,event.top,get_uvalue=info
  widget_control,info.Ipath,GET_VALUE=Ifile
  widget_control,info.Epath,GET_VALUE=Efile
  widget_control,info.Output,GET_VALUE=out_name
  if Ifile eq '' || Efile eq '' then begin
    print,'Please input the parameter'
    return
  endif else begin
    print,'Path of input Image file:',Ifile
    print,'Path of input Endmember file:',Efile
    print,'output:',out_name
  endelse

  widget_control,event.top,/destroy

  t=systime(1)
    
  envi_open_file,Ifile,R_fid=fid
  enddata=read_csv(Efile,header=end_names,count=band_num)
  end_num=n_elements(end_names)-1
  end_names=end_names[1:*]
  endmember=fltarr(end_num,band_num)
  for i=0,end_num-1 do endmember[i,*]=enddata.(i+1)
  sizeX=size(endmember)
  n_endmember=sizeX[1]
  nb_endmember=sizeX[2]
  
  envi_file_query, fid, fname=fname, ns=ns, nl=nl,nb=nb,dims=dims
  map_info=envi_get_map_info(fid=fid)
  Imgdata=fltarr(ns,nl,nb)
  if nb ne nb_endmember then begin
    print, 'error: band number of image and endmember are inconsistent!'
    return
  endif

  if nb lt n_endmember then begin
    print, 'error: band number should be larger than number of endmember!'
    return
  endif

  for k=0,nb-1 do begin
    Imgdata[*,*,k]=Envi_Get_Data(Fid=fid, dims=dims, pos=k)
  endfor

  F_result=FCLS_unmixing(ImgData,endmember, end_names=end_names, residual=residual)
  out_data=[[[F_result]],[[residual]]]
  
  envi_write_envi_file,out_data,out_name=out_name,bnames=[end_names,'residual err'],map_info=map_info
  print,'time cost:  ', systime(1)-t, 's'
  
END