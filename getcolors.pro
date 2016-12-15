FUNCTION GetColors, Value, BreakVal, BreakName, open=open,      $
    legend=legend, format=format, ctable=ctable, cpart=cpart,   $
    usedcolors=usedcolors, badforeground=badforeground,         $
    badbackground=badbackground, charsize=charsize,             $
    flip=flip, noedge=noedge, ncolors=ncolors, step=step, _extra=_extra

    @compile_opt.pro                        ; On error, return to caller

nBreak = n_elements(BreakVal)   ; # break values between colors

InitVar, badforeground  , /key
InitVar, badbackground  , /key
InitVar, open           , /key
InitVar, legend         , /key
InitVar, flip           , /key
InitVar, noedge         , /key
InitVar, ncolors        , !d.n_colors

CASE noedge OF
0: InitVar, cpart       , [0,1]
1: cpart = [1.0,ncolors-2.0]/(ncolors-1.0)
ENDCASE

IF nBreak GE ncolors THEN message,  $
    '#breaks ('+strcompress(nBreak,/remove_all)+') must be less then #colors ('$
    +strcompress(ncolors,/remove_all)+')'

InitVar, BreakName, BreakVal

nCols = nBreak+1                        ; # colors to be used
norm = cpart*(ncolors-1.0)

CASE open OF
0B: Col =  gridgen(nCols  ,norm=norm)           ; Color index array
1B: Col = (gridgen(nCols+2,norm=norm))[1:nCols] ; Color index array
ENDCASE
Col = round(Col)

IF flip THEN Col = reverse(Col)

BreakVal = BreakVal[sort(BreakVal)]

bad = where(1B-finite(Value))           ; Locate bad values

IF bad[0] NE -1 THEN BEGIN
    badv = Value[bad]                   ; Save all bad values

    ; Set bad values to minimum value. If everything is bad use the lowest break value
    
    ValueMin = min(Value, /nan)         ; Pick up minimum value
    IF NOT finite(ValueMin) THEN ValueMin = BreakVal[0]

    Value[bad] = ValueMin               ; Replace bad values by minimum value
ENDIF

; Set up the colors array by assigning colors to all elements in between
; successive break values.
 
Colors = Col[0]+long(0*Value)           ; Set everything to first color

FOR I=1,nBreak-1 DO BEGIN
    A = where ( BreakVal[I-1] LT Value AND Value LE BreakVal[I] , N )
    IF N NE 0 THEN Colors[A] = Col[I]
ENDFOR

A = where ( Value GT BreakVal[nBreak-1] , N )
IF N NE 0 then Colors[A] = Col[nBreak]

IF bad[0] NE -1 THEN BEGIN
    CASE 1 OF
    badbackground: Colors[bad] = !p.background
    badforeground: Colors[bad] = !p.color
    ELSE         : Colors[bad] = -1
    ENDCASE

    Value [bad] = badv
ENDIF

; Load color table if specified

IF IsType(ctable, /defined) THEN loadct, ctable[0]

usedcolors = Colors[ uniq(Colors) ]

IF legend THEN BEGIN

    ; Draw vertical legend at left side of screen

    DyOff = .02
    Dy = (1.-2*DyOff)/nCols  &  Dx = .015   ; Size of box for each color

    ; Loop from bottom to top of screen

    Xleft = 0.01  &  Xright = Xleft+Dx  &  Ytop = DyOff

    FOR I=0,nBreak DO BEGIN                 ; Loop over nCols = nBreak+1 colors
        Ybottom = Ytop  &  Ytop = Ytop+Dy
        polyfill, color=Col[I], /normal,        $
            [Xleft  ,Xleft,Xright,Xright ],     $
            [Ybottom,Ytop ,Ytop  ,Ybottom]
    ENDFOR

    ; Convert BreakName array to string array. If BreakName is of type float then
    ; the specified FORMAT is used in the conversion

    CASE IsType(BreakName, /generic_int) OF
    1: BreakLabel = strcompress(BreakName, /remove_all)
    0: BEGIN
        IF NOT keyword_set(format) THEN $
            message, 'Keyword FORMAT must be assigned a valid float format value'
        BreakLabel = string(BreakName, format=format)
        BreakLabel = strcompress(BreakLabel, /remove_all)
    END
    ENDCASE

    Xleft = Xright+.75*Dx
    InitVar, charsize, 1.0

    InitVar, Step, 1+nCols/8
    step = step > 1

    FOR I=0,nBreak-1,Step DO BEGIN
        Ybottom = DyOff+(I+1)*Dy-0.0125*charsize
        xyouts, Xleft, Ybottom, /normal, BreakLabel[I], charsize=charsize, _extra=_extra
    ENDFOR

ENDIF

RETURN, Colors  &  END

