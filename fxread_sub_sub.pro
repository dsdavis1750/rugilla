PRO FXREAD_SUB_SUB,UNIT,DATA,CURLEV,SKIP,LINE,LLEN,DELDIMS,INOFFSET,OUTOFFSET

;+
; NAME:
;       FXREAD_SUB_SUB
; Purpose     :
;       A subroutine required by FXREAD_SUB_SUB
; Use         :
;       FXREAD_SUB_SUB,UNIT,DATA,CURLEV,SKIP,LINE,LLEN,DELDIMS,INOFFSET,OUTOFFSET
; Inputs      :
;       UNIT = the previously opened unit for the FITS file
;	DATA = the entire output data array
;	CURLEV = dimension of array currently being read
;	SKIP = array showing how many bytes to be skipped for each dimension
;	LINE = array into which to read each vector
;	LLEN = length of LINE in bytes
;	DELDIMS = dimensions of input - dimensions of output
;	INOFFSET = filepointer to position in input file
;	OUTOFFSET = "pointer" to current position in output array
; Calls       :
;       Itself
; Common      :
;       None.
; Written     :
;       K. Kuntz, GSFC, Feb 2004.
; Modified    :
; Version     :
;       Version 1   22-Feb-2004

	IF(CURLEV EQ 0)THEN BEGIN
		READU,UNIT,LINE
		DATA(OUTOFFSET:OUTOFFSET+DELDIMS(0)-1)=LINE
		INOFFSET=INOFFSET+SKIP(0)+LLEN
		POINT_LUN,UNIT,INOFFSET
                OUTOFFSET=OUTOFFSET+DELDIMS(0)
	ENDIF ELSE BEGIN
		FOR I=0,DELDIMS(CURLEV)-1 DO BEGIN
		FXREAD_SUB_SUB,UNIT,DATA,CURLEV-1,SKIP,LINE,LLEN,DELDIMS,INOFFSET,OUTOFFSET
		ENDFOR
		INOFFSET=INOFFSET+SKIP(CURLEV)
		POINT_LUN,UNIT,INOFFSET
	ENDELSE

	END
