
	PRO FXREAD_SUB, FILENAME, DATA, HEADER, OUTDIMS, $
		NANVALUE=NANVALUE,  NOSCALE=NOSCALE, NOUPDATE=NOUPDATE,	$
		ERRMSG=ERRMSG
;+
; NAME: 
;	FXREAD_SUB
; Purpose     : 
;	Read sub array from basic FITS files.
; Explanation : 
;	Read a subarray of the primary array from a disk FITS file.  
; Use         : 
;	FXREAD_SUB, FILENAME, DATA, HEADER, OUTDIMS
; Inputs      : 
;	FILENAME = String containing the name of the file to be read.
;	OUTDIMS = An array containing the limits of each dimension to be
;		returned in the form [[minx1,maxx1],[minx2,maxx2],
;		[minx3,maxx3]...[minxn,maxxn]]
; Outputs     : 
;	DATA	 = Data array to be read from the file.
; Opt. Outputs: 
;	HEADER	 = String array containing the header for the FITS file.
; Keywords    : 
;	NANVALUE = Value signalling data dropout.  All points corresponding to
;		   IEEE NaN (not-a-number) are set to this value.  Ignored
;		   unless DATA is of type float or double-precision.
;	NOSCALE	 = If set, then the output data will not be scaled using the
;		   optional BSCALE and BZERO keywords in the FITS header.
;		   Default is to scale, if and only if BSCALE and BZERO are
;		   present and nontrivial.
;	NOUPDATE = If set, then the optional BSCALE and BZERO keywords in the
;		   optional HEADER array will not be changed.  The default is
;		   to reset these keywords to BSCALE=1, BZERO=0.  Ignored if
;		   NOSCALE is set.
;	ERRMSG   = If defined and passed, then any error messages will be
;		   returned to the user in this parameter rather than
;		   depending on the MESSAGE routine in IDL.  If no errors are
;		   encountered, then a null string is returned.  In order to
;		   use this feature, ERRMSG must be defined first, e.g.
;
;			ERRMSG = ''
;			FXREAD, ERRMSG=ERRMSG, ...
;			IF ERRMSG NE '' THEN ...
;
; Calls       : 
;	GET_DATE, IEEE_TO_HOST, FXADDPAR, FXHREAD, FXPAR, WHERENAN, 
;       FXREAD_SUB_SUB
; Common      : 
;	None.
; Restrictions: 
;	Groups are not supported.
;
; Side effects: 
;	If the keywords BSCALE and BZERO are present in the FITS header, and
;	have non-trivial values, then the returned array DATA is formed by the
;	equation
;
;			DATA = BSCALE*original + BZERO
;
;	However, this behavior can overridden by using the /NOSCALE keyword.
;
;	If the data is scaled, then the optional HEADER array is changed so
;	that BSCALE=1 and BZERO=0.  This is so that these scaling parameters
;	are not applied to the data a second time by another routine.  Also,
;	history records are added storing the original values of these
;	constants.  Note that only the returned array is modified--the header
;	in the FITS file itself is untouched.
;
;	If the /NOUPDATE keyword is set, however, then the BSCALE and BZERO
;	keywords are not changed.  It is then the user's responsibility to
;	ensure that these parameters are not reapplied to the data.  In
;	particular, these keywords should not be present in any header when
;	writing another FITS file, unless the user wants their values to be
;	applied when the file is read back in.  Otherwise, FITS readers will
;	read in the wrong values for the data array.
;	
; Category    : 
;	Data Handling, I/O, FITS, Generic.
; Prev. Hist. : 
;       K. Kuntz, Feb 2004, based on FXREAD by W. Thompson which was, 
;	in turn, based in part on READFITS by W. Landsman, and
;			       STSUB by M. Greason and K. Venkatakrishna.
; Written     : 
;	K. Kuntz, GSFC, Feb 2004.
; Modified    : 
; Version     :
;       Version 1   22-Feb-2004
;-
;
	ON_ERROR, 0
        READ_OK=0
;
;  Get the UNIT number, and open the file.
;
	OPENR, UNIT, FILENAME, /BLOCK, /GET_LUN, ERROR=ERROR
        IF ERROR NE 0 THEN BEGIN
	    MESSAGE='Error opening '+FILENAME
	    IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
		ERRMSG = MESSAGE
		RETURN
	    END ELSE MESSAGE, MESSAGE
        ENDIF
;
;  Read in the FITS header.
;
	FXHREAD,UNIT,HEADER,STATUS
	IF STATUS NE 0 THEN BEGIN
		FREE_LUN,UNIT
		MESSAGE = 'Unable to read FITS header'
		IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
			ERRMSG = MESSAGE
			RETURN
		END ELSE MESSAGE, MESSAGE
	ENDIF
;
;  Extract the keywords BITPIX, NAXIS, NAXIS1, ...
;
	BITPIX = FXPAR(HEADER,'BITPIX')
	NAXIS = FXPAR(HEADER,'NAXIS')
	IF NAXIS EQ 0 THEN BEGIN
		FREE_LUN,UNIT
		MESSAGE = 'File does not contain a primary array'
		IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
			ERRMSG = MESSAGE
			RETURN
		END ELSE MESSAGE, MESSAGE
	ENDIF
	INDIMS = FXPAR(HEADER,'NAXIS*')
;
;  Try to catch an error
;
	IF(NAXIS NE N_ELEMENTS(OUTDIMS)/2)THEN BEGIN
        	FREE_LUN,UNIT
		MESSAGE = 'Input dimension incompatible with primary array'
                IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
                        ERRMSG = MESSAGE
                        RETURN
                END ELSE MESSAGE, MESSAGE
        ENDIF   
	FOR I=0,NAXIS-1 DO BEGIN
		IF((OUTDIMS(0,I) LT 0) OR (OUTDIMS(1,I) GE INDIMS(I)))THEN BEGIN
                	FREE_LUN,UNIT
                	MESSAGE = 'Extraction dimension: '+STRING(I)+$
			' incompatible with primary array'
                	IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
                        	ERRMSG = MESSAGE
                        	RETURN
                	END ELSE MESSAGE, MESSAGE
		ENDIF
	ENDFOR
	DELDIMS=INTARR(NAXIS)
        FOR I=0,NAXIS-1 DO BEGIN
        	DELDIMS(I)=OUTDIMS(1,I)-OUTDIMS(0,I)+1
        ENDFOR
;
;  Determine the array type from the keyword BITPIX.
;
	CASE BITPIX OF
		  8:	IDLTYPE = 1	; Byte
		 16:	IDLTYPE = 2	; Integer*2
		 32:	IDLTYPE = 3	; Integer*4
		-32:	IDLTYPE = 4	; Real*4
		-64:	IDLTYPE = 5	; Real*8
	ENDCASE
;
;  Make the array.
;
	DATA = MAKE_ARRAY(DIMENSION=DELDIMS,TYPE=IDLTYPE,/NOZERO)
;
;  Beginning of the stuff created by kuntz
;  Figure out how many pixels need to be skipped in each dimension
;
;	PRINT,'Input dimensions:'
;       PRINT,INDIMS
	SKIP=ULONG64(INDIMS-DELDIMS)
	SKIP(NAXIS-1)=0	;it will never get used
	FOR I=1,NAXIS-1 DO SKIP(I:NAXIS-1)=SKIP(I:NAXIS-1)*INDIMS(I-1)
	SKIP=SKIP*ABS(BITPIX)/8		;skip now in bytes
;
;  And the offset from the beginning of the data
;
	INDELTA=ULONG64(OUTDIMS(0,0))
	FOR I=1,NAXIS-1 DO BEGIN
		FACTOR=ULONG64(1)
		FOR J=1,I DO FACTOR=FACTOR*INDIMS(J-1)
		INDELTA=INDELTA+OUTDIMS(0,I)*FACTOR
	ENDFOR
	OUTOFFSET = ULONG64(0)
;
;  Find the start of the data to be read in.
;
	POINT_LUN,-UNIT,INOFFSET		;Current position
;	WHILE(INDELTA GT 2^11)DO BEGIN
;		INDELTA = INDELTA - 2^11
;        	INOFFSET = INOFFSET + (2^11)*ABS(BITPIX)/8
;	ENDWHILE
        INOFFSET = INOFFSET + INDELTA*ABS(BITPIX)/8
	POINT_LUN,UNIT,INOFFSET		;Position to beginning of wanted data
	CURLEV = NAXIS-1
	LINE = MAKE_ARRAY(DELDIMS(0),TYPE=IDLTYPE,/NOZERO)
        LLEN = DELDIMS(0)*ABS(BITPIX)/8
	ON_IOERROR,QUIT
	FXREAD_SUB_SUB,UNIT,DATA,CURLEV,SKIP,LINE,LLEN,DELDIMS,INOFFSET,OUTOFFSET
;
;  Convert the data from IEEE to host format, keeping track of any IEEE NaN
;  values.  
;
	IF (N_ELEMENTS(NANVALUE) EQ 1) AND (IDLTYPE GE 4) AND	$
		(IDLTYPE LE 6) THEN W = WHERENAN(DATA,COUNT) ELSE $
		COUNT = 0
	IEEE_TO_HOST,DATA
;
;  If the parameters BZERO and BSCALE are non-trivial, then adjust the array by
;  these values.
;
	IF NOT KEYWORD_SET(NOSCALE) THEN BEGIN
		BZERO  = FXPAR(HEADER,'BZERO')
		BSCALE = FXPAR(HEADER,'BSCALE')
		GET_DATE,DTE
		IF (BSCALE NE 0) AND (BSCALE NE 1) THEN BEGIN
			DATA = BSCALE*DATA
			IF NOT KEYWORD_SET(NOUPDATE) THEN BEGIN
				FXADDPAR,HEADER,'BSCALE',1.
				FXADDPAR,HEADER,'HISTORY',DTE +		$
					' applied BSCALE = '+ STRTRIM(BSCALE,2)
			ENDIF
		ENDIF
		IF BZERO NE 0 THEN BEGIN
			DATA = DATA + BZERO
			IF NOT KEYWORD_SET(NOUPDATE) THEN BEGIN
				FXADDPAR,HEADER,'BZERO',0.
				FXADDPAR,HEADER,'HISTORY',DTE +		$
					' applied BZERO = '+ STRTRIM(BZERO,2)
			ENDIF
		ENDIF
	ENDIF
;
;  Store NANVALUE everywhere where the data corresponded to IEE NaN.
;
	IF COUNT GT 0 THEN DATA[W] = NANVALUE
;
;  Close the file and return.
;
        READ_OK=1
QUIT:   ON_IOERROR,NULL
	FREE_LUN, UNIT
        IF NOT READ_OK THEN BEGIN
	    MESSAGE='Error reading file '+FILENAME
	    IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
		ERRMSG = MESSAGE
		RETURN
	    END ELSE MESSAGE, MESSAGE
	ENDIF
	IF N_ELEMENTS(ERRMSG) NE 0 THEN ERRMSG = ''
	RETURN
	END

