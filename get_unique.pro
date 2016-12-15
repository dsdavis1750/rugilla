FUNCTION get_unique,list

  indx=0
  unik=list(indx)
  plc=where(list NE unik(indx))
  while(plc(0) NE -1)do begin
    list=list(plc)
    unik=[unik,list(0)]
    plc=where(list NE list(0))
  endwhile

  ord=sort(unik)
  unik=unik(ord)
  return,unik
  end

