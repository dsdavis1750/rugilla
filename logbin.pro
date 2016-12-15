pro logbin

emin=100
emax=300000
nbins = 5
ebegin = 100
eend = emax
ebegin = emin
eend = ebegin*exp((0+1) * alog(emax/ebegin) / nbins)
print,ebegin,eend
ebegin = eend
for i = 1,nbins-1 do begin
eend = ebegin*exp((i+1) * alog(emax/ebegin) / nbins)
print,ebegin,eend
ebegin = eend
endfor

return
end
