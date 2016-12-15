!x.title='!17'
!y.title='!17'
plot,[0,1],[0,1]
xmm_src2reg,'EMOS1','mos1-eml.fits','ccf.cif','mlf.ps',1500.,0.8,40.,'sky','sources_sky.reg'
xmm_src2reg,'EMOS1','mos1-eml.fits','ccf.cif','mlf.ps',1500.,0.8,40.,'det','sources_det1.reg'
spawn,'mv sas.txt mos1-region.txt'
xmm_src2reg,'EMOS2','mos2-eml.fits','ccf.cif','mlf.ps',1500.,0.8,40.,'det','sources_det2.reg'
spawn,'mv sas.txt mos2-region.txt'
exit
