function labnh,lon,lat
  common labnh_common,labnh,vmin,dv,bzero,bscale,h

  if (n_elements(labnh) eq 0) then begin
     ;; read in map data
     labnh=mrdfits('labnh.fit',0,h)
     vmin=fxpar(h,"CRVAL3")
     dv=fxpar(h,"CDELT3")
     bzero=fxpar(h,"BZERO")
     bscale=fxpar(h,"BSCALE")
  end 

  ;; Galaxy: integrate from -400 to 400 km/s
  zmin=fix((-400000-vmin)/dv)
  zmax=fix((+400000-vmin)/dv)-1

  adxy,h,lon,lat,pixx,pixy


  ;; should implement nearest neighbor sampling here
  nh=1.82d15*dv*total(bzero+bscale*labnh[pixx,pixy,zmin:zmax])

  return,nh

end 

