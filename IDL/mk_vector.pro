function mk_vector,vec,width=width,headlength=headlength,headwidth=headwidth, $
  origin=origin,_extra=_extra

  scl=sqrt(total(vec^2))
  u=vec/scl
  
  dummy=min(abs(vec),dim)
  v1=shift([double(1),0,0],dim)
  a=[v1[1]*u[2]-v1[2]*u[1],v1[2]*u[0]-v1[0]*u[2],v1[0]*u[1]-v1[1]*u[0]]
  b=[ a[1]*u[2]- a[2]*u[1], a[2]*u[0]- a[0]*u[2], a[0]*u[1]- a[1]*u[0]]
  a=reform(a/sqrt(total(a^2)))
  b=reform(b/sqrt(total(b^2)))
  
  hl=1.0      ; head proportion (length/width)
  hw=0.5      ; head width
  w=0.1       ; width of line
  
  if 2*hl*hw gt scl then begin
    hw=scl/hl*hw
    w =scl/hl*w
  endif

  if keyword_set(width) then w=width
  if keyword_set(headwidth) then begin
    hw=headwidth
    if keyword_set(headlength) then hl=headlength/hw
  endif else if keyword_set(headlength) then hw=headlength/hl
  
  mesh_obj,6,vert_in,poly,[[0,0,0],[w/2,0,0],[w/2,0,scl-hl*hw], $
    [hw/2,0,scl-hl*hw],[0,0,scl]],p1=20
  vert=vert_in
  vert=(a # vert[0,*])+(b # reform(vert[1,*]))+(u # reform(vert[2,*]))
  
  if n_elements(origin) eq 3 then for i=0,2 do vert[i,*]=vert[i,*]+origin[i]

  return,obj_new("IDLgrPolygon",data=vert,polygons=poly,_extra=_extra)
end
