function symbol_obj,type,colour,position

  radius=0.1
  scl=((8*!dpi)/(3*sqrt(6)))^(1/3.)

  if type eq 3 then begin
    mesh_obj,6,vert,poly,[[0,0,1],[sqrt(2),0,1],[sqrt(2),0,-1],[0,0,-1]] $
      *radius*0.5*(4*!dpi/3.)^(1/3.),p1=4,p4=!dpi/4.,p5=9/4.*!dpi
  endif else if type eq 2 then begin
    ;mesh_obj,6,vert,poly,[[0,0,sqrt(3)-sqrt(2)],[1,0,sqrt(3)-sqrt(2)], $
    ;  [0,0,sqrt(3)]]*radius*scl
    mesh_obj,6,vert,poly,[[0,0,-sqrt(2)/4],[1,0,-sqrt(2)/4], $
      [0,0,3/4.*sqrt(2)]]*radius*scl
  endif else begin
    mesh_obj,4,vert,poly,replicate(radius,21,21)
  endelse
  if n_elements(position) ge 3 then for j=0,2 do $
    vert[j,*]=vert[j,*]+position[j]
  return,obj_new("IDLgrPolygon",data=vert,polygons=poly,color=colour, $
    /shading)
end
