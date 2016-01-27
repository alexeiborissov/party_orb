PRO SCATTER3D, PostScript=postscript

    ; Set the PostScript keyword to send draw the plot in a PostScript file
    ; instead of on the display.
    IF Keyword_Set(postscript) THEN PS_Start, FONT=1, Charsize=3.0
    
    ; Create the random data. Set the seed so you see what I see.
    seed = 1L
    x = RANDOMU(seed, 32)
    y = RANDOMU(seed, 32)
    z = EXP(-3 * ((x - 0.5)^2 + (y - 0.5)^2))
    
    ; Load a color table and create colors for the scatterplot.
    cgLoadCT, 33
    zcolors = BYTSCL(z)
    
    ; Set the 3D coordinate space with axes.
    cgSurf, DIST(5), /NODATA, /SAVE, XRANGE=[0,1], $
       YRANGE=[0,1], ZRANGE=[0, 1], XSTYLE=1, $
       YSTYLE=1, ZSTYLE=1, CHARSIZE=2.0, $
       POSITION=[0.1, 0.1, 0.95, 0.95, 0.1, 0.95], $
       XTICKLEN=1, YTICKLEN=1, XGRIDSTYLE=1, YGRIDSTYLE=1
    cgAXIS, XAXIS=1, /T3D, CHARSIZE=2.0
    cgAXIS, YAXIS=1, /T3D, CHARSIZE=2.0
    
    ; Plot the random points in 3D space with a filled circle shape.
    phi = Findgen(32) * (!PI * 2 / 32.)
    phi = [ phi, phi(0) ]
    cgPlotS, x, y, z, PSYM=16, COLOR=zcolors, SYMSIZE=2.5, /T3D
    
    ; Connect the data points to the XY plane of the plot.
    FOR j=0,31 DO cgPlotS, [x(j), x(j)], [y(j), y(j)], [0, z(j)], $
       COLOR=zcolors(j), /T3D
    
    ; Close the PostScript file and clean-up, if required.
    IF Keyword_Set(postscript) THEN PS_End;, /PNG
        
END