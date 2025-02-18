from pyaxidraw import axidraw
ad = axidraw.AxiDraw()
ad.interactive()
ad.options.units = 2
ad.options.pen_pos_up = 100
ad.options.pen_pos_down = 0
print(f'Connected to AxiDraw: {ad.connect()}')
