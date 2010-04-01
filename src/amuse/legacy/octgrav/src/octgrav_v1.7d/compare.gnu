set log
p[1.0e-4:1][1.0e-4:1]'oct_0.7' u 1:($2/7.63) title 'oct_0.7' w d, 0.5, \
   'quad_0.7' u 1:($2/7.63) title 'quad_0.7'  w d, 0.1, 0.01, \
'oct_0.5' u 1:($2/7.63) title 'oct_0.5' w d lt 1 lw 2, \
'quad_0.5' u 1:($2/7.63) title 'quad_0.5' w d lt 3 lw 2, \
'oct_0.3' u 1:($2/7.63) title 'oct_0.3'  w d lt 1 lw 4, \
'quad_0.3' u 1:($2/7.63) title 'quad_0.3' w d lt 3 lw 4 \

