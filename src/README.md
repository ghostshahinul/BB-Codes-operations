# QEC signed integer coning with cellulation (MATLAB)

Add this folder to your MATLAB path. Main entry points:

```matlab
[HX_new,HZ_new,data] = coneOneZ_overZ_cellulated(HX,HZ_ldi,zi,4);
[HX_all,HZ_all,allData] = coneManyZ_overZ_cellulated(HX,HZ_ldi,1:size(HZ_ldi,1),4);
```

All arithmetic is over the integers Z. Do not reduce mod 2.
The mapping-cone minus is already absorbed into D1, so assemble with `D1.'`, not `-D1.'`.

Important checks:

```matlab
max(abs(data.F0*data.D1 + HX*data.F1),[],'all')
max(abs(data.D0*data.D1),[],'all')
max(abs(HX_new*HZ_new.'),[],'all')
```

Cellulation target defaults to 4. New cellulation edges have zero columns in F0.
