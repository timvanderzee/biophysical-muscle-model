function f_overlap = return_f_overlap(x, parms)
% Function returns f_overlap
f_overlap = 0;

x_no_overlap = x - parms.myofilaments.thick_filament_length;
x_overlap = parms.myofilaments.thin_filament_length - x_no_overlap;
max_x_overlap = parms.myofilaments.thick_filament_length -  ...
    parms.myofilaments.bare_zone_length;

if (x_overlap<0)
    f_overlap=0;
end

if ((x_overlap>0)&(x_overlap<=max_x_overlap))
    f_overlap = x_overlap/max_x_overlap;
end

if (x_overlap>max_x_overlap)
    f_overlap=1;
end

protrusion = parms.myofilaments.thin_filament_length - ...
    (x + parms.myofilaments.bare_zone_length);

if (protrusion > 0)
    x_overlap = (max_x_overlap - protrusion);
    f_overlap = x_overlap / max_x_overlap;
end

end