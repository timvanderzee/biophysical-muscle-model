clear all; close all; clc

x = linspace(3e2,2e3);

parms.myofilaments.thick_filament_length = 1530/2;
parms.myofilaments.thin_filament_length = 1120;
parms.myofilaments.bare_zone_length = 80;

for i = 1:length(x)
    parms.hs_length = x(i);

    f_overlap(i) = return_f_overlap(x(i), parms);
end

%%
figure(1)
plot(x, f_overlap, 'linewidth', 2)
ylim([0 1.1])

xline(parms.myofilaments.thin_filament_length + parms.myofilaments.bare_zone_length,'r--')
xline(parms.myofilaments.thin_filament_length - parms.myofilaments.bare_zone_length,'r--')

xline(1300,'k--')

%%

function f_overlap = return_f_overlap(x, parms)
% Function returns f_overlap

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