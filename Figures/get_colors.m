function[colors] = get_colors()

acolors = lines(6);
acolors(1,:) = brighten(acolors(1,:), .2);
acolors(2,:) = brighten(acolors(2,:), -.7);
acolors(3,:) = brighten(acolors(3,:), .5);
acolors(4,:) = brighten(acolors(4,:), .5);
acolors(5,:) = brighten(acolors(5,:), -.5);
acolors(6,:) = brighten(acolors(6,:), .7);

colors = [acolors(6,:); acolors(1:end,:)];

end