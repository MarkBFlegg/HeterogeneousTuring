function [] = paraSpacePlot(d, beta, eta)
%PARASPACEPLOT Summary of this function goes here
%   Detailed explanation goes here

d= 1/d;
b = beta;
a = eta;

con1 = @(a, b, d) (b-a) < (a+b).^3;
con3 = @(a, b, d) (a+b).^3 < d * (b-a);
con4 = @(a, b, d) (d* (b-a) - (a+b).^3).^2 > 4 * d * (a+b).^4;

%color_1 = [233 196 106]/255;
color_1 = [255 255 0]/255;

color_2 = [38 70 83]/255;

%color_3 = [42 157 143]/255;
color_3 = [0 0 255]/255;

color_4 = [244 162 97]/255;
%color_5 = [231 111 81]/255;
color_5 = [255 0 0 ]/255;

    
clf
a = linspace(0, 2*a, 1000);
b = linspace(0, 2*b, 1000);
[aa, bb] = meshgrid(a, b);

hold on 

pattern = con1(aa, bb, d) .* con3(aa, bb, d) .* con4(aa, bb, d);
pattern = pattern.* reshape(color_1, [1 1 3]);
pattern = pattern + ~con1(aa, bb, d) .* con3(aa, bb, d) .* con4(aa, bb, d) .* reshape(color_5, [1 1 3]);
pattern = pattern + con1(aa, bb, d) .* (~con3(aa, bb, d) | ~con4(aa, bb, d)) .* reshape(color_3, [1 1 3]);



imagesc(a, b, pattern);

%colormap pink;

%{
% condition 1
pattern_1 = con1(aa, bb, d);
ai = edgeDetect(pattern_1);
mask = (ai ~= 0);
b_plot = b(mask);
ai = ai(mask);
a_plot = a(ai);
plot(a_plot, b_plot, "black", "LineWidth", 3);


% condition 4
pattern_4 = con4(aa, bb, d);
ai = edgeDetect(pattern_4);
ai = ai(:, 1);
mask = (ai ~= 0);
b_plot = b(mask);
ai = ai(mask);
a_plot = a(ai);
plot(a_plot, b_plot, "black", "LineWidth", 3);

%}


xlim([a(1) a(end)]);
ylim([b(1) b(end)]);

title(sprintf("Parameter Space for Schnakenberg System\n d = %i", 1/d), "FontSize", 13, "interpreter", "latex");
xlabel("$\eta$", "FontSize", 13, "interpreter", "latex");
ylabel("$\beta$", "FontSize", 13, "interpreter", "latex")

h = zeros(3, 1);
h(1) = plot(NaN,NaN,'s', 'Color', color_1, 'MarkerFaceColor', color_1, "MarkerSize", 15);
h(2) = plot(NaN,NaN,'s', 'Color', color_5, 'MarkerFaceColor', color_5, "MarkerSize", 15);
h(3) = plot(NaN,NaN,'s', 'Color', color_3, 'MarkerFaceColor', color_3, "MarkerSize", 15);
hlegend = legend(h, "Turing patterns", "Unstable without diffusion", "Stable with diffusion", "FontSize", 10, "interpreter", "latex");
set(hlegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));

shg;
end

