function makePretty

FontName = 'Helvetica';
set(groot, 'defaultAxesFontName', FontName);
set(groot, 'defaultTextFontName', FontName);

FontSize = 30;
set(groot, 'defaultAxesFontSize', FontSize);
set(groot, 'defaultTextFontSize', FontSize);

LineWidth = 3;
set(groot, 'defaultAxesLineWidth', LineWidth);
set(groot, 'defaultLineLineWidth', LineWidth);

MarkerSize = 10;
set(groot, 'defaultLineMarkerSize', MarkerSize);

set(groot, 'defaultAxesBox', 'on');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultColorbarTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultPolaraxesTickLabelInterpreter', 'latex');