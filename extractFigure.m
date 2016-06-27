function im=extractFigure(figName)

% Read in file
f=load([figName '.fig'],'-mat');
im=f.hgS_070000.children.children(1).properties.CData;
