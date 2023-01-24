clearvars;close all;

vel = readRSF('vel.rr',376,1151);
figure;
imagesc(vel);