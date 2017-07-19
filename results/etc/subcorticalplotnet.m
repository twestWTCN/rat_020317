close all; clear all
% 1 is control; 2 is lesion
lowbeta(:,:,1) = [0  0.06    0.13    0.08;
    0.06  0   0.14    0.1;
    0.04    0.06    0   0.07;
    0.05    0.06    0.17    0];
lowbeta(:,:,2) = [0 0.14 0.14 0.14;
    0.06    0   0.16    0.1;
    0.14    0.05    0   0.04;
    0.05    0.06    0.16    0];

highbeta(:,:,1) = ...
    [0  0.15    0.19    0.16;
    0.12    0   0.18    0.1;
    0.08    0.08    0   0.05;
    0.1 0.05    0.17    0];
highbeta(:,:,2) = ...
    [0  0.06    0.1 0.05;
    0.06    0   0.12    0.08;
    0.05    0.05    0   0.06;
    0.08    0.06    0.1 0];
subplot(2,3,1)
imagesc(squeeze(lowbeta(:,:,1)))
subplot(2,3,2)
imagesc(squeeze(lowbeta(:,:,2)))
subplot(2,3,3)
imagesc(squeeze(lowbeta(:,:,2)-lowbeta(:,:,1)))
subplot(2,3,4)
imagesc(squeeze(highbeta(:,:,1)))
subplot(2,3,5)
imagesc(squeeze(highbeta(:,:,2)))
subplot(2,3,6)
imagesc(squeeze(highbeta(:,:,2)-highbeta(:,:,1)))
