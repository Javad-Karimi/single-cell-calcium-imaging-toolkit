function mytform = register_to_reference(fixed, moving)
% Input:
%   fiexed: reference image
%   moving: the image to be registered

[optimizer, metric] = imregconfig('monomodal');
mytform = imregtform(moving, fixed,'rigid', optimizer, metric);
Rfixed = imref2d(size(fixed));
movingRegistered = imwarp(moving,mytform,'OutputView',Rfixed,'Interp','nearest'); 
movingRegisteredE = imwarp(movingE,mytform,'OutputView',Rfixed,'Interp','nearest'); 

f = imfuse(fixed,moving);
figure; imshow(f,'InitialMagnification','fit');
title('before-registration');
f = imfuse(fixed,movingRegistered);
figure; imshow(f,'InitialMagnification','fit');
title('after-registration');

figure;
ax1 = subplot(1,2,1);
imagesc(fixed), axis square; colormap('gray');
title('fixed');
ax2 = subplot(1,2,2);
imagesc(movingRegistered), axis square; colormap('gray');
title('movingRegistered');
linkaxes([ax1, ax2], 'xy');

figure;
ax1 = subplot(1,2,1);
imagesc(fixedE), axis square; colormap('gray');
title('fixed');
ax2 = subplot(1,2,2);
imagesc(movingRegisteredE), axis square; colormap('gray');
title('movingRegistered');
linkaxes([ax1, ax2], 'xy');

end

