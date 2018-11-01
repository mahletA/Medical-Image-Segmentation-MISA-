%% Plotting Dice indices


sim = dice(seg2, truth);

imshowpair(seg2 ==1 , truth==1)
title(['CSF Dice Index = ' num2str(sim(1))])
saveas(gcf, '../Results/csf5.png');

imshowpair(seg2 ==2 , truth==2)
title(['GM Dice Index = ' num2str(sim(2))])
saveas(gcf, '../Results/gm5.png');

imshowpair(seg2 ==3 , truth==3)
title(['WM Dice Index = ' num2str(sim(3))])
saveas(gcf, '../Results/wm5.png');