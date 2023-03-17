figure

mestAll = [mEstLM_AW3_mEst_1417(2, 7:end-5),mEstLM_AW3_mEst_1941(2, 7:end-5)];
fcAll = [SO2D_1417(5:end-5); SO2D_1941(5:end-5)]; 
plot(fcAll, mestAll, 'o')

figure
mestAll = [mEstLM_AW3_mEst_1417(2, 7:end-5),mEstLM_AW3_mEst_1941(2, 7:end-5)];
fcAll = [nSO2D_1417(5:end-5); nSO2D_1941(5:end-5)]; 
plot(fcAll, mestAll, 'o')

figure
mestAll = [mEstLM_AW3_mEst_1417(2, 7:end-5)+5.56e14,mEstLM_AW3_mEst_1941(2, 7:end-5)+2.6921e16];
fcAll = [nSO2D_1417(5:end-5); nSO2D_1941(5:end-5)]; 
plot(fcAll, mestAll, 'o')


figure
mestAll = [mEstLM_AW3_mEst_1417(2, 7:end-5)+5.56e14,mEstLM_AW3_mEst_1941(2, 7:end-5)+2.6921e16+offset];
fcAll = [nSO2D_1417(5:end-5); nSO2D_1941(5:end-5)]; 
plot(fcAll, mestAll, 'o')

offset = 7.34e17 
SO2D_1941(5:end-5)

SO2D_1417(5:end-5)



p = plot(SO2D_1941(5:end-5), mEstLM_AW3_mEst_1941(2, 7:end-5), 'o'); hold on; 
p.MarkerFaceColor = markerfacecolor_AW3; 
p = plot(SO2D_1417(5:end-5), mEstLM_AW3_mEst_1417(2, 7:end-5), 'om'); 

p = plot(SO2D_1941(5:end-5), mEstLM_AW3_mEst_1941(2, 7:end-5)+2.6921e16, 'o'); hold on; 
p.MarkerFaceColor = markerfacecolor_AW3; 
p = plot(SO2D_1417(5:end-5), mEstLM_AW3_mEst_1417(2, 7:end-5)+5.56e14, 'om'); 




p.MarkerFaceColor = markerfacecolor_AW3; 
p.MarkerEdgeColor = 'k'; 
xlabel('DOAS SO_2 fit coefficient')
ylabel('mEst_S_O_2')
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w'); Fig.Box = 'on';
p = plot(SO2D_1941(23), mEstLM_AW3_mEst_1941(2, 25), 'ok'); 
p.MarkerEdgeColor = 'k'; 
skyoffset= 5.2511317e17+2.6921e16;
yline(skyoffset)
p = plot(nSO2D_1941(5:end-5), mEstLM_AW3_mEst_1941(2, 7:end-5)+skyoffset, 's'); hold on; 
p.MarkerFaceColor = markerfacecolor_AW3; 
p.MarkerEdgeColor = 'k'; 



p = plot(SO2D_1417(5:end-5), mEstLM_AW3_mEst_1417(2, 7:end-5), 'o'); hold on; 
p.MarkerFaceColor = 'm'; 
p.MarkerEdgeColor = 'k'; 
xlabel('DOAS SO_2 fit coefficient')
ylabel('mEst_S_O_2')
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w'); Fig.Box = 'on';
p = plot(SO2D_1417(23), mEstLM_AW3_mEst_1417(2, 25), 'ok'); 
p.MarkerEdgeColor = 'k'; 
skyoffset= 5.2511317e17+5.56e14;
yline(skyoffset)
p = plot(nSO2D_1417(5:end-5), mEstLM_AW3_mEst_1417(2, 7:end-5)+skyoffset, 's'); hold on; 
p.MarkerFaceColor = 'm'; 
p.MarkerEdgeColor = 'k'; 











plot(SO2D_1941, mEstLM_AW3_mEst_1941(2, 3:end), 'ok') ; hold on 




figure; plot(SO2D_1941, mEstLM_AW3_mEst_1941(2, 3:end), 'ok') ; hold on 
xlabel('DOAS fit coefficient')
ylabel('mEst_S_O_2')
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';


plot(SO2D_1941(5:end-5), mEstLM_AW3_mEst_1941(2, 7:end-5), 'om')
xlabel('DOAS fit coefficient')
ylabel('mEst_S_O_2')
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';


plot(SO2D_1941(23), mEstLM_AW3_mEst_1941(2, 25), 'oy')


plot(nSO2D_1941(5:end-5), mEstLM_AW3_mEst_1941(2, 7:end-5), 'xm')

plot(nSO2D_1941(5:end-5), mEstLM_AW3_mEst_1941(2, 7:end-5)-mEstLM_AW3_mEst_1941(2, 25), 'xg')

yline(mEstLM_AW3_mEst_1941(2, 1)