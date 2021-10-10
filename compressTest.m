file = "prague-astronomical-clock-detail-871291743639AGq.jpg";

figure
subplot_i = 1;
for Nretain = [ 5, 10, 50, 100 ]
    img = SVDcompress(file, Nretain);
    
    subplot(2,2,subplot_i)
    imagesc(img)
    axis equal tight off
    title(['N = ', num2str(Nretain)])
    
    subplot_i = subplot_i + 1;
end


file = "boat-in-caribbean-14884763094mZ.jpg";

figure
subplot_i = 1;
for Nretain = [ 5, 10, 50, 100 ]
    img = SVDcompress(file, Nretain);
    
    subplot(2,2,subplot_i)
    imagesc(img)
    axis equal tight off
    title(['N = ', num2str(Nretain)])
    
    subplot_i = subplot_i + 1;
end