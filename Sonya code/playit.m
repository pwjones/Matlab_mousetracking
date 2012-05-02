

function playit(mouse,start,width,height)

figure;
for aa=start
    hold off
    plot(mouse(1,aa).cdata(:,2),mouse(1,aa).cdata(:,1),'k.','MarkerSize',15);
    hold on
    plot(mean(mouse(1,aa).cdata(:,2)),mean(mouse(1,aa).cdata(:,1)),'ro');
    xlim([0 width]);
    ylim([0 height]);
    
    set(gca,'YDir','reverse')
    pause(0.005);
    
end

end














