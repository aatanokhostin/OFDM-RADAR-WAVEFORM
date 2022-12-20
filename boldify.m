function boldify(fig)
if nargin<1
    fig=gcf; 
end

Pa=get(fig,'Children');

for i=1:length(Pa)
    
    if strcmp(get(Pa(i),'Type'),'axes')
        set(Pa(i),'FontSize',18);
        set(Pa(i),'LineWidth',2);

        set(get(Pa(i),'XLabel'),'FontSize',18);

        set(get(Pa(i),'YLabel'),'FontSize',18);

        set(get(Pa(i),'ZLabel'),'FontSize',18);

        set(get(Pa(i),'Title'),'FontSize',18);
    end
    
    Pc=get(Pa(i),'Children');
    for j=1:length(Pc)
        chtype=get(Pc(j), 'Type');
        if strcmp(chtype(1:4), 'text')
            set(Pc(j), 'FontSize', 18);
        elseif strcmp(chtype(1:4), 'line')
            set(Pc(j), 'LineWidth', 2);
            set(Pc(j), 'MarkerSize', 12);
        end
    end
end
