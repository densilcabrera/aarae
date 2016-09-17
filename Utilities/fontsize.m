function fontsize
if ismac
    font = 12; %Default font for mac
else
    font = 8; %Default font for windows
end
aarae = gcf;
aaraechildren = get(aarae,'Children');
for i = 1:length(aaraechildren)
    if ~isempty(get(aaraechildren(i),'tag'))
        set(aaraechildren(i),'FontSize',font)
        childrenofchild = get(aaraechildren(i),'Children');
        for j = 1:length(childrenofchild)
            if ~isempty(get(childrenofchild(j),'tag'))
                set(childrenofchild(j),'FontSize',font)
                grandkids = get(childrenofchild(j),'Children');
                for k = 1:length(grandkids)
                    if ~isempty(get(grandkids(k),'tag'))
                        set(grandkids(k),'FontSize',font)
                    end
                end
            end
        end
    end
end
