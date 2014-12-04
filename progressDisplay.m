function progressDisplay( current_index ,total )

%When 100% removes previous 3 chars.
if current_index==total
    fprintf('\b\b\b')
    return;
end

percentage =  floor(100*current_index/total);
if percentage>floor(100*(current_index-1)/total)
    if percentage==1
    elseif percentage<11
        fprintf('\b\b')
    else
        fprintf('\b\b\b')
    end
    if percentage<100
        fprintf(num2str(percentage));
        fprintf('%s','%');
    end
end

end