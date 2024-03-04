function c = hamming(s1, s2)
    y1 = double(s1);
    y2 = double(s2);
    c = 0;
    for i=1:length(y1)
        if(y1(i)~=y2(i))
            c=c+1;
        end
    end
    
end