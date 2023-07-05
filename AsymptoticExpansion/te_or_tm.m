
function [te_pks,tm_pks]=te_or_tm(pks)
    odds=pks(1:2:end);
    evens=pks(2:2:end);  

    d1=pks(2)-pks(1);
    d2=pks(3)-pks(2);

    if (d1>d2) 
        te_pks=odds;
        tm_pks=evens;
    else
        tm_pks=odds;
        te_pks=evens;
    end

   
end