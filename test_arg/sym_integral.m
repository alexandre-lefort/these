function [res_min, res_max] = sym_integral(expr, val_tab,symbol)
    
    res_min = 0;
    res_max = 0;
    eval_int = eval(subs(expr,symbol,val_tab));
    res_min_tab = zeros(length(val_tab),1);
    res_max_tab = zeros(length(val_tab),1);
    
    %%figure; plot(val_tab,log(eval_int));
    for ii=1:(length(val_tab)-1)
        x1 = val_tab(ii);
        x2 = val_tab(ii+1);
        dx = x2 - x1;
        res1 = eval_int(ii)*dx;
        res2 = eval_int(ii+1)*dx;
        
        res_min = res_min + min(res1, res2);
        res_max = res_max + max(res1, res2);
        res_min_tab(ii) = res_min;
        res_max_tab(ii) = res_max;
    end
    
    %%figure;plot(val_tab,res_min_tab(end-1),'r',val_tab,res_max_tab(end-1),'b');