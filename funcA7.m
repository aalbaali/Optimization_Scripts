function val = funcA7(x,n)

[n,m] = size(x);
if (m ~= 1)
    if(n ~= 1) %it's not a row matrix
        disp("Not a vector");
        return;
    else
        x = x';  %a row vector is entered instead of column
        %fix it but don't exit
        fprintf("a row vector was entered: [");
        for i=1:m-1
            fprintf("%u, ",x(i));
        end
        fprintf("%u]\n",x(m));
    end
end
   



n = size(x,1);
m = n+2;



for i = 1:n
    @F

        