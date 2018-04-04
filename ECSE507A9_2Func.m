% This function is for A9 of ECSE507 W2018

function [objFun, con1, con2] = ECSE507A9_2Func(x)
    gamma = 1;


    x1 = x(1);
    x2 = x(2);

    objFun = -(x1+1)^2-(x2+1)^2;
    con1 = x1^2+x2^2-2;
    con2 = x1-gamma;

end