
%% ******************************************************************** %%
%                 Create a Reference Input ( two inputs )                          
% *********************************************************************** %
function ref = ref_func(t, varargin)
    n = 1;
    if length(varargin) == 1
        n = varargin{1};
    end
    switch(n)
        case 1
           if(0 <= t) && (t < 1)
                ref = 10*cos(2*pi*t);
           elseif (1 <= t) && (t < 2)
                ref = 10*0.5*t^2*(1-t);
           else
                ref = 10*0.5*cos(4*pi*t)+1;
           end
        case 2
           if(0 <= t) && (t < 1)
                ref = 10*[cos(2*pi*t);
                       1.2*t^2*(1-t)];
           elseif (1 <= t) && (t < 2)

                ref = 10*[0.5*t^2*(1-t);
                          cos(2*pi*t)];
           else
                ref = 10*[0.5*cos(4*pi*t)+1;
                       0.2*sin(4*pi*t)-0.5];   
           end
        case 3
           if(0 <= t) && (t < 1)
                ref = 10*[cos(2*pi*t);
                        0.5*t^2*(1-t);
                       1.2*t^2*(1-t)];
           elseif (1 <= t) && (t < 2)

                ref = 10*[0.5*t^2*(1-t);
                           1.2*t^2*(1-t);
                           cos(2*pi*t)];
           else
                ref = 10*[0.5*cos(4*pi*t)+1;
                              1.2*t^2*(1-t);
                       0.2*sin(4*pi*t)-0.5];   
           end
    end
end
