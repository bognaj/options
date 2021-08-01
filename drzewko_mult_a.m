function y = drzewko_mult_a(T, S, sigma, K, r, dt, call_put)
    % wyliczenie u i d
    u = exp(sigma*sqrt(dt));
    d = 1/u;
    
    % wymiary macierzy do przechowywania drzewka z cenami akcji
    h = 2*T+1;
    l = T+1;
    % tworzenie drzewka z cenami akcji
    St = zeros(h, l);
    St(l, 1) = S;
    for col = 2:l
        for row = 1:h
            if St(min(h, row+1), col-1) > 0 && St(max(1, row-1), col-1) == 0
                St(row, col) = u*St(min(h, row+1), col-1);
            elseif St(min(h, row+1), col-1) == 0 && St(max(1, row-1), col-1) > 0
                St(row, col) = d*St(max(1, row-1), col-1);
            elseif St(min(h, row+1), col-1) > 0 && St(max(1, row-1), col-1) > 0
                St(row, col) = u*St(min(h, row+1), col-1);
            end
        end
    end

    % tworzenie drzewka prawdopodobieństw
    
    hq = 2*T - 1;
    lq = T;
    
    f = exp(r*dt);
    qt = zeros(hq, lq);
    
    for col = 1:lq
        for row = (lq-col)+1:2:hq-(lq-col)
            qt(row, col) = (f*St(min(row+1, h), col) - St(min(row+2, h), min(col+1, l)))/(St(row, min(col+1, l)) - St(min(row+2, h), min(col+1, l)));
        end  
    end
    
    if call_put == "call"
        f = @(x) max(x - K, 0);
    else
        f = @(x) max(K - x, 0);
    end
    
    % tworzymy drzewko na cenę opcji
    Xt = zeros(h, l);
    Xt(1:2:h, end) = f(St(1:2:h, end));  
    
    for col = l-1:-1:1
        for row = (l-col)+1:2:h-(l-col)
            Xt(row, col) = max(exp(-r*dt) * (qt(row-1, col)*Xt(max(1, row-1), col+1) + (1-qt(row-1, col))*Xt(min(row+1, h), col+1)), f(St(row, col)));
        end
    end
    
    St
    qt
    Xt
    y = Xt(l, 1);
end