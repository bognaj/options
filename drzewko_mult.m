function y = drzewko_mult(T, S, sigma, K, r, dt, call_put)
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

    % dla drzewka multiplikatywnego prawdopodobieństwo q jest stałe
    q = (exp(r*dt) - d)/(u - d);
    
    % tworzymy drzewko na cenę opcji
    Xt = zeros(h, l);
    
    if call_put == "call"
        f = @(x) max(x - K, 0);
        Xt(1:2:h, end) = f(St(1:2:h, end));
    elseif call_put == "put"
        f = @(x) max(K - x, 0);
        Xt(1:2:h, end) = f(St(1:2:h, end));
    elseif call_put == "lookback"
        for k = 1:2:h
            n = (k+1)/2;
            min_val = St(k, end) * d^(l-n);
            Xt(k, end) = St(k, end) - min_val;
        end
    end
    
    for col = l-1:-1:1
        for row = (l-col)+1:2:h-(l-col)
            Xt(row, col) = exp(-r*dt) * (q*Xt(max(1, row-1), col+1) + (1-q)*Xt(min(row+1, h), col+1));
        end
    end
    St
    q
    Xt
    y = Xt(l, 1);
end