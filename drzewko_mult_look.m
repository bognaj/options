function y = drzewko_mult_look(T, S, sigma, r, dt)
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
    
    % znajdujemy minima "rodziców" dla pierwszych dwóch kroków
    if T >= 2
        minima = cell(h, l);
        minima{l, 1} = S;
    
        minima{l+1, 2} = St(l+1, 2);
        minima{l-1, 2} = S;
    end
    
    % oraz dla następnych
    for col = 3:l
        for row = 1:h
            a = min(min(minima{max(row-1, 1), col-1}), St(row, col));
            b = min(min(minima{min(row+1, h), col-1}), St(row, col));
            minima{row, col} = [a, b];
        end
    end
    
    % wyliczamy wypłatę w momencie T 
    
    Xt = cell(h, l);
    
    for i = 1:2:h
        Xt{i, end} = St(i, end) - minima{i, end};
    end
    
    % wyceniamy opcję korzystając z macierzy minimów oraz wyliczonego
    % prawdopodobieństwa q
    for col = l-1:-1:1
        for row = (l-col)+1:2:h-(l-col)
            if length(Xt{max(1, row-1), col+1}) == 2 && length(Xt{min(row+1, h), col+1}) == 2
                Xt{row, col}(1) = exp(-r*dt) * (q*Xt{max(1, row-1), col+1}(2) + (1-q)*Xt{min(row+1, h), col+1}(1));
                Xt{row, col}(2) = exp(-r*dt) * (q*Xt{max(1, row-1), col+1}(2) + (1-q)*Xt{min(row+1, h), col+1}(2));
            elseif length(Xt{max(1, row-1), col+1}) == 2 && length(Xt{min(row+1, h), col+1}) == 1
                Xt{row, col} = exp(-r*dt) * (q*Xt{max(1, row-1), col+1}(2) + (1-q)*Xt{min(row+1, h), col+1}(1));
            elseif length(Xt{max(1, row-1), col+1}) == 1 && length(Xt{min(row+1, h), col+1}) == 2
                Xt{row, col} = exp(-r*dt) * (q*Xt{max(1, row-1), col+1}(1) + (1-q)*Xt{min(row+1, h), col+1}(1));
            else
                Xt{row, col} = exp(-r*dt) * (q*Xt{max(1, row-1), col+1}(1) + (1-q)*Xt{min(row+1, h), col+1}(1));
            end
        end
    end
    
    y = Xt{l, 1};
end
    
    
    
    