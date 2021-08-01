% "Wyceń europejską i amerykańską opcję sprzedaży na akcję nie wypłacającą dywidendy z ceną wykonania K=100 i terminem wygaśnięcia T=3 mies.
% Załóż, że wolna od ryzyka stopa procentowa r=10% rocznie, zmienność cen akcji sigma = 20% rocznie oraz
% że dynamika cen akcji jest opisana drzewkiem multiplikatywnym z u=exp(sigma*dt^0.5) i d=1/u. 
% Sprawdź zbieżność do ceny Blacka-Scholesa dla drzewek n=5, 4, ..., 100 krokowych oraz S_0 = 95, ..., 105."

T = 3;
K = 100;
sigma = 0.2;
r = 0.1;
dt = 1/12;
n = 4:100;

% opcja europejska

figure (1)

for S = 95:1:105
    X = [];
    for i = n
        dtt = dt*T/i;
        X = [X, drzewko_mult(i, S, sigma, K, r, dtt, "put")];
    end
    P = K*exp(-r*T/12)*normcdf((-log(S/K)-(r-sigma^2/2)*T/12)/(sigma*sqrt(T/12))) - S*normcdf((-log(S/K)-(r+sigma^2/2)*T/12)/(sigma*sqrt(T/12)));
    subplot(3, 4, S-94)
    plot(n, X)
    hold on
    plot(n, ones(1, length(n))*P, '--', 'linewidth', 1.2)
    xlabel('n')
    ylabel('cena')
    title("S_0 = " + S)
end

% opcja amerykańska

figure (2)

for S = 95:1:105
    X = [];
    for i = n
        dtt = dt*T/i;
        X = [X, drzewko_mult_a(i, S, sigma, K, r, dtt, "put")];
    end
    P = K*exp(-r*T/12)*normcdf((-log(S/K)-(r-sigma^2/2)*T/12)/(sigma*sqrt(T/12))) - S*normcdf((-log(S/K)-(r+sigma^2/2)*T/12)/(sigma*sqrt(T/12)));
    subplot(3, 4, S-94)
    plot(n, X)
    hold on
    plot(n, ones(1, length(n))*P, '--', 'linewidth', 1.2)
    xlabel('n')
    ylabel('cena')
    title("S_0 = " + S)
end
