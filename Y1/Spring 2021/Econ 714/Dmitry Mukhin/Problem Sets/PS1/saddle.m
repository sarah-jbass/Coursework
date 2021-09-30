function saddlepath = saddle(kgrid, cgrid, kss, css, D, psigma, palpha, pbeta, pdelta,ns)
ng = length(kgrid);
saddlepath = 0*kgrid; %Pre-allocates an empty grid for the saddle path
for i=1:ng
    K = kgrid(i) + 0*cgrid;
    C = cgrid; % feed these K and C into the iterate function
    for t=1:ns % begin iterating forwards
        [K,C] = iterate(K, C, D, psigma, palpha, pbeta, pdelta); %feed new K and C values into iterate function
    end
    [~,ind] = min(abs(C-css)+abs(K-kss));
    saddlepath(i) = C(ind);
end
end
