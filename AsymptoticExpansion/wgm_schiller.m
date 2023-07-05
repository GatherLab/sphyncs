%Asymptotic expansion of morphological resonance frequencies in Mie
%scattering (S. Schiller, Applied Optics, 1993 Vol 32)
%Isla Barnard May 2016
%Changed by SCaixeiro 2023- cleaned up code corrected typo
% given mode number RANGE, refractive index contrast NUMBER, a
%nd mode type (polarisation: TE or TM)

%returns (ALTERED) RESONANCE SIZE PARAMETER: 'n_cell/RSP'; given fixed
%n_bead

%accepts -vector of MODE RANGE
%        -double precision refractive index contrast, and refractive index
%        of the bead

function isla_rsp=wgm_schiller(all_modes,m,mode_type)

z=-2.3381074105;  %First AIRY ZERO (NIST)

switch mode_type
    case 'TE'
        % e_k cofficients zero
        p=1;
        %% e coefficients
        e_4 = @(m) 0;
        e_5 = @(m) 0;
        e_6 = @(m,z) 0;
        e_7 = @(m) 0;
        e_8 = @(m,z) 0;
    case 'TM' % e_k coefficients non zero
        p=1/(m^2);
        %% e coefficients
        e_k= @(m) (m^2 -1);
        
        e_4=@(m) e_k(m) * (-8 + (12 * m^4) + m^8)/m^8;
        e_5=@(m) e_k(m) * (7000 * m^(-6)) * (- 28 - (m^2) + (56 * m^4) - (16 * m^6) -(7 * m^8) + (2 * m^10));
        e_6=@(m,z) e_k(m) * m^(-8) * ( (5 * (-200 - (32 * m^2) + (526 * m^4) - (226 * m^6) -(99 * m^8) + (62 * m^10) + (4 * m^12))) + (2 * (-400 + (272 * m^2) + (744 * m^4) - (424 * m^6) - (366 * m^8) - (2 * m^10) + (m^12))*z^3));
        e_7=@(m) e_k(m) * (-269500 * m^(-8)) * (-232 + (160 * m^2) + (543 * m^4) - (477 * m^6) - (186 * m^8) + (165 * m^10) - (15 * m^12) + (4 * m^14));
        e_8=@(m,z) e_k(m) * m^(-10) * z * ((-10 * (-459200 + (286000 * m^2) + (1360312 * m^4) - (1305476 * m^6) - (433952 * m^8) + (717562 * m^10) - (209039 * m^12) - (21542 * m^14) + (7060 * m^16))) + (3 * (336000 - (441600 * m^2) - (626496 * m^4) + (891008 * m^6) + (306416 * m^8) - (505696 * m^10) - (72488 * m ^12) - (7664 * m^14) + (2395 * m^16))*(z^3)));
end %TE, TM have different uk

%% d_coefficients
d_0=@(p) -p;
d_1=@(m,z) ( (2^(1/3)) * 3 * ((m^2) -1) * (z^2) ) / (20 * m);
d_2=@(p,m,z) (-( 2^(2/3) * m.^2 * p * ( -3 + (2 * p.^2) )* z ) / 6);
d_3=@(p,m,z) (( 350 * m^4 * (1-p) * p * (-1 + p+ p^2) )+ ( (m^2 - 1)^2 * (10+ z^3)))/ ( 700 * m );
d_4=@(m,z) (-(2^(1/3) * m^2 * z^2 * (4 - m^2+ e_4(m)))/20);
d_5=@(m,z) ( z * (40 * (-1 + (3 * m^2) - (3 * m^4) + (351 * m^6)) - 479*(((m^2 -1)^3) * z^3) - e_5(m)))  /     (2^(4/3)* 63000*m);
d_6=@(m,z) ((5 * m^2 * (-13 - (16 * m^2) + (4 * m^4))) + (2 * m^2 * (128 - (4*m^2) + m^4) * z^3) - e_6(m,z)) / 1400;
d_7=@(m,z) (z^2*((100*(-551 + (2204 * m^2) - (3306 * m^4) - (73256 * m^6) + (10229 * m^8) )) - (20231*((m^2 -1)^4)*z^3) + e_7(m) ))/((2^(2/3))*16170000*m);
d_8=@(m,z) (m^2 * z * ( (10 * (11082 + (44271 * m^2) - (228 * m^4) + (7060 * m^6))) - ( 3 * (52544 + (48432 * m^2) -(11496*m^4) + (2395 * m^6)) * z^3)) + e_8(m,z))/((2^(10/3))* 141750);

d=zeros(1,9);
d(1)=d_0(p);
d(2)=d_1(m,z);
d(3)=d_2(p,m,z);
d(4)=d_3(p,m,z);
d(5)=d_4(m,z);
d(6)=d_5(m,z);
d(7)=d_6(m,z);
d(8)=d_7(m,z);
d(9)=d_8(m,z);

%%
n=all_modes;

d_expanded=repmat(d,size(all_modes,2),1); %expanded vector of d coefficients
mu=@(n) (n + (1/2));

u=mu(n);    %mus

t=linspace(0,8,9); %EXPANSION up to 8 coefficients
[T,U]=meshgrid(t,u);

c=@(x,u,m) (power(u,x./3)).* power((m^2 -1),((x+1)./2));

coeff=c(T,U,m);

isla_rsp= ((u./m) - ((z/m)*((u./2).^(1/3)))) + sum(d_expanded'./coeff');


end
