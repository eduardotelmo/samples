% Tomografia de transmissão ultra-sônica
% Doutorado em Engenharia Elétrica
% Universidade Federal da Bahia
% Eduardo Telmo Fonseca Santos
% eduardot@ifba.edu.br
% 13/02/2010

clc
clear
close all

% Velocidade do som (cm/s)
vs = 345 * 100; % 347.97 * 100 @ 28oC ?

data_set = 1
folder_tt = ['C:\Traveltimes\' num2str(data_set) '\'];
tt_sem_bola = load([folder_tt 'tt_sem_bola']); tt_sem_bola = tt_sem_bola.tt * (16.0/12.0e+6);
tt_com_bola = load([folder_tt 'tt_com_bola']); tt_com_bola = tt_com_bola.tt * (16.0/12.0e+6);

figure
imagesc(tt_sem_bola * 1e+6, [250 410])
xlabel('Receptores')
ylabel('Fontes')
colorbar

figure
imagesc(tt_com_bola * 1e+6, [250 410])
xlabel('Receptores')
ylabel('Fontes')
colorbar

ns = 8;
nr = 8;

nx = 8;
nz = nx;

lx = 9.7;
lz = 11;

dx = lx/nx
dz = lz/nz

% Coordenada da primeira fonte/receptor
x0 = 0.5;
z0 = 0.5;

xs = zeros(1,ns) + x0;
zs = ([1:ns] - 1) * ((lz-2*z0)/(ns-1)) + z0

xr = lx * ones(1,nr);
zr = ([1:nr] - 1) * ((lz-2*z0)/(nr-1)) + z0

tt = zeros(ns,nr);
dist = zeros(ns,nr);

% Distancias e tempos de transitos
for i=1:ns,
    for j=1:nr,
        dist(i,j) = sqrt( (xs(i)-xr(j))^2 + (zs(i)-zr(j))^2 );
        tts(i,j) = dist(i,j) / vs;
    end
end

figure
plot(xs, zs, 'r*');
hold on
plot(xr, zr, 'b<');
xlabel('X')
ylabel('Y')
axis('ij');
for i=1:8,
    for j=1:8,        
        plot([xs(i) xr(j)], [zs(i) zr(j)], 'k');
    end
end

for i=1:8,
    text(xs(i) - 0.3, zs(i), num2str(i));   
    text(xr(i) + 0.11, zr(i), num2str(i));   
end

% Calcula matriz tomografica (raio reto)
figure
st = zeros(nx,nz);
g = zeros(ns*nr,nx*nz);
cnt = 1;
for i=1:ns,
    for j=1:nr,
        s = zeros(nx,nz);
        x1 = ceil(max([xs(i)/dx 1]));
        z1 = ceil(max([zs(i)/dz 1]));
        x2 = ceil(max([xr(j)/dx 1]));
        z2 = ceil(max([zr(j)/dz 1]));
        sc = 0;
        len = 1;
        [s,x,z,sc] = bres2d(z1,x1,z2,x2,s,1,sc,len);
        lx = length(x);
        len = dist(i,j) / lx;
        s = s * len;
        [s,x,z,sc] = bres2d(z1,x1,z2,x2,s,1,sc,len);
        imagesc(s, [0 2])
        colorbar
        colormap(gray)
        st = st + s;        
        g(cnt,1:nx*nz) = reshape(s,1,nx*nz);
        cnt = cnt + 1;
        pause(0)
    end
end

m = (1/vs) * ones(nx*nz,1);
d = g * m;

figure
imagesc(st);
axis image
colorbar
colormap(gray)

hold on
plot(xs, zs, 'r*');
hold on
plot(xr, zr, 'b<');

figure
imagesc(g);
axis image
colorbar
colormap(gray)

figure
plot(d)
tv = reshape(tts,ns*nr,1);
hold on
plot(tv, 'r:')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = 40;
vmin = 320;
vmax = 350;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = tt_sem_bola;
dmed = t * (16.0/12.0e+6);
dmed = reshape(dmed,ns*nr,1);

% figure
% plot(dmed)

gtg = g'*g;
I = eye(size(gtg));
mest = pinv(gtg + lambda * I) * g' * dmed;
mest = reshape(mest,nx,nz);
mest = (1./mest) / 100.0;

mest(mest<=vmin) = vmin;
mest(mest>=vmax) = vmax;

figure
imagesc(mest, [vmin vmax])
colormap(gray)
colorbar

% pause(0.75)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xx = 1:8;
yy = 1:8;
f = mest;

figure
[xg,yg] = meshgrid(xx,yy);
surf(xg,yg,real(f))
axis vis3d
shading interp
camlight
lighting phong
material metal
xlabel('X');
ylabel('Y');
zlabel('Velocidade');
colorbar
rotate3d on

% pause(0.75)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = tt_com_bola;
dmed = t;
dmed = reshape(dmed,ns*nr,1);

% figure
% plot(dmed)

gtg = g'*g;
I = eye(size(gtg));
mest = pinv(gtg + lambda * I) * g' * dmed;
mest = reshape(mest,nx,nz);
mest = (1./mest) / 100.0;

mest(mest<=vmin) = vmin;
mest(mest>=vmax) = vmax;

figure
imagesc(mest, [vmin vmax])
colormap(gray)
colorbar

% pause(0.75)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xx = 1:8;
yy = 1:8;
f = mest;

figure
[xg,yg] = meshgrid(xx,yy);
surf(xg,yg,real(f))
axis vis3d
shading interp
camlight
lighting phong
material metal
xlabel('X');
ylabel('Y');
zlabel('Velocidade');
colorbar
rotate3d on

% pause(0.75)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

