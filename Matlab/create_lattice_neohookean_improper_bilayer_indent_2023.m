function create_lattice_neohookean_improper_bilayer_indent_2023

clear all
clc

a = 1.;                     % lattice length
ac = num2str(a);

nx = fix(2022);              % number of lattice along x direction
ny = fix(1000);              % number of lattice along y direction

Lx = nx*a;                  % width along x direction
Ly = ny*a;                  % width along y direction

Lxc = num2str(Lx);
Lyc = num2str(Ly);

nh = 4;                     % top layer
nhc = num2str(nh);          % 


nxy = (ny+1)*nx;            % total atom number for peridical boundary condition along x
xy = zeros(nxy,2);          % atomic coordinate


% create the atom coordinate and bond connection in one unit cell
% create angle in one unit cell
ind = 0;
b0 = 1:nx;                % first layer
b1 = [nx+2:2*nx,nx+1];    % second layer
b2 = [2:nx,1];
b3 = [nx+1:2*nx];         % second layer

% create bond the first atom has small number

% consider periodical boundary condition

bond_u1 = zeros(2*nx,2);
bond_u1(:,1) = [b0,b2]';
bond_u1(:,2) = [b1,b3]';
nbon1 = length(bond_u1);
ty_bond_u1 = ones(nbon1,1);

% edge bond

bond_u2 = zeros(2*nx,2);
bond_u2(:,1) = [b0,b0]';
bond_u2(:,2) = [b2,b3]';
nbon2 = length(bond_u2);
ty_bond_u2 = 2*ones(nbon2,1);


% combine two bonds
bond_u = [bond_u1;bond_u2];
ty_bond_u = [ty_bond_u1;ty_bond_u2];
nbon = nbon1 + nbon2;
bond = zeros(ny*nbon+nx,2);
ty_bond = zeros(ny*nbon+nx,1);
indb = 0;

% create improper

improper_u = ones(nx,5);
improper_u(:,2) = [b0]';
improper_u(:,3) = [b1]';
improper_u(:,4) = [b3]';
improper_u(:,5) = [b2]';
nimp = nx;

improper = zeros(ny*nimp,5);
indi = 0;

% coordinate, bond and angle, improper angle
% for substrate

for i = 1:ny-nh
    
    xy(ind+1:ind+nx,1) = a*[1:nx];
    xy(ind+1:ind+nx,2) = (i-1)*a;
    ind = ind+nx;


    bond(indb+1:indb+nbon,:) = bond_u + nx*(i-1);
    ty_bond(indb+1:indb+nbon) = ty_bond_u;
    indb = indb+nbon;
       
    
    improper(indi+1:indi+nimp,2:5) =improper_u(:,2:5) + nx*(i-1);
    improper(indi+1:indi+nimp,1) =improper_u(:,1);
    indi = indi + nimp;
    
end

% for film

ty_bond_u = ty_bond_u + 2;
improper_u(:,1) = improper_u(:,1) + 1;

for i = ny-nh+1:ny
    
    xy(ind+1:ind+nx,1) = a*[1:nx];
    xy(ind+1:ind+nx,2) = (i-1)*a;
    ind = ind+nx;


    bond(indb+1:indb+nbon,:) = bond_u + nx*(i-1);
    ty_bond(indb+1:indb+nbon) = ty_bond_u;
    indb = indb+nbon;
    
    
    improper(indi+1:indi+nimp,2:5) = improper_u(:,2:5) + nx*(i-1);
    improper(indi+1:indi+nimp,1) = improper_u(:,1);
    indi = indi + nimp;
    
end

% add the last layer of atom

xy(ind+1:ind+nx,1) = a*[1:nx];
xy(ind+1:ind+nx,2) = ny*a;

% create atom type

ty_atom = ones(nxy,1);

film_thickness = (ny-nh)*a - 1e-5;

la = xy(:,2) > film_thickness;

ty_atom(la) = 2*ones(sum(la),1);

% transform the lattice to orignal center

xc = nx*a/2 + a;
yc = ny*a/2;
xy(:,1) = xy(:,1)-xc;
xy(:,2) = xy(:,2) - yc;

display(yc)

% add the last bond (film top surface)
bond(indb+1:indb+nx,1) = ind+1:ind+nx;
bond(indb+1:indb+nx,2) = [ind+2:ind+nx,ind+1];
ty_bond(indb+1:indb+nx) = 4*ones(nx,1);

% find top and bottom edge bond and change its type to 3
n1 = bond(:,1);
n2 = bond(:,2);

yn = xy(n1,2) + xy(n2,2);

% bottom 
l = yn < -ny*a+1e-3;
ty_bond(l) = 5;

% film-substrate interface

l = abs(yn - (ny-2*nh)*a) < 1e-3;
ty_bond(l) = 6;

% top
l = yn > ny*a-1e-3;
ty_bond(l) = 7;


R = 2000;              % indenter radius
dd = 6*a;

Rc = num2str(R);
ddc = num2str(dd);

y0 = yc + dd + R;     % height of the indenter 

ntheta = fix(pi/(a/R));    % number of segment
theta = [0:2*ntheta-1]/2/ntheta;

x_r = R*sin(theta*2*pi);
y_r = R*cos(theta*2*pi) + y0;

l = y_r > y0-R+50;

x_r(l) = [];
y_r(l) = [];

xy_r = [x_r;y_r]';
ty_atom_r = 3*ones(length(xy_r),1);

% cut the top half circle

xy = [xy;xy_r];
ty_atom = [ty_atom;ty_atom_r];


% create atom type

angle = [];
dbond = [];                     % no dihedral angle

fn = ['lattice-bilayer-indent-improper-a',ac,'-R',Rc,'-d',ddc,'-Lx',Lxc,'-Ly',Lyc,'-hc',nhc];

N_i = [length(xy),length(bond),length(angle),length(dbond),length(improper)];

T_i = [3,7,0,0,2];


M_atom = [2,50,1];
x_atom = [xy,zeros(N_i(1),1)];

% charge for additional variable later

q = zeros(N_i(1),1);

b_s = [-nx*a/2 nx*a/2; -ny*a/2*1.01 ny*a/2*1.35; -50. 50.];       % box

write_lammps_bond_improper(fn,N_i,T_i,M_atom,ty_atom,x_atom,ty_bond,bond,angle,dbond,improper,q,b_s);


end