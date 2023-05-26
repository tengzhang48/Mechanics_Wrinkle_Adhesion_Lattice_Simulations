function create_lattice_neohookean_indent_improper_2023

clear all
clc

a = 1.;                     % lattice length
ac = num2str(a);

nx = fix(1000);              % number of lattice along x direction
ny = fix(800);              % number of lattice along y direction

Lx = nx*a;                  % width along x direction
Ly = ny*a;                  % width along y direction

Lxc = num2str(Lx);
Lyc = num2str(Ly);

nxy = (ny+1)*(nx+1);        % total atom number for non-peridical boundary condition
xy = zeros(nxy,2);          % atomic coordinate


% create the atom coordinate and bond connection in one unit cell
% create angle in one unit cell
ind = 0;
b0 = 1:nx;                % first layer
b1 = nx+3:2*nx+2;           % second layer

% create bond the first atom has small number

bond_u1 = zeros(2*nx,2);
bond_u1(:,1) = [b0,b0+1]';
bond_u1(:,2) = [b1,b1-1]';
nbon1 = length(bond_u1);
ty_bond_u1 = ones(nbon1,1);

% edge bond
bond_u2 = zeros(2*nx+1,2);
bond_u2(:,1) = [b0,b0,nx+1,]';
bond_u2(:,2) = [b0+1,b1-1,2*nx+2]';
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
improper_u(:,4) = [b1-1]';
improper_u(:,5) = [b0+1]';
nimp = nx;

improper = zeros(ny*nimp,5);
indi = 0;

% coordinate, bond and angle
for i = 1:ny
    
    xy(ind+1:ind+nx+1,1) = a*[1:nx+1];
    xy(ind+1:ind+nx+1,2) = (i-1)*a;
    ind = ind+nx+1;

    bond(indb+1:indb+nbon,:) = bond_u + (nx+1)*(i-1);
    ty_bond(indb+1:indb+nbon) = ty_bond_u;
    indb = indb+nbon;
    
    
    improper(indi+1:indi+nimp,2:5) =improper_u(:,2:5) + (nx+1)*(i-1);
    improper(indi+1:indi+nimp,1) =improper_u(:,1);
    indi = indi + nimp;
    
end

% add the last layer of atom

xy(ind+1:ind+nx+1) = a*[1:nx+1];
xy(ind+1:ind+nx+1,2) = ny*a;

% transform the lattice to orignal center

xc = nx*a/2 + a;
yc = ny*a/2;
xy(:,1) = xy(:,1)-xc;
xy(:,2) = xy(:,2) - yc;

display(yc)

% add the last bond
bond(indb+1:indb+nx,1) = ind+1:ind+nx;
bond(indb+1:indb+nx,2) = ind+2:ind+nx+1;
ty_bond(indb+1:indb+nx) = 2*ones(nx,1);

% find edge bond and change its type to 3
n1 = bond(:,1);
n2 = bond(:,2);

xn = xy(n1,1) + xy(n2,1);
l = xn > nx*a-1e-3;
ty_bond(l) = 3;
l = xn < -nx*a+1e-3;
ty_bond(l) = 3;

yn = xy(n1,2) + xy(n2,2);
l = yn > ny*a-1e-3;
ty_bond(l) = 3;
l = yn < -ny*a+1e-3;
ty_bond(l) = 3;

% create atom type

ty_atom = ones(nxy,1);

% Create indenter

% R = 40;              % indenter radius
% y0 = yc + 8*a + R;     % height of the indenter 

R = 200;              % indenter radius
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
ty_atom_r = 2*ones(length(xy_r),1);

% cut the top half circle


xy = [xy;xy_r];


ty_atom = [ty_atom;ty_atom_r];

% angle = [ty_angle,angle];

angle = [];

dbond = [];                     % no dihedral angle

fn = ['lattice-indent-improper-a',ac,'-R',Rc,'-d',ddc,'-Lx',Lxc,'-Ly',Lyc];

N_i = [length(xy),length(bond),length(angle),length(dbond),length(improper)];
T_i = [2,3,0,0,1];
M_atom = [1, 1];
x_atom = [xy,zeros(N_i(1),1)];

% charge for additional variable later

q = zeros(N_i(1),1);

b_s = [-(nx+5)*a/2 (nx+5)*a/2; -(ny+5)*a/2 ny*a/2 + 80*a; -50. 50.];       % box

write_lammps_bond_improper(fn,N_i,T_i,M_atom,ty_atom,x_atom,ty_bond,bond,angle,dbond,improper,q,b_s);


end
