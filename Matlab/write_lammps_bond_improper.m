% write the information into a data file for lammps simulation

function write_lammps_bond_improper(fn,N_i,T_i,M_atom,ty_atom,x_atom,ty_bond, bond,angle,dbond,improper,q,b_s)

N_atom = N_i(1);                % number of atom
N_bond = N_i(2);                % number of bond
N_angle = N_i(3);               % number of angle
N_dihedral = N_i(4);            % number of dihedral
N_improper = N_i(5);            % number of improper

T_atom = T_i(1);                % type of atom
T_bond = T_i(2);                % type of bond
T_angle = T_i(3);               % type of angle
T_dihedral = T_i(4);            % type of dihedral
T_improper = T_i(5);            % type of improper

xlo = b_s(1,1); xhi = b_s(1,2);
ylo = b_s(2,1); yhi = b_s(2,2);
zlo = b_s(3,1); zhi = b_s(3,2);

% open the file with write permission

filename = [fn,'-bond.lam'];                % create file name
fid = fopen(filename, 'w');

% LAMMPS Description

fprintf(fid,'LAMMPS Description\n\n');

% atom, bond, angles, dihedrals, impropers

fprintf(fid,'%d\t', N_atom);
fprintf(fid,'atoms\n');
fprintf(fid,'%d\t', N_bond);
fprintf(fid,'bonds\n');
fprintf(fid,'%d\t', N_angle);
fprintf(fid,'angles\n');
fprintf(fid,'%d\t', N_dihedral);
fprintf(fid,'dihedrals\n');
fprintf(fid,'%d\t', N_improper);
fprintf(fid,'impropers\n\n');

fprintf(fid,'%d\t', T_atom);
fprintf(fid,'atom types\n');
fprintf(fid,'%d\t', T_bond);
fprintf(fid,'bond types\n');
fprintf(fid,'%d\t', T_angle);
fprintf(fid,'angle types\n');
fprintf(fid,'%d\t', T_dihedral);
fprintf(fid,'dihedral types\n');
fprintf(fid,'%d\t', T_improper);
fprintf(fid,'improper types\n\n');

% size of box

fprintf(fid,'%8.3f\t%8.3f\t', xlo,xhi);
fprintf(fid,'xlo xhi\n');
fprintf(fid,'%8.3f\t%8.3f\t', ylo,yhi);
fprintf(fid,'ylo yhi\n');
fprintf(fid,'%8.3f\t%8.3f\t', zlo,zhi);
fprintf(fid,'zlo zhi\n\n');

% Masses

fprintf(fid,'Masses\n\n');
for i = 1:T_atom
    fprintf(fid,'%d %8.3f\n',i,M_atom(i));
end
fprintf(fid,'\n');

% Atoms

fprintf(fid,'Atoms\n\n');
for i = 1:N_atom
    x = x_atom(i,1);
    y = x_atom(i,2);
    z = x_atom(i,3);
    fprintf(fid,'%d %d %d %12.9e %8.3f %8.3f %8.5f\n',i,i,ty_atom(i),q(i),x,y,z);
    
end

% bonds
if N_bond > 0
    fprintf(fid,'\nBonds\n\n');
    
    for i = 1:N_bond
        j = bond(i,1);
        k = bond(i,2);
        fprintf(fid,'%d %d %d %d\n',i,ty_bond(i), j,k);
    end
end

% angles
if N_angle > 0
    fprintf(fid,'\nAngles\n\n');
    
    for i = 1:N_angle-1
        
        fprintf(fid,'%d %d %d %d %d\n',i,angle(i,:));
        
    end
    
    fprintf(fid,'%d %d %d %d %d',i+1,angle(i+1,:));
end

% dihedral angles

if N_dihedral > 0
    fprintf(fid,'\nDihedrals\n\n');
    
    for i = 1:N_dihedral
        fprintf(fid,'%d %d %d %d %d %d\n',i, 1, dbond(i,:));
    end
    
end

if N_improper > 0
    fprintf(fid,'\nImpropers\n\n');
    
    for i = 1:N_improper
        fprintf(fid,'%d %d %d %d %d %d\n',i, improper(i,:));
    end
    
end


fclose(fid);

end
