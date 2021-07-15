# script to translate the spin number in lattice output files to real space positions for the Fe8 system
BEGIN{
Lx=L
Ly=L
Lz=L
num_in_cell=1
split("1.2676399478585705 0.39071827411333043 -0.07965196320688474",a)
split("-0.44183577847448624 1.274681524880277 -0.4848259797892014",b)
split("-0.25857698036738 -0.09182746019752541 0.9603430499711932",c)
}

FNR > 10 {
    i=int(int(int($1 / num_in_cell) / Lz) / Ly)
    j=int(int(int($1 / num_in_cell) / Lz) % Ly)
    k=int(int($1 / num_in_cell) % Lz)
    x=i*a[1]+j*b[1]+k*c[1]
    y=i*a[2]+j*b[2]+k*c[2]
    z=i*a[3]+j*b[3]+k*c[3]
    print x","y","z","$2
}
