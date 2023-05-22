
# package containing many of the scripts I've written - for ease of use

# list of dependencies:
# ase; numpy; scipy; matplotlib; spglib; prettytable; pymatgen; 

##################################################################################
# list of functions: 
# get_irr_kpts(atoms,kpts,is_shift=[0,0,0])
# get_line(x,y,extra=0.1)
# greplines(cmd)
# chk_output(cmd)
# get_wf_environ(path)
# get_wf_implicit(path)
# get_wf_explicit(path)
# get_chgcar(path,locpot=False,pavg=False)
# pos_swap(atoms,ind1,ind2)
# scale_metal(atoms,a_old,a_new,old_sym,new_sym=None)
# get_dos(path)
# get_parab used by the fed function
# fed(path, label, fig)
# get_n0(path,use_pbe = True):
# get_omega(path):
# match_cell(ref_atoms,change_atoms,lower_vac,anchor_atom=None)
# fmax(atoms):
# set_pot(atoms,calc,desired_U,tolerance=0.02):
# get_closest(ref,atoms,ind):
# reindex_atoms(ref_atoms,reindex_atoms,manual_skip_atoms=[]):
# const_U_relax(atoms,calc,desired_U,tolerance=0.02,ediffg=0.05):
# const_U_dimer(atoms,calc,desired_U,ediffg=0.05):
# const_U_FBL(atoms,calc,desired_U,ind1,ind2,z_cutoff=None,ediffg=0.05):
##################################################################################



##################################################################################
# Here are a few module wide constants that get used in various functions. You can
# update them with e.g. "common._she_U = 4.6" or "common._tolerance_U = 0.01".

# SHE reference potential
_she_U=4.43

# tolerance criteria used when optimizing NELECT for a given potential, in Volts
_tolerance_U=0.02

# Planck's constant in units of eV s
_h = 4.135667696e-15

#Boltzmann's constant in units of eV
_kb = 8.617333262e-5
##################################################################################


def get_irr_kpts(atoms,kpts,is_shift=[0,0,0]):
    # returns the number of irreducible kpoints
    # given an atoms object and a desired k-pt grid

    import numpy as np
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga
    from pymatgen import Structure
    import os
    atoms.write('./tmp.cif')
    mesh = sga(structure=Structure.from_file('./tmp.cif')).get_ir_reciprocal_mesh(mesh=kpts)
    os.system('rm ./tmp.cif')
    return len(mesh)

    # kpts = np.array(kpts)
    # lattice = np.array([atoms.cell[0],atoms.cell[1],atoms.cell[2]])
    # positions = atoms.get_positions()
    # numbers = [1,]*len(positions)
    # cell = (lattice,positions,numbers)

    # mapping,grid = spglib.get_ir_reciprocal_mesh(kpts, cell, is_shift=is_shift)
    # return len(np.unique(mapping))


def get_line(x,y,extra=0.1,extramin=0.0,extraplus=0.0,return_mae=False):
    import sys,subprocess
    import numpy as np
    # returns a np array for the x and y axis of a line
    # given some data
    xax = np.linspace(min(x)-extra-extramin,max(x)+extra+extraplus,10)
    a,b = np.polyfit(x,y,1)
    yax = a*xax+b
    if not return_mae:
        return xax,yax,a,b
    mae = 0
    for i in range(len(x)):
        mae += abs(y[i]-(a*x[i]+b))
    mae = mae/len(x)
    return xax,yax,a,b,mae


def greplines(cmd):
    import sys,subprocess
    # easier subprocess use - auto split by newline character
    ver = sys.version_info[0]
    if ver == 2:
        try:
            return subprocess.check_output(cmd,shell=True).split('\n')[:-1]
        except subprocess.CalledProcessError: 
            return ''
    if ver == 3:
        try:
            return subprocess.check_output(cmd,shell=True).decode('utf-8').split('\n')[:-1]
        except subprocess.CalledProcessError: 
            return ''

def chk_output(cmd):
    import sys,subprocess
    # easier subprocess use
    ver = sys.version_info[0]
    if ver == 2:
        try:
            return subprocess.check_output(cmd,shell=True)
        except subprocess.CalledProcessError: 
            return ''
    if ver == 3:
        try:
            return subprocess.check_output(cmd,shell=True).decode('utf-8')
        except subprocess.CalledProcessError: 
            return ''

def param_set(param,val):
    # read in INCAR file, set parameter to desired level
    incar = open('INCAR','r')
    lines = incar.readlines()
    incar.close()
    check = 0
    for i,line in enumerate(lines):
        if param+' =' in line:
            del lines[i]
            check = 1
            if val == 'del':
                break
            lines.insert(i,' ' + param + ' = ' + str(val) + '\n')

    # if the parameter is not in the INCAR file, add it to the bottom
    if check == 0 and val != 'del':
        lines.insert(len(lines),' ' + param + ' = ' + str(val) + '\n')

    incar = open('INCAR','w')
    incar.writelines(lines)
    incar.close()

def get_wf_environ(path):
    import os,sys,subprocess
    import numpy as np
    from ase.io import read
    # check if environ was used
    try:
    # grep log file for fermi level
        test = subprocess.check_output("grep 'Environ Module' "+path+"/log",shell=True)
        out = subprocess.check_output("grep 'the Fermi energy is' "+path+"/log | tail -n 1",shell=True).split('\n')[-2]
        fermi = float(out.split()[-2])

        out = subprocess.check_output("grep ' due to the parabolic pbc-correction' "+path+"/log | tail -n 1",shell=True).split('\n')[-2]
        shift = float(out.split()[-2])

        return -1*(fermi+shift)
    except:
        if os.path.isfile(path+'/wf.out'):
            f = open(path+'/wf.out','r')
            lines = f.readlines()
            f.close()
            return float(lines[0].rstrip())
        else:
            pot = pickle.load(open(path+'/elpot.pkl','rb'))
# load pickle file and take only the potential data
            pot = pot[b'data']

# get planar average by taking mean across x and y axes
            pavg = np.mean(np.mean(pot,axis=0),axis=0)

# assumes the system is positioned with (roughly) atoms in the center of the cell
# with vacuum above and below
            lower_bound = int(3*len(pavg)/4)
            upper_bound = int(9*len(pavg)/10)
            vac = np.mean(pavg[lower_bound:upper_bound])

# grep log file for fermi level
# using absolutely ridiculous python 3 subprocess command
            out = subprocess.check_output("grep 'the Fermi energy is' "+path+"/log | tail -n 1",shell=True).split('\n')[-2]
            fermi = float(out.split()[-2])

            f = open(path+'/wf.out','w')
            f.writelines('%6f\n'%(vac-fermi))
            f.close()
            return vac-fermi

def get_wf_implicit(path):
    import os,sys,subprocess
    import numpy as np
    from ase.io import read
    if not os.path.exists('%s/OUTCAR'%path):
        print('No OUTCAR found -- use vasprun.xml instead')
        fermi = float(greplines('grep fermi %s/vasprun.xml'%path)[0].split()[-2])
    else:
        out1 = greplines('grep fermi '+path+'/OUTCAR | tail -n 1')
        fermi = float(out1[0].split()[2])

    try:
        out2 = greplines('grep FERMI_SHIFT '+path+'/opt.log | tail -n 1')
        shift = float(out2[0].split(' = ')[-1])
    except:
        out2 = greplines('grep FERMI_SHIFT '+path+'/vasp.out | tail -n 1')
        shift = float(out2[0].split(' = ')[-1])

    return -1*(fermi+shift)

def get_wf_explicit(path):
    if os.path.isfile(path+'/wf.out'):
        f = open(path+'/wf.out')
        lines = f.readlines()
        return float(lines[0].rstrip())

def get_chgcar(path,locpot=False,pavg=False):
    # get CHGCAR data as a numpy array
    # optionally, return the planar average of this density
    # if the file is named "LOCPOT", the flag 'locpot' will be 
    # set to True
    if os.path.isfile(path[:-6]+'pavg.txt') and pavg:
        # already did the analysis, just read and return the file
        f = open(path[:-6]+'pavg.txt')
        lines = f.readlines()
        f.close()
        pavg = [float(line.rstrip()) for line in lines]
        return np.array(pavg)

    f = open(path,'r')
    lines = f.readlines()
    f.close()

    # get unit cell volume
    a = np.array([float(lines[2].split()[i]) for i in range(3)])
    b = np.array([float(lines[3].split()[i]) for i in range(3)])
    c = np.array([float(lines[4].split()[i]) for i in range(3)])
    V = np.dot(np.cross(b,c),a)

    # get header length and find nx, ny, nz
    for i,line in enumerate(lines):
        if len(line.strip()) == 0:
            break
    i += 1
    line = lines[i]
    vox_n = np.array([int(line.split()[j]) for j in range(3)])
    n_tot = vox_n[0]*vox_n[1]*vox_n[2]
    i += 1

    # find footer length
    out = greplines('grep -n augmentation ' + path).split('\n')[0]
    k = len(lines)-int(out.split(':')[0])

    density = np.genfromtxt(path,skip_header=i,skip_footer=k,invalid_raise=False)
    density = np.reshape(density,density.size)
    if density.size < n_tot:
        for chg in lines[-k-1].split():
            density = np.append(density,float(chg))
    density = np.reshape(density,vox_n,order='F')
    if not locpot:
        density /= n_tot
    assert (density.shape == vox_n).all()

    if pavg:
        pavg = np.mean(np.mean(density,axis=0),axis=0)
        f = open(path[:-6]+'pavg.txt','w')
        lines = [str(i)+'\n' for i in pavg]
        f.writelines(lines)
        f.close()
        return pavg
    return density

def pos_swap(atoms,ind1,ind2):
    # swap position of two atoms
    x1,y1,z1 = [atoms[ind1].x,atoms[ind1].y,atoms[ind1].z]
    x2,y2,z2 = [atoms[ind2].x,atoms[ind2].y,atoms[ind2].z]

    atoms[ind1].x = x2
    atoms[ind1].y = y2
    atoms[ind1].z = z2

    atoms[ind2].x = x1
    atoms[ind2].y = y1
    atoms[ind2].z = z1

    return atoms

def scale_metal(atoms,a_old,a_new,old_sym,new_sym=None):
    # useful for changing lattice constant of metal, or 
    # changing the metal being used in a trajectory

    # generally will only work for othorhombic cells
    ratio = a_new/a_old

    if new_sym is not None:
        for atom in atoms:
            if atom.symbol == old_sym:
                atom.symbol = new_sym
    
    atoms.set_cell(atoms.cell*ratio)
    translate = [(1-ratio)*atoms.cell[i][i] for i in range(3)]
    for atom in atoms:
        # if atom.symbol == old_sym or atom.symbol == new_sym:
        atom.x *= ratio
        atom.y *= ratio
        atom.z *= ratio
        # else:
            # atom.x -= translate[0]
            # atom.y -= translate[1]
            # atom.z -= translate[2]
    return atoms

# now for a few functions that get_dos depends on
# taken from VTST, split_dos.py from Henkelman group
def read_dosfile():
    f = open("DOSCAR", 'r')
    lines = f.readlines()
    f.close()
    index = 0
    natoms = int(lines[index].strip().split()[0])
    index = 5
    nedos = int(lines[index].strip().split()[2])
    efermi = float(lines[index].strip().split()[3])
    # print natoms, nedos, efermi

    return lines, index, natoms, nedos, efermi

###READ POSCAR or CONTCAR and save pos
def read_posfile():
    from ase.io import read

    try:
        atoms = read('POSCAR')
    except IOError:
        print("[__main__]: Couldn't open input file POSCAR, atomic positions will not be written...\n")
        atoms = []
   
    return atoms

### WRITE DOS0 CONTAINING TOTAL DOS ###
def write_dos0(lines, index, nedos, efermi):
    
    fdos = open("DOS0", 'w')
    index +=1
    line = lines[index+2].strip().split()
    ncols = int(len(line))
    fdos.write('# %d \n' % (ncols))

    for n in range(1,nedos):
        index +=1
        e = float(lines[index].strip().split()[0])
        e_f = e-efermi
        fdos.write('%15.8f ' % (e_f))
        
        for col in range(1, ncols):
            dos = float(lines[index].strip().split()[col])
            fdos.write('%15.8f ' % (dos))
            col +=1  
        fdos.write('\n ')          
        n +=1  
    return index

### LOOP OVER SETS OF DOS, NATOMS ###
def write_nospin(lines, index, nedos, natoms, ncols, efermi):
    import numpy as np
    
    atoms = read_posfile()
    if len(atoms) < natoms:
        pos = np.zeros((natoms, 3))
    else:
        pos = atoms.get_positions()

    for i in range(1,natoms+1):
        si = str(i)
        
    ## OPEN DOSi FOR WRITING ##
        fdos = open("DOS"+si, 'w')
        index +=2
        ia = i -1
        fdos.write('# %d \n' % (ncols))
        fdos.write('# %15.8f %15.8f %15.8f \n' % (pos[ia,0], pos[ia,1], pos[ia,2]))
    
    ### LOOP OVER NEDOS ###
        for n in range(1,nedos):
            index +=1
            e = float(lines[index].strip().split()[0])
            e_f = e-efermi
            fdos.write('%15.8f ' % (e_f))
            
            for col in range(1, ncols):
                dos = float(lines[index].strip().split()[col])
                fdos.write('%15.8f ' % (dos))
                col +=1
            fdos.write('\n ')              
            n +=1
        i+=1
    fdos.close()

def write_spin(lines, index, nedos, natoms, ncols, efermi):
    import numpy as np
    #pos=[]
    atoms = read_posfile()
    if len(atoms) < natoms:
        pos = np.zeros((natoms, 3))
    else:
        pos = atoms.get_positions()

    nsites = (ncols -1)/2
    
    for i in range(1,natoms+1):
        si = str(i)
    ## OPEN DOSi FOR WRITING ##
        fdos = open("DOS"+si, 'w')
        index +=2
        ia = i-1
        fdos.write('# %d \n' % (ncols))
        fdos.write('# %15.8f %15.8f %15.8f \n' % (pos[ia,0], pos[ia,1], pos[ia,2]))
    
    ### LOOP OVER NEDOS ###
        for n in range(1,nedos):
            index +=1   
            e = float(lines[index].strip().split()[0])
            e_f = e-efermi
            fdos.write('%15.8f ' % (e_f))
            
            for col in range(1, nsites):
                dos_up = float(lines[index].strip().split()[col])
                dos_down = float(lines[index].strip().split()[col+1])*-1
                fdos.write('%15.8f %15.8f ' % (dos_up, dos_down))
                col +=1
            fdos.write('\n ')
            n +=1
        i+=1
        fdos.close()

# convert DOSCAR into useable pickle files
# bottom portion written by LDC
def get_dos(path='./'):
    import numpy as np
    from ase.io import read
    import sys,os,pickle
    top = os.getcwd()
    os.chdir(path)
    lines, index, natoms, nedos, efermi = read_dosfile()
    index = write_dos0(lines, index, nedos, efermi)

    ##Test if there a spin calculation was performed ##
    line = lines[index+2].strip().split()
    ncols = int(len(line)) 
    if ncols==7 or ncols==19 or ncols==9 or ncols==33:
        write_spin(lines, index, nedos, natoms, ncols, efermi)
        is_spin=True
    else: 
        write_nospin(lines, index, nedos, natoms, ncols, efermi)
        is_spin=False
    # print("Spin unrestricted calculation: ", is_spin)

    atoms = read('POSCAR')
    dosatoms = range(len(atoms))
    content = [0]*len(dosatoms)

    for i in range(len(dosatoms)):
        fname = 'DOS'+str(dosatoms[i]+1)
        with open(fname) as f:
            content[i] = f.readlines()[2:-1]  #read in energy and dos
        f.close()

    o2s = [0]*len(dosatoms) #number of oxygens in total; easy enough
    for i in range(len(o2s)):  
        o2s[i] = [0]*10  #10 different channels: energy + 1s + 3p + 5d orbitals
        for j in range(len(o2s[i])):
            o2s[i][j] = []  
            for k in range(len(content[i])):
                o2s[i][j].append(float(content[i][k].split()[j]))  #ith oxygen, jth channel, kth line in content

    o2s = np.array(o2s)  #change the array to an numpy array so we can do fun stuff

    total_dos = [0]*len(atoms)

    for i in range(len(o2s)):
        total_dos[i] = [0]*len(o2s[0][0])
        for j in range(1,len(o2s[0])):
            for k in range(len(o2s[0][0])):
                total_dos[i][k] += o2s[i][j][k]

    f = open('dos.pickle', 'wb')
    pickle.dump(total_dos, f)
    f.close()

    f = open('energy.pickle', 'wb')
    pickle.dump(o2s[0][0], f)
    f.close()

    for i in range(len(dosatoms)):
        os.system('rm DOS'+str(dosatoms[i]+1))
    os.system('rm DOS0')
    print('pDOS: %s\nEnergy: %s'%(path+'dos.pickle',path+'energy.pickle'))
    os.chdir(top)


# get_parab function is used by the next function, fed
def get_parab(x0,x1,y0,y1,side):
    import numpy as np
    if side == 0:
        # left side of parabola
        a = (y0-y1)/(x0**2-x1**2+2*x1*(x1-x0))
        b = -2*a*x1
        c = y0-a*x0**2+2*a*x1*x0
    elif side == 1:
        # right side of parabola
        a = (y0-y1)/(x0**2-x1**2-2*x0*(x0-x1))
        b = -2*a*x0
        c = y0-a*x0**2+2*a*x0**2
    xax = np.linspace(x0,x1,50)
    yax = a*np.square(xax)+b*xax+c
    return xax,yax


# given a figure object, a reaction path, and a label,
# returns a free energy diagram (FED) of the reaction path.
###
# reaction path should be a list of tuples, with the first
# element of the tuple being either 'min' or 'ts' for the 
# respective types (energy minimum or transition state).
# The second element of the tuple is the energy relative to
# the first entry. 
# e.g., fig = plt.figure(figsize=(6*1.618,6))
# label = 'CO + COH -> OCCOH'
# path = [('min',0),('ts',1.5),('min',1.0)]
def fed(path,label,fig):
    import matplotlib.pyplot as plt
    import numpy as np
    from common import get_line

    es = [i[0] for i in path]
    types = [i[1] for i in path] # 'min' for minimum, 'ts' for barrier

    p1 = plt.plot(0,es[0],'-',label=label)
    xtrack = 0
    for i in range(len(path)):
        if i == len(path)-1:
            # end of FED
            plt.plot([xtrack,xtrack+1],[es[i]]*2,'-',color=p1[0].get_color())
            continue
        if types[i] == 'min':
            # current energy is a local minimum
            plt.plot([xtrack,xtrack+1],[es[i]]*2,'-',color=p1[0].get_color())
            xtrack += 1
            if types[i+1] == 'min':
                # next energy is also a min, connect with a dashed line
                xax,yax,a,b = get_line([xtrack,xtrack+2],[es[i],es[i+1]],extra=0.0)
                plt.plot(xax,yax,'--',color=p1[0].get_color())
                xtrack += 2
            # otherwise do nothing
        elif types[i] == 'ts':
            # current energy is a barrier, so last one was a min.
            # need to connect with a parabola
            # known info: point at left edge (connect to min)
            # point at top (barrier)
            # slope at top (= 0)
            y0 = es[i-1]; y1 = float(es[i])
            x0 = xtrack; x1 = xtrack+1
            xax,yax = get_parab(x0,x1,y0,y1,side=0)
            # first half of parab
            plt.plot(xax,yax,'-',color=p1[0].get_color())
            xtrack += 1
            
            #second half
            y0 = es[i]; y1 = es[i+1]
            x0 = xtrack; x1 = xtrack+1
            xax,yax = get_parab(x0,x1,y0,y1,side=1)
            plt.plot(xax,yax,'-',color=p1[0].get_color())
            xtrack += 1

    plt.xticks([])

def get_zval(symbol,path,use_pbe):
    pbe_dict = {'Pd': 10.0, 'Sb': 5.0, 'Cr': 12.0, 'Se': 6.0, 'Sn': 14.0, 'Li': 3.0, 'He': 2.0, 'Fr': 9.0, 'Cs': 9.0, 'Nd': 11.0, 'Ac': 11.0, 'Ho': 9.0, 'Eu': 8.0, 'Ni': 10.0, 'Po': 16.0, 'Am': 17.0, 'Be': 2.0, 'Sr': 10.0, 'Al': 3.0, 'Mg': 2.0, 'Ir': 9.0, 'Ge': 14.0, 'Sm': 11.0, 'Cu': 11.0, 'Ra': 10.0, 'Hf': 10.0, 'Co': 9.0, 'Fe': 8.0, 'Ga': 13.0, 'Ba': 10.0, 'Te': 6.0, 'U': 14.0, 'Tb': 9.0, 'I': 7.0, 'Er': 9.0, 'N': 5.0, 'Rn': 8.0, 'Ca': 10.0, 'Nb': 13.0, 'S': 6.0, 'Tl': 13.0, 'F': 7.0, 'O': 6.0, 'Ta': 11.0, 'Pb': 14.0, 'H': 1.0, 'Zn': 12.0, 'Na': 7.0, 'Pu': 16.0, 'Gd': 9.0, 'Pt': 10.0, 'Sc': 11.0, 'V': 13.0, 'Lu': 9.0, 'Dy': 9.0, 'Pa': 13.0, 'Si': 4.0, 'Ag': 11.0, 'Kr': 8.0, 'Pm': 11.0, 'Tc': 13.0, 'Ar': 8.0, 'Rb': 9.0, 'Au': 11.0, 'W': 12.0, 'Ne': 8.0, 'At': 7.0, 'Np': 15.0, 'Tm': 9.0, 'As': 5.0, 'Hg': 12.0, 'K': 9.0, 'Br': 7.0, 'Os': 8.0, 'Yb': 8.0, 'Cd': 12.0, 'Cm': 18.0, 'Pr': 11.0, 'Ru': 14.0, 'Mo': 14.0, 'In': 13.0, 'Cl': 7.0, 'La': 11.0, 'Ce': 12.0, 'C': 4.0, 'Th': 12.0, 'B': 3.0, 'Y': 11.0, 'Mn': 13.0, 'Bi': 15.0, 'Re': 7.0, 'Xe': 8.0, 'Rh': 9.0, 'Zr': 12.0, 'Ti': 12.0, 'P': 5.0}
    if use_pbe:
        return pbe_dict[symbol]
    # else:
        # out = greplines('grep ZVAL %s/POTCAR'%path)[0]
        # return float(out.split()[5])


def get_n0(path,use_pbe = True,atoms_file=None):
    from ase.io import read
    import os
    # if use_pbe:
        # print('assuming standard PBE POTCAR')
    files = ['CONTCAR','XDATCAR','OUTCAR','vasprun.xml','POSCAR']
    geometry_present = False
    if atoms_file is None:
        for file in files:
            if os.path.exists('%s/%s'%(path,file)):
                geometry_present=True
                break
        atoms_file = file
        if not geometry_present:
            print('no geometry file present, abort')
            exit()
    n0 = 0
    atoms = read('%s/%s'%(path,atoms_file))
    for atom in atoms:
        n0 += get_zval(atom.symbol,path,use_pbe)
    return n0

def get_omega(path):
    from ase.io import read
    import os
    # calculate the grand canonical energy given a path
    # run save_space.py first
    n0 = get_n0(path)

    if not os.path.exists('%s/OUTCAR'%path):
        # print('no OUTCAR in directory, trying vasprun.xml instead')
        e = read('%s/vasprun.xml'%path).get_potential_energy()
        nel = float(greplines('grep NELECT %s/vasprun.xml'%path)[0].split()[-1][:-4])
        fermi = float(greplines('grep fermi %s/vasprun.xml'%path)[0].split()[-2])
    else:
        e = read('%s/OUTCAR'%path).get_potential_energy()
        nel = float(greplines('grep NELECT %s/OUTCAR'%path)[0].split()[2])
        out1 = greplines('grep fermi %s/OUTCAR | tail -n 1'%path)
        fermi = float(out1[0].split()[2])
    q = nel-n0
    # print(path,q,fermi)

    return e-q*fermi

def match_cell(ref_atoms,change_atoms,lower_vac,anchor_atom=None):
    from ase.io import read
    new_cell = ref_atoms.cell
    change_atoms.cell = new_cell
    if anchor_atom is not None:
        # rather than change positions according to a vacuum,
        # reference everything to the z-position of an atom
        # in the reference atoms object
        zdiff = change_atoms[anchor_atom].z - ref_atoms[anchor_atom].z
        for atom in change_atoms:
            atom.z -= zdiff
        return change_atoms

    # otherwise, use the lower vacuum constraint to reposition
    zs = [atom.z for atom in change_atoms]
    zs.sort()
    zdiff = zs[0] - lower_vac
    for atom in change_atoms:
        atom.z -= zdiff
    return change_atoms

def fmax(atoms):
    import os,sys,pickle,math
    from ase.io import read
    # given an atoms object, return the maximum force on
    # an unconstrained atom
    if atoms.constraints != []:
        unconstrained = [atom.index for atom in atoms if atom.index not in list(atoms.constraints[0].index)]
    else:
        unconstrained = [atom.index for atom in atoms]
    ftemp = []
    for ind in unconstrained:
        f = atoms.get_forces()[ind]
        ftemp.append((f[0]**2+f[1]**2+f[2]**2)**0.5)
    ftemp.sort()
    return ftemp[-1]

def set_pot(atoms,calc,desired_U):
    import os,sys,pickle,math
    from ase.io import read
    # determine NELECT required to have potential=desired_U
    calc.bool_params['lcharg'] = False
    calc.int_params['ichain'] = 0
    calc.bool_params['lwave'] = True
    calc.int_params['nsw'] = 0
    calc.exp_params['ediff'] = 1.0e-4
    calc.int_params['ibrion'] = 1

    # previous optimization was done, use that as starting point
    if os.path.isfile('nelect_data.pkl') and os.stat('nelect_data.pkl').st_size != 0:
        nel_data = pickle.load(open('./nelect_data.pkl','rb'))
        calc.float_params['nelect'] = nel_data['nelect'][-1]
        atoms.set_calculator(calc)
    else:
        nel_data = {}
        nel_data['nelect'] = []
        nel_data['potential'] = []
        nel_data['energy'] = []
        print('Running the first single point to get PZC')
        atoms.get_potential_energy()
        
        # store info from the first single point
        nel_data['potential'].append(get_wf_implicit('./')-_she_U)
        nel_data['energy'].append(atoms.get_potential_energy())
        nel_out = float(greplines('grep NELECT OUTCAR')[0].split()[2])
        nel_data['nelect'].append(nel_out)
        pickle.dump(nel_data,open('nelect_data.pkl','wb'))

    # no need to run further optimization if you're already at the desired potential
    if abs(nel_data['potential'][-1]-desired_U) < _tolerance_U:
        # if len(nel_data['potential']) > 2:
            # if abs(nel_data['potential'][-1]-nel_data['potential'][-2]) < 0.001:
                # # strange bug -- gets stuck, so update nelect by a tiny amount
                # # to force a new single point calculation
                # calc.float_params['nelect'] += 1e-3
                # return
        return

    if len(nel_data['nelect']) < 2:
        # only one data point - do another single point with slightly more
        # electrons to get an initial gradient for newton's method
        # initial guess for C: 1 e/V
        calc.float_params['nelect'] = nel_data['nelect'][-1]+(nel_data['potential'][-1]-desired_U)
        atoms.set_calculator(calc)
        atoms.get_potential_energy()

        nel_data['potential'].append(get_wf_implicit('./')-_she_U)
        nel_data['energy'].append(atoms.get_potential_energy())
        nel_out = float(greplines('grep NELECT OUTCAR')[0].split()[2])
        nel_data['nelect'].append(nel_out)
        pickle.dump(nel_data,open('nelect_data.pkl','wb'))

    #start the optimization, initialize vars
    while abs(nel_data['potential'][-1]-desired_U) > _tolerance_U:
        # Newton's method to optimize NELECT
        grad_numer = nel_data['potential'][-2]-nel_data['potential'][-1]
        grad_denom = nel_data['nelect'][-2]-nel_data['nelect'][-1]
        if abs(grad_denom) < 0.0001:
            diff = 0.001
        else:
            grad = grad_numer/grad_denom
            y = nel_data['potential'][-1]-desired_U
            diff = abs(y)**2/(y*grad)

        # don't take too big of a step ..
        # can happen if two subsequent steps are too close together
        if diff > 5.0:
            diff = 0.75
        elif diff < -5.0:
            diff = -0.75

        # update nelect
        new_nel = nel_data['nelect'][-1] - diff

        #check if nelect is nan
        if math.isnan(new_nel):
            print('Error: Check NELECT (nan)')
            sys.exit()

        # check guess from newton's method
        calc.float_params['nelect'] = new_nel
        atoms.set_calculator(calc)
        atoms.get_potential_energy()

        nel_data['potential'].append(get_wf_implicit('./')-_she_U)
        nel_data['energy'].append(atoms.get_potential_energy())
        nel_out = float(greplines('grep NELECT OUTCAR')[0].split()[2])
        nel_data['nelect'].append(nel_out)
        pickle.dump(nel_data,open('nelect_data.pkl','wb'))

    calc.bool_params['lwave']=True

def get_closest(ref,atoms,ind,mic=True):
    from ase.geometry import get_distances
    # find the index of the closest atom between two states
    # making sure that the symbol is the same
    if not mic:
        pbc=(False,False,False)
    else:
        pbc=(True,True,True)
    dists = []
    for atom in atoms:
        if atom.symbol != ref[ind].symbol:
            continue
        dist = get_distances(p1=(atom.x,atom.y,atom.z),
                            p2=(ref[ind].x,ref[ind].y,ref[ind].z),
                            cell=atoms.cell,
                            pbc=pbc)[1][0][0]
        dists.append((atom.index,dist))
    dists.sort(key=lambda x:x[-1])
    return dists[0][0]

def reindex_atoms(ref_atoms,reindex_atoms,manual_skip_atoms=[]):
    from ase.io import read
    # used to reindex atoms in reindex_atoms to match those in 
    # ref_atoms. Necessary for e.g. NEB interpolation.
    for atom in reindex_atoms:
        if atom.index in manual_skip_atoms:
            continue
        closest_ind = get_closest(ref_atoms,reindex_atoms,atom.index)
        if atom.index == closest_ind:
            continue
        else:
            pos_swap(reindex_atoms,closest_ind,atom.index)
    return reindex_atoms

def const_U_relax(atoms,calc,desired_U,ediffg=0.05):
    ############################################################################
    # Script to perform a geometry optimization at constant potential. This
    # routine, along with the other const_U routines in this package, are 
    # designed to handle checkpointing smoothly.
    ############################################################################
    #
    ############################################################################
    # Expects an atoms object, a calculator object, and a desired potential.
    # All const_U routines in this package use ASE -- see the ASE website for
    # more details on how to set up atoms and calculator objects.
    #
    # You should set the specific calculator parameters that you want in the 
    # calculator object that is passed to this routine (e.g. ENCUT, kpts, ...).
    ############################################################################

    import os,sys,pickle,math
    
    # ensure the force cutoff is set properly
    calc.float_params['ediffg'] = -1*ediffg
    atoms.set_calculator(calc)

    converged = 0
    i = 0
    while converged == 0:
        i += 1
        if i > 10:
            print('Stuck in a loop -- bug report?')
            exit()
        # first optimize NELECT
        set_pot(atoms,calc,desired_U)
        calc.int_params['nsw'] = 300
        calc.bool_params['lwave'] = True
        nel_data = pickle.load(open('./nelect_data.pkl','rb'))

        # geometry optimize using VASP optimizer
        print('Starting geometry optimization, iteration %i'%i)
        sys.stdout.flush()
        atoms.get_potential_energy() # calls VASP

        # update NELECT history
        nel_data['potential'].append(get_wf_implicit('./')-_she_U)
        nel_data['energy'].append(atoms.get_potential_energy())
        nel_out = float(greplines('grep NELECT OUTCAR')[0].split()[2])
        nel_data['nelect'].append(nel_out)
        pickle.dump(nel_data,open('nelect_data.pkl','wb'))

        # restart from CONTCAR
        os.system('cp CONTCAR POSCAR')
        atoms.write('iter%02d.traj'%i)

        # check convergence criteria: max forces, and current potential
        if fmax(atoms) < ediffg and abs(float(get_wf_implicit('./'))-_she_U - desired_U) < _tolerance_U:
            converged = 1
        else:
            print('Not yet converged')
            print('U = %.2f V vs SHE'%(float(get_wf_implicit('./'))-_she_U))
            print('max force = %.2f eV/A'%fmax(atoms))
        sys.stdout.flush()

    print('\nFinished!\n')


def const_U_dimer(atoms,calc,desired_U,ediffg=0.05,iopt=2):
    ############################################################################
    # Script to locate transition state at constant potential using
    # the Dimer method. See https://theory.cm.utexas.edu/vtsttools/dimer.html
    # for more details on the Dimer method.
    ############################################################################
    #
    ############################################################################
    # Expects an atoms object, a calculator object, and a desired potential.
    #
    # All const_U routines in this package use ASE -- see the ASE website for
    # more details on how to set up atoms and calculator objects.
    #
    # You should set the specific values of IOPT etc that you want in the 
    # calculator object that is passed to this routine.
    #
    # Ideally, you should also already have a MODECAR file created, as 
    # convergence is bad without a good initial guess in my experience.
    ############################################################################

    import os,sys,pickle,math

    # set required flags for Dimer method, if they're not already set
    calc.float_params['ediffg'] = -1*ediffg
    atoms.set_calculator(calc)

    converged = 0
    i = 0
    while converged == 0:
        i += 1
        if i > 10:
            print('Stuck in a loop -- bug report?')
            exit()

        # first optimize NELECT
        calc.int_params['ichain'] = 0

        if i > 1:
            calc.float_params['nelect'] += 1e-4 
            # avoid a weird bug going through the loops too fast
            # should probably figure out at some point ...
            # it's related to atoms.get_potential_energy() only calling
            # vasp to run if the atoms object has not changed at all
            # ... so if you change the calculator by some very small amount,
            # vasp will run instead of instantly returning the energy/forces
            # ¯\_(ツ)_/¯ 
            atoms.set_calculator(calc)
        set_pot(atoms,calc,desired_U)

        # ICHAIN 2 == Dimer method
        calc.int_params['ichain'] = 2

        # IBRION 3 == molecular dynamics
        # POTIM 0 == zero time step
        # this will ensure VASP uses VTST optimizers
        calc.int_params['ibrion'] = 3
        calc.float_params['potim'] = 0

        # IOPT is set to two by default in this function, which 
        # corresponds to a CG method. VTST sets it to 1 by default
        # which is L-BFGS -- in my experience this can hang forever
        # as it searches for a lower energy step.
        # Another recommended setting would be 7, which is FIRE. 
        # details: https://theory.cm.utexas.edu/vtsttools/optimizers.html
        calc.int_params['iopt'] = iopt

        # ediff = 1e-8 needed for accurate estimation of forces
        # this gets set to 1e-4 during the NELECT optimization routine
        calc.exp_params['ediff'] = 1.0e-8
        calc.int_params['nsw'] = 300
        calc.bool_params['lwave'] = True
        nel_data = pickle.load(open('./nelect_data.pkl','rb'))

        # Run Dimer calculation 
        print('Starting dimer optimization, iteration %i'%i)
        sys.stdout.flush()
        atoms.get_potential_energy()

        # update NELECT history
        nel_data['potential'].append(get_wf_implicit('./')-_she_U)
        nel_data['energy'].append(atoms.get_potential_energy())
        nel_out = float(greplines('grep NELECT OUTCAR')[0].split()[2])
        nel_data['nelect'].append(nel_out)
        pickle.dump(nel_data,open('nelect_data.pkl','wb'))

        # CENTCAR and CONTCAR should be similar...
        # but CENTCAR is technically the write one to restart from
        os.system('cp CENTCAR POSCAR')
        atoms.write('iter%02d.traj'%i)

        # check convergence criteria: max forces, and current potential
        if fmax(atoms) < ediffg and abs(float(get_wf_implicit('./'))-_she_U - desired_U) < _tolerance_U:
            converged = 1
        else:
            print('Not yet converged')
            print('U = %.2f V vs SHE'%(float(get_wf_implicit('./'))-_she_U))
            print('max force = %.2f eV/A'%fmax(atoms))
        sys.stdout.flush()

    print('\nFinished!\n')

def const_U_FBL(atoms,calc,desired_U,ind1,ind2,z_cutoff=None,ediffg=0.05):
    ############################################################################
    # Script to perform a geometry optimization with a fixed bond length 
    # constraint. This can be seen as an alternative to the dimer method, 
    # provided your reaction pathway is roughly one dimensional.
    ############################################################################
    #
    ############################################################################
    # Expects an atoms object, a calculator object, a desired potential, and
    # the indices of the two atoms to be fixed during geometry optimization.
    #
    # All const_U routines in this package use ASE -- see the ASE website for
    # more details on how to set up atoms and calculator objects.
    #
    # You should set the specific calculator parameters that you want in the 
    # calculator object that is passed to this routine (e.g. ENCUT, kpts, ...).
    #
    # Optional argument: z_cutoff specifices a z coordinate below which all
    # atoms will be fixed during geometry optimization. Alternately, you can
    # just fix the atoms you want and pass that atoms object into this routine.
    ############################################################################
    import os,sys,pickle,math,time
    from ase.constraints import FixBondLength,FixAtoms
    from ase.optimize import QuasiNewton,BFGS

    calc.float_params['ediffg'] = -1*ediffg
    atoms.set_calculator(calc)

    converged = 0
    i = 0
    while converged == 0:
        i += 1
        if i > 10:
            print('Stuck in a loop -- bug report?')
            exit()
        # first optimize NELECT
        set_pot(atoms,calc,desired_U)
        # calc.int_params['nsw'] = 300
        calc.bool_params['lwave'] = True
        nel_data = pickle.load(open('./nelect_data.pkl','rb'))

        # geometry optimize 
        print('Starting geometry optimization, iteration %i'%i)
        sys.stdout.flush()
        # need to set up constraints
        c = atoms.constraints
        fbl = FixBondLength(ind1,ind2)
        # only add FBL constraint if it's not already in the constraints list:
        if fbl not in c:
            c.append(fbl)

        # optionally, you can specify a z to fix atoms
        if z_cutoff is not None:
            fix_inds = [atom.index for atom in atoms if atom.z < z_cutoff]
            c.append(FixAtoms(indices=fix_inds))
        atoms.set_constraint(c)
        atoms.set_calculator(calc)

        # Now run optimization useing ASE optimizer.
        # In my experience, linesearch methods fail to converge,
        # which is why this uses regular old BFGS.
        dyn = BFGS(atoms,trajectory='./qn.traj',logfile='./qn.log')
        dyn.run(fmax=ediffg)

        # update NELECT history
        nel_data['potential'].append(get_wf_implicit('./')-_she_U)
        nel_data['energy'].append(atoms.get_potential_energy())
        nel_out = float(greplines('grep NELECT OUTCAR')[0].split()[2])
        nel_data['nelect'].append(nel_out)
        pickle.dump(nel_data,open('nelect_data.pkl','wb'))

        # restart from CONTCAR
        os.system('cp CONTCAR POSCAR')
        atoms.write('iter%02d.traj'%i)

        # check convergence criteria: max forces, and current potential
        if fmax(atoms) < ediffg and abs(float(get_wf_implicit('./'))-_she_U - desired_U) < _tolerance_U:
            converged = 1
        else:
            print('Not yet converged')
            print('U = %.2f V vs SHE'%(float(get_wf_implicit('./'))-_she_U))
            print('max force = %.2f eV/A'%fmax(atoms))
        sys.stdout.flush()

    print('\nFinished!\n')


def match_pbcs(fs_atoms,is_atoms,moving_atoms=[],tolerance=1.0):
    # This function tries to match atoms across a reaction
    # coordinate so that interpolation can be used without
    # minimum image convention (which sometimes breaks)
    #
    # takes as input the FS and IS atoms, and a list of 
    # atoms which you expect to move across the reaction 
    # coordinate
    #
    # returns an updated FS atoms

    assert len(is_atoms)==len(fs_atoms)

    for i in range(len(is_atoms)):
        if i in moving_atoms:
            continue
        
        is_x = is_atoms[i].x
        fs_x = fs_atoms[i].x

        is_y = is_atoms[i].y
        fs_y = fs_atoms[i].y

        while abs(fs_x-is_x) > tolerance:
            print('X difference: %.2f'%(abs(fs_x-is_x)))
            # need to move fs atoms in x
            if fs_x > is_x:
                fs_atoms[i].x -= (fs_atoms.cell[0][0]+fs_atoms.cell[1][0])
            elif fs_x < is_x: 
                fs_atoms[i].x += (fs_atoms.cell[0][0]+fs_atoms.cell[1][0])
            is_x = is_atoms[i].x
            fs_x = fs_atoms[i].x


        while abs(fs_y-is_y) > tolerance:
            print('Y difference: %.2f'%(abs(fs_y-is_y)))
            # need to move fs atoms in x
            if fs_y > is_y:
                fs_atoms[i].y -= (fs_atoms.cell[1][1])
            elif fs_y < is_y: 
                fs_atoms[i].y += (fs_atoms.cell[1][1])
            is_y = is_atoms[i].y
            fs_y = fs_atoms[i].y
    return fs_atoms
