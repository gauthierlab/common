
# package containing many of the scripts i've written - for ease of use

# list of dependencies:
# ase; numpy; scipy; matplotlib; spglib; prettytable; pymatgen; 

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
# get_n0(atoms)
# get_omega(path)
# match_cell(ref_atoms,change_atoms,lower_vac,anchor_atom=None)


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


def get_line(x,y,extra=0.1,extramin=0.0,extraplus=0.0):
    import sys,subprocess
    import numpy as np
    # returns a np array for the x and y axis of a line
    # given some data
    xax = np.linspace(min(x)-extra-extramin,max(x)+extra+extraplus,10)
    a,b = np.polyfit(x,y,1)
    yax = a*xax+b
    return xax,yax,a,b


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
        if atom.symbol == old_sym or atom.symbol == new_sym:
            atom.x *= ratio
            atom.y *= ratio
            atom.z *= ratio
        else:
            atom.x -= translate[0]
            atom.y -= translate[1]
            atom.z -= translate[2]
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

def get_n0(atoms):
    # assumes standard VASP PBE pseudopotentials
    n0 = 0
    for atom in atoms:
        if atom.symbol == 'Cu':
            n0 += 11
        elif atom.symbol == 'O':
            n0 += 6
        elif atom.symbol == 'C':
            n0 += 4
        elif atom.symbol == 'H':
            n0 += 1
    return n0

def get_omega(path):
    from ase.io import read
    # calculate the grand canonical energy given a path
    # run save_space.py first
    e = read('%s/lastimage.traj'%path).get_potential_energy()

    f = open('%s/nel.txt'%path,'r')
    nel = float(f.readlines()[0].rstrip())
    f.close()

    q = nel-get_n0(read('%s/lastimage.traj'%path))

    f = open('%s/fermi.txt'%path,'r')
    fermi = float(f.readlines()[0].rstrip())
    f.close()

    print(path,q,fermi)
    return e-q*fermi

def match_cell(ref_atoms,change_atoms,lower_vac,anchor_atom=None):
    from ase.io import read
    new_cell = ref_atoms.cell
    change_atoms.cell = new_cell
    if anchor_atom is not None:
        # rather than change positions according to a vacuum,
        # reference everything to the z-position of an atom
        # in the reference atoms object
        zdiff = change_atoms[anchor_atom] - ref_atoms[anchor_atom]
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
