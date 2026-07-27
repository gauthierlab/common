"""Common utilities for constant-potential DFT with VASP/VASPsol++ and ASE.

Shared across research groups. Import as `import common` (with
/Users/jgauth32/PythonModules or equivalent on PYTHONPATH), or install with
`pip install -e .`.

Use `help(common)` or `help(common.some_func)` for the function list --
there is deliberately no hand-maintained index here, since it drifts.

Module-level defaults (override globally with e.g. `common._she_U = 4.6`,
or per-call via the `she_U=`/`tolerance_U=` keyword arguments on set_pot
and the const_U_* drivers):
  _she_U       -- SHE reference potential (default 4.43 V)
  _tolerance_U -- potential convergence tolerance (default 0.02 V)

Dependencies: ase, numpy, scipy, matplotlib, pymatgen.
"""

# SHE reference potential
_she_U=4.43

# tolerance criteria used when optimizing NELECT for a given potential, in Volts
_tolerance_U=0.02

# Planck's constant in units of eV s
_h = 4.135667696e-15

#Boltzmann's constant in units of eV
_kb = 8.617333262e-5


def get_irr_kpts(atoms,kpts,is_shift=[0,0,0]):
    """Returns the number of irreducible kpoints given an atoms object and a desired k-pt grid."""
    import numpy as np
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga
    from pymatgen.io.ase import AseAtomsAdaptor
    import os
    aaa = AseAtomsAdaptor
    struc = aaa.get_structure(atoms)
    mesh = sga(struc).get_ir_reciprocal_mesh(mesh=kpts)
    # atoms.write('./tmp.cif')
    # mesh = sga(structure=Structure.from_file('./tmp.cif')).get_ir_reciprocal_mesh(mesh=kpts)
    # os.system('rm ./tmp.cif')
    return len(mesh)

    # kpts = np.array(kpts)
    # lattice = np.array([atoms.cell[0],atoms.cell[1],atoms.cell[2]])
    # positions = atoms.get_positions()
    # numbers = [1,]*len(positions)
    # cell = (lattice,positions,numbers)

    # mapping,grid = spglib.get_ir_reciprocal_mesh(kpts, cell, is_shift=is_shift)
    # return len(np.unique(mapping))


def get_line(x,y,extra=0.1,extramin=0.0,extraplus=0.0,return_mae=False):
    """Returns a np array for the x and y axis of a line given some data."""
    import sys,subprocess
    import numpy as np
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


def _fgrep(filepath, keyword):
    """Return list of lines containing keyword; safe for paths with spaces."""
    matches = []
    try:
        with open(filepath, errors='replace') as f:
            for line in f:
                if keyword in line:
                    matches.append(line.rstrip('\n'))
    except (OSError, IOError):
        pass
    return matches

def greplines(cmd):
    """Easier subprocess use - auto split by newline character."""
    import sys,subprocess
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
    """Easier subprocess use."""
    import sys,subprocess
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
    """Read in INCAR file, set parameter to desired level."""
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
    """Compute the work function for an Environ-based calculation, if used."""
    import os,sys,subprocess,pickle
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
        if not os.path.exists('%s/vasprun.xml'%path):
            # no vasprun either
            if os.path.exists('%s/fermi.txt'%path):
                fermi = float(open('%s/fermi.txt'%path).read().strip())
            else:
                fermi = float(greplines('cat %s/fermi.txt'%path)[0])
            if not os.path.exists('%s/INCAR'%path):
                return -1*fermi  # INCAR stripped; all constant-U jobs here are VASPsol++
            solcheck = _fgrep('%s/INCAR'%path, 'ISOL')
        else:
            print('No OUTCAR found -- use vasprun.xml instead')
            fermi = float(_fgrep('%s/vasprun.xml'%path, 'fermi')[0].split()[-2])
            solcheck = _fgrep('%s/vasprun.xml'%path, 'ISOL')
    else:
        fermi = float(_fgrep('%s/OUTCAR'%path, 'fermi')[-1].split()[2])
        solcheck = _fgrep('%s/OUTCAR'%path, 'ISOL')

    if solcheck != []:
        # using VASPsol++, no need to check for FERMI_SHIFT
        shift = 0.0
        return -1*fermi

    try:
        out2 = _fgrep('%s/opt.log'%path, 'FERMI_SHIFT')
        shift = float(out2[-1].split(' = ')[-1])
    except (IndexError, ValueError):
        try:
            out2 = _fgrep('%s/vasp.out'%path, 'FERMI_SHIFT')
            shift = float(out2[-1].split(' = ')[-1])
        except (IndexError, ValueError):
            print('Error: could not find FERMI_SHIFT in %s/opt.log or %s/vasp.out' % (path, path))
            import sys
            sys.exit()

    return -1*(fermi+shift)

def get_wf_explicit(path):
    import os
    if os.path.isfile(path+'/wf.out'):
        f = open(path+'/wf.out')
        lines = f.readlines()
        return float(lines[0].rstrip())

def get_chgcar(path,locpot=False,pavg=False):
    """Get CHGCAR data as a numpy array.

    Optionally, return the planar average of this density. If the file is
    named "LOCPOT", the flag 'locpot' will be set to True.
    """
    import os
    import numpy as np
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
    out = greplines('grep -n augmentation ' + path)[0].split('\n')[0]
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
    """Swap position of two atoms."""
    x1,y1,z1 = [atoms[ind1].x,atoms[ind1].y,atoms[ind1].z]
    x2,y2,z2 = [atoms[ind2].x,atoms[ind2].y,atoms[ind2].z]

    atoms[ind1].x = x2
    atoms[ind1].y = y2
    atoms[ind1].z = z2

    atoms[ind2].x = x1
    atoms[ind2].y = y1
    atoms[ind2].z = z1

    return atoms

def scale_metal(atoms,a_old,a_new,old_sym=None,new_sym=None):
    """Useful for changing lattice constant of metal, or changing the metal
    being used in a trajectory.

    Generally will only work for orthorhombic cells.
    """
    ratio = a_new/a_old

    if new_sym is not None:
        assert old_sym is not None, "Need to specify which symbol to replace with new_sym using keyword 'old_sym'"
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


def fed(path,label=None,fig=None):
    """Given a figure object, a reaction path, and a label, returns a free
    energy diagram (FED) of the reaction path.

    The reaction path should be a list of tuples, with the first element of
    each tuple being the energy relative to the first entry, and the second
    element being either 'min' or 'ts' for the respective types (energy
    minimum or transition state).

    e.g.,
        fig = plt.figure(figsize=(6*1.618,6))
        label = 'CO + COH -> OCCOH'
        path = [(0,'min'),(1.5,'ts'),(1.0,'min')]
    """
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
    pbe_dict = {'Pd': 10.0, 'Sb': 5.0, 'Cr': 12.0, 'Se': 6.0, 'Sn': 14.0, 'Li': 3.0, 'He': 2.0, 'Fr': 9.0, 'Cs': 9.0, 'Nd': 11.0, 'Ac': 11.0, 'Ho': 9.0, 'Eu': 8.0, 'Ni': 10.0, 'Po': 16.0, 'Am': 17.0, 'Be': 2.0, 'Sr': 10.0, 'Al': 3.0, 'Mg': 2.0, 'Ir': 9.0, 'Ge': 14.0, 'Sm': 11.0, 'Cu': 11.0, 'Ra': 10.0, 'Hf': 10.0, 'Co': 9.0, 'Fe': 8.0, 'Ga': 13.0, 'Ba': 10.0, 'Te': 6.0, 'U': 14.0, 'Tb': 9.0, 'I': 7.0, 'Er': 9.0, 'N': 5.0, 'Rn': 8.0, 'Ca': 10.0, 'Nb': 13.0, 'S': 6.0, 'Tl': 13.0, 'F': 7.0, 'O': 6.0, 'Ta': 11.0, 'Pb': 14.0, 'H': 1.0, 'Zn': 12.0, 'Na': 7.0, 'Pu': 16.0, 'Gd': 9.0, 'Pt': 10.0, 'Sc': 11.0, 'V': 13.0, 'Lu': 9.0, 'Dy': 9.0, 'Pa': 13.0, 'Si': 4.0, 'Ag': 11.0, 'Kr': 8.0, 'Pm': 11.0, 'Tc': 13.0, 'Ar': 8.0, 'Rb': 9.0, 'Au': 11.0, 'W': 12.0, 'Ne': 8.0, 'At': 7.0, 'Np': 15.0, 'Tm': 9.0, 'As': 5.0, 'Hg': 12.0, 'K': 9.0, 'Br': 7.0, 'Os': 8.0, 'Yb': 8.0, 'Cd': 12.0, 'Cm': 18.0, 'Pr': 11.0, 'Ru': 8.0, 'Mo': 14.0, 'In': 13.0, 'Cl': 7.0, 'La': 11.0, 'Ce': 12.0, 'C': 4.0, 'Th': 12.0, 'B': 3.0, 'Y': 11.0, 'Mn': 13.0, 'Bi': 15.0, 'Re': 7.0, 'Xe': 8.0, 'Rh': 9.0, 'Zr': 12.0, 'Ti': 12.0, 'P': 5.0}
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
    """Calculate the grand canonical energy given a path.

    Run save_space.py first.
    """
    from ase.io import read
    import os
    n0 = get_n0(path)

    if not os.path.exists('%s/OUTCAR'%path):
        # print('no OUTCAR in directory, trying vasprun.xml instead')
        try:
            e = read('%s/vasprun.xml'%path).get_potential_energy()
        except:
            e = read('%s/lastimage.traj'%path).get_potential_energy()

        if os.path.exists('%s/vasprun.xml'%path):
            nel = float(_fgrep('%s/vasprun.xml'%path, 'NELECT')[0].split()[-1][:-4])
            fermi = float(_fgrep('%s/vasprun.xml'%path, 'fermi')[0].split()[-2])
        else:
            # newer ASE doesn't write necessary files to vasprun.xml
            # maybe you saved necessary info to a text file to save space?
            nel = float(open('%s/nel.txt'%path).read().strip())
            fermi = float(open('%s/fermi.txt'%path).read().strip())
    else:
        try:
            e = read('%s/OUTCAR'%path).get_potential_energy()
        except:
            e = read('%s/vasprun.xml'%path).get_potential_energy()
        nel = float(_fgrep('%s/OUTCAR'%path, 'NELECT')[0].split()[2])
        fermi = float(_fgrep('%s/OUTCAR'%path, 'fermi')[-1].split()[2])
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
    """Given an atoms object, return the maximum force on an unconstrained atom."""
    import numpy as np
    from ase.constraints import FixAtoms
    fixed = set()
    for c in atoms.constraints:
        if isinstance(c, FixAtoms):
            fixed.update(c.index)
    unconstrained = [atom.index for atom in atoms if atom.index not in fixed]
    forces = atoms.get_forces()
    return max(np.linalg.norm(forces[i]) for i in unconstrained)

def _next_nelect(nel_data, desired_U):
    """Pure helper: given the NELECT/potential history so far, return the
    proposed next NELECT value for set_pot's Newton-iteration loop.

    Mirrors the logic previously inlined in set_pot's while loop:
      - if the last two 'nelect' entries are equal, perturb by +0.1 to
        escape (the two-point gradient estimate would be undefined).
      - elif fewer than 3 potential entries: estimate the gradient from
        the last two points; if the NELECT denominator is too small to
        trust, perturb by +0.1 instead.
      - else: estimate the gradient via a 3-point linear fit
        (common.get_line) over the last three (nelect, potential) points.
      - apply a Newton step using the estimated gradient, clamped to
        +/-0.75 if the raw step would exceed +/-5.0.
    """
    last_nelect = nel_data['nelect'][-1]

    if nel_data['nelect'][-1] == nel_data['nelect'][-2]:
        # NELECT didn't change -- perturb by a meaningful amount
        # based on the last known gradient to escape
        return last_nelect + 0.1

    if len(nel_data['potential']) < 3:
        # only two points to estimate gradient
        grad_numer = nel_data['potential'][-2]-nel_data['potential'][-1]
        grad_denom = nel_data['nelect'][-2]-nel_data['nelect'][-1]
        if abs(grad_denom) < 0.0001:
            # gradient is unreliable, take a small fixed step instead
            return last_nelect + 0.1
        # if grad_denom is not almost zero, calculate a gradient normally
        grad = grad_numer/grad_denom
    else:
        # use the last 3 points rather than last two to estimate gradient
        x,y,grad,intercept = get_line(nel_data['nelect'][-3:],nel_data['potential'][-3:])

    y = nel_data['potential'][-1]-desired_U
    diff = abs(y)**2/(y*grad)  # how much to change NELECT by this step

    # don't take too big of a step ..
    # can happen if two subsequent steps are too close together
    # again, shouldn't be possible, but just make sure
    if diff > 5.0:
        diff = 0.75
    elif diff < -5.0:
        diff = -0.75

    return last_nelect - diff

def set_pot(atoms,calc,desired_U,she_U=None,tolerance_U=None):
    """Determine NELECT required to have potential=desired_U.

    Optional arguments:
    she_U       -- SHE reference potential, in V. Default (None) resolves
                   to the module-level common._she_U at call time.
    tolerance_U -- potential convergence tolerance, in V. Default (None)
                   resolves to the module-level common._tolerance_U at
                   call time.
    """
    import os,sys,pickle,math
    from ase.io import read
    import numpy as np
    if she_U is None:
        she_U = _she_U
    if tolerance_U is None:
        tolerance_U = _tolerance_U
    calc.bool_params['lcharg'] = False
    calc.int_params['ichain'] = 0
    calc.int_params['iopt'] = 0
    calc.bool_params['lwave'] = True
    calc.int_params['nsw'] = 0
    calc.exp_params['ediff'] = 1.0e-4
    calc.int_params['ibrion'] = 1
    if calc.bool_params['lsol'] == False:
        print('Continuum solvation is disabled -- are you sure this is correct?')

    # previous optimization was done, use that as starting point
    if os.path.isfile('nelect_data.pkl') and os.stat('nelect_data.pkl').st_size != 0:
        nel_data = pickle.load(open('./nelect_data.pkl','rb'))
        nel_data['she_U'] = she_U
        nel_data['tolerance_U'] = tolerance_U
        calc.float_params['nelect'] = nel_data['nelect'][-1]
        atoms.set_calculator(calc)
    else:
        nel_data = {}
        nel_data['nelect'] = []
        nel_data['potential'] = []
        nel_data['energy'] = []
        nel_data['she_U'] = she_U
        nel_data['tolerance_U'] = tolerance_U
        print('Running the first single point to get PZC')
        atoms.set_calculator(calc)
        atoms.get_potential_energy()

        # store info from the first single point
        nel_data['potential'].append(get_wf_implicit('./')-she_U)
        nel_data['energy'].append(atoms.get_potential_energy())
        nel_out = float(greplines('grep NELECT OUTCAR')[0].split()[2])
        nel_data['nelect'].append(nel_out)
        pickle.dump(nel_data,open('nelect_data.pkl','wb'))

    # no need to run further optimization if you're already at the desired potential
    if abs(nel_data['potential'][-1]-desired_U) < tolerance_U:
        return

    if len(nel_data['nelect']) < 2:
        # only one data point - do another single point with slightly more
        # electrons to get an initial gradient for newton's method
        # initial guess for C: 1 e/V in a ~60 square angstrom cell
        area = np.linalg.norm(np.cross(atoms.cell[0],atoms.cell[1]))
        calc.float_params['nelect'] = nel_data['nelect'][-1]+area/60.0*(nel_data['potential'][-1]-desired_U)
        atoms.set_calculator(calc)
        atoms.get_potential_energy()

        nel_data['potential'].append(get_wf_implicit('./')-she_U)
        nel_data['energy'].append(atoms.get_potential_energy())
        nel_out = float(greplines('grep NELECT OUTCAR')[0].split()[2])
        nel_data['nelect'].append(nel_out)
        pickle.dump(nel_data,open('nelect_data.pkl','wb'))

    #start the optimization, initialize vars
    max_iter = 20
    n_iter = 0
    while abs(nel_data['potential'][-1]-desired_U) > tolerance_U:
        n_iter += 1
        if n_iter > max_iter:
            print('Error: set_pot did not converge after %d iterations' % max_iter)
            print('Last potential: %.4f V, desired: %.4f V' % (nel_data['potential'][-1], desired_U))
            sys.exit()

        # Newton's method to optimize NELECT (see _next_nelect for the
        # branch logic: equal-nelect escape, two-point gradient with a
        # tiny-denominator escape, or 3-point linear-fit gradient)
        new_nel = _next_nelect(nel_data, desired_U)

        #check if nelect is nan
        if math.isnan(new_nel):
            print('Error: Check NELECT (nan)')
            sys.exit()

        # run VASP with the proposed NELECT and record the result
        calc.float_params['nelect'] = new_nel
        atoms.set_calculator(calc)
        atoms.get_potential_energy()

        nel_data['potential'].append(get_wf_implicit('./')-she_U)
        nel_data['energy'].append(atoms.get_potential_energy())
        nel_out = float(greplines('grep NELECT OUTCAR')[0].split()[2])
        nel_data['nelect'].append(nel_out)
        pickle.dump(nel_data,open('nelect_data.pkl','wb'))

    calc.bool_params['lwave']=True

def get_closest(ref,atoms,ind,mic=True):
    """Find the index of the closest atom between two states, making sure
    that the symbol is the same.
    """
    from ase.geometry import get_distances
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
    """Used to reindex atoms in reindex_atoms to match those in ref_atoms.
    Necessary for e.g. NEB interpolation.
    """
    from ase.io import read
    for atom in reindex_atoms:
        if atom.index in manual_skip_atoms:
            continue
        closest_ind = get_closest(ref_atoms,reindex_atoms,atom.index)
        if atom.index == closest_ind:
            continue
        else:
            pos_swap(reindex_atoms,closest_ind,atom.index)
    return reindex_atoms

# files worth snapshotting between runs. WAVECAR/CHGCAR are deliberately
# excluded -- they are large and are not useful after the geometry moves.
_BACKUP_FILES = ['INCAR','KPOINTS','POSCAR','CONTCAR','OUTCAR','OSZICAR',
                 'XDATCAR','vasprun.xml','opt.log','nelect_data.pkl']

# extra files written by the VTST dimer method (ICHAIN=2)
_DIMER_BACKUP_FILES = _BACKUP_FILES + ['CENTCAR','DIMCAR','MODECAR','NEWMODECAR']


def _next_backup_path(prefix='backup'):
    """Return the first unused '<prefix>_NN' directory name in the cwd.

    Numbering continues past the highest existing index rather than
    counting directories, so a gap (e.g. backup_01 deleted by hand) does
    not cause a collision with an existing backup.
    """
    import os,glob
    i = 0
    for d in glob.glob(prefix+'_[0-9][0-9]'):
        try:
            i = max(i,int(d.rsplit('_',1)[-1])+1)
        except ValueError:
            continue
    while os.path.exists('%s_%02d'%(prefix,i)):
        i += 1
    return '%s_%02d'%(prefix,i)


def make_backup(files=None,extra=(),prefix='backup',path=None,quiet=False):
    """Copy the given files into a fresh backup directory, and return its path.

    Files that do not exist are skipped silently -- the file lists are
    supersets covering both plain relaxations and dimer runs, so a missing
    DIMCAR (say) is normal rather than an error.

    Optional arguments:
    files  -- list of filenames to copy. Default (None) is _BACKUP_FILES.
    extra  -- additional filenames to copy, appended to `files`. Handy for
              per-iteration artifacts such as 'iter03.traj'.
    prefix -- directory prefix, so backups land in <prefix>_00, _01, ...
              Default 'backup'.
    path   -- explicit destination directory. Default (None) picks the next
              free '<prefix>_NN'. Must not already exist.
    quiet  -- suppress the one-line summary printed on success.
    """
    import os,shutil,sys
    if files is None:
        files = _BACKUP_FILES
    if path is None:
        path = _next_backup_path(prefix)
    os.makedirs(path)

    n = 0
    for f in list(files)+list(extra):
        if os.path.isfile(f):
            shutil.copy2(f,os.path.join(path,os.path.basename(f)))
            n += 1
    if not quiet:
        print('Backed up %i files to %s'%(n,path))
        sys.stdout.flush()
    return path


def const_U_relax(atoms,calc,desired_U,ediffg=0.05,optimizer='vasp',iopt=2,ediff=None,she_U=None,tolerance_U=None):
    """Script to perform a geometry optimization at constant potential. This
    routine, along with the other const_U routines in this package, are
    designed to handle checkpointing smoothly.

    Expects an atoms object, a calculator object, and a desired potential.
    All const_U routines in this package use ASE -- see the ASE website for
    more details on how to set up atoms and calculator objects.

    You should set the specific calculator parameters that you want in the
    calculator object that is passed to this routine (e.g. ENCUT, kpts, ...).

    Optional arguments:
    optimizer   -- 'vasp' (default) uses VASP's internal optimizer, exactly
                   as before. 'vtst' hands the geometry steps to the VTST
                   optimizer compiled into VASP via IBRION=3/POTIM=0/IOPT
                   (same pattern as const_U_dimer, but with ICHAIN=0).
    iopt        -- VTST optimizer selector, only used when optimizer='vtst'.
                   2 = CG (default), 7 = FIRE.
                   https://theory.cm.utexas.edu/vtsttools/optimizers.html
    ediff       -- electronic convergence for the geometry steps. Default
                   (None) keeps current behavior on the 'vasp' path and uses
                   1e-5 on the 'vtst' path.
    she_U       -- SHE reference potential, in V. Default (None) resolves
                   to the module-level common._she_U at call time. Passed
                   through to set_pot.
    tolerance_U -- potential convergence tolerance, in V. Default (None)
                   resolves to the module-level common._tolerance_U at
                   call time. Passed through to set_pot.
    """

    import os,sys,pickle,math

    if she_U is None:
        she_U = _she_U
    if tolerance_U is None:
        tolerance_U = _tolerance_U

    if optimizer not in ('vasp','vtst'):
        raise ValueError("optimizer must be 'vasp' or 'vtst', got %r" % (optimizer,))

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

        if optimizer == 'vtst' and i > 1:
            # force ASE to re-run VASP by clearing cached results
            # (same workaround as const_U_dimer)
            if hasattr(atoms, 'calc') and hasattr(atoms.calc, 'results'):
                atoms.calc.results = {}

        # first optimize NELECT
        # (set_pot resets IBRION=1, ICHAIN=0, IOPT=0, EDIFF=1e-4, NSW=0)
        set_pot(atoms,calc,desired_U,she_U=she_U,tolerance_U=tolerance_U)

        if optimizer == 'vtst':
            # hand geometry steps to the VTST optimizer compiled into VASP:
            # IBRION 3 == molecular dynamics, POTIM 0 == zero time step so
            # VASP itself never moves the ions; ICHAIN 0 == no chain method
            # (plain relaxation, unlike the dimer); IOPT selects the VTST
            # optimizer (2 = CG, 7 = FIRE). VTST respects the POSCAR
            # selective dynamics flags, so FixAtoms constraints still apply.
            # details: https://theory.cm.utexas.edu/vtsttools/optimizers.html
            calc.int_params['ichain'] = 0
            calc.int_params['ibrion'] = 3
            calc.float_params['potim'] = 0
            calc.int_params['iopt'] = iopt
            # force-based optimizers want forces converged tighter than the
            # 1e-4 that set_pot leaves behind, but the VASPsol++ LPB solver
            # bases its convergence cutoff on EDIFF and struggles much below
            # 1e-5, so default to 1e-5 rather than the dimer's 1e-8
            calc.exp_params['ediff'] = 1e-5 if ediff is None else ediff
        elif ediff is not None:
            calc.exp_params['ediff'] = ediff

        calc.int_params['nsw'] = 300
        calc.bool_params['lwave'] = True
        nel_data = pickle.load(open('./nelect_data.pkl','rb'))

        # geometry optimize using VASP optimizer
        print('Starting geometry optimization, iteration %i'%i)
        sys.stdout.flush()
        atoms.get_potential_energy() # calls VASP

        # update NELECT history
        nel_data['potential'].append(get_wf_implicit('./')-she_U)
        nel_data['energy'].append(atoms.get_potential_energy())
        nel_out = float(greplines('grep NELECT OUTCAR')[0].split()[2])
        nel_data['nelect'].append(nel_out)
        pickle.dump(nel_data,open('nelect_data.pkl','wb'))

        # restart from CONTCAR
        os.system('cp CONTCAR POSCAR')
        atoms.write('iter%02d.traj'%i)

        # check convergence criteria: max forces, and current potential
        if fmax(atoms) < ediffg and abs(float(get_wf_implicit('./'))-she_U - desired_U) < tolerance_U:
            converged = 1
        else:
            print('Not yet converged')
            print('U = %.2f V vs SHE'%(float(get_wf_implicit('./'))-she_U))
            print('max force = %.2f eV/A'%fmax(atoms))
        sys.stdout.flush()

    print('\nFinished!\n')


def const_U_dimer(atoms,calc,desired_U,ediff=1e-5,ediffg=0.05,iopt=2,she_U=None,tolerance_U=None,backup=True,backup_files=None):
    """Script to locate transition state at constant potential using
    the Dimer method. See https://theory.cm.utexas.edu/vtsttools/dimer.html
    for more details on the Dimer method.

    Expects an atoms object, a calculator object, and a desired potential.

    All const_U routines in this package use ASE -- see the ASE website for
    more details on how to set up atoms and calculator objects.

    You should set the specific values of IOPT etc that you want in the
    calculator object that is passed to this routine.

    Ideally, you should also already have a MODECAR file created, as
    convergence is bad without a good initial guess in my experience.

    Optional arguments:
    she_U       -- SHE reference potential, in V. Default (None) resolves
                   to the module-level common._she_U at call time. Passed
                   through to set_pot.
    tolerance_U -- potential convergence tolerance, in V. Default (None)
                   resolves to the module-level common._tolerance_U at
                   call time. Passed through to set_pot.
    backup      -- if True (default), snapshot each NSW=300 dimer run into
                   a fresh backup_NN directory before the restart files are
                   overwritten. Set False to disable.
    backup_files -- list of filenames to snapshot. Default (None) uses
                   common._DIMER_BACKUP_FILES. Pass a shorter list to save
                   disk if vasprun.xml/OUTCAR get unwieldy.
    """

    import os,sys,pickle,math

    if she_U is None:
        she_U = _she_U
    if tolerance_U is None:
        tolerance_U = _tolerance_U
    if backup_files is None:
        backup_files = _DIMER_BACKUP_FILES

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
        calc.int_params['iopt'] = 0

        if i > 1:
            # force ASE to re-run VASP by clearing cached results
            if hasattr(atoms, 'calc') and hasattr(atoms.calc, 'results'):
                atoms.calc.results = {}
        set_pot(atoms,calc,desired_U,she_U=she_U,tolerance_U=tolerance_U)

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

        # a tight ediff is needed for accurate estimation of forces
        # (it gets set to 1e-4 during the NELECT optimization routine),
        # but the VASPsol++ LPB solver bases its convergence cutoff on
        # EDIFF and struggles much below 1e-5, hence the 1e-5 default
        calc.exp_params['ediff'] = ediff
        calc.int_params['nsw'] = 300
        calc.bool_params['lwave'] = True
        nel_data = pickle.load(open('./nelect_data.pkl','rb'))

        # Run Dimer calculation 
        print('Starting dimer optimization, iteration %i'%i)
        sys.stdout.flush()
        atoms.get_potential_energy()

        # update NELECT history
        nel_data['potential'].append(get_wf_implicit('./')-she_U)
        nel_data['energy'].append(atoms.get_potential_energy())
        nel_out = float(greplines('grep NELECT OUTCAR')[0].split()[2])
        nel_data['nelect'].append(nel_out)
        pickle.dump(nel_data,open('nelect_data.pkl','wb'))

        atoms.write('iter%02d.traj'%i)

        # Snapshot the run *before* the restart shuffle below overwrites
        # POSCAR/MODECAR, so backup_NN holds a self-consistent record of
        # this pass: the geometry and mode it started from, plus the
        # OUTCAR/DIMCAR/CENTCAR/NEWMODECAR it produced. One backup per
        # iteration is exactly one NSW=300 dimer run at fixed NELECT --
        # i.e. the point where the potential gets reset and the search
        # resumes. The next call to set_pot clobbers OUTCAR/DIMCAR with
        # NSW=0 single points, which is why this cannot wait.
        if backup:
            make_backup(backup_files,extra=['iter%02d.traj'%i])

        # CENTCAR and CONTCAR should be similar...
        # but CENTCAR is technically the write one to restart from
        os.system('cp CENTCAR POSCAR')
        os.system('cp MODECAR oldMODECAR')
        os.system('cp NEWMODECAR MODECAR')

        # check convergence criteria: max forces, and current potential
        if fmax(atoms) < ediffg and abs(float(get_wf_implicit('./'))-she_U - desired_U) < tolerance_U:
            converged = 1
        else:
            print('Not yet converged')
            print('U = %.2f V vs SHE'%(float(get_wf_implicit('./'))-she_U))
            print('max force = %.2f eV/A'%fmax(atoms))
        sys.stdout.flush()

    print('\nFinished!\n')

def const_U_FBL(atoms,calc,desired_U,ind1,ind2,z_cutoff=None,ediffg=0.05,optimizer='ase',ase_optimizer='bfgs',she_U=None,tolerance_U=None):
    """Script to perform a geometry optimization with a fixed bond length
    constraint. This can be seen as an alternative to the dimer method,
    provided your reaction pathway is roughly one dimensional.

    Expects an atoms object, a calculator object, a desired potential, and
    the indices of the two atoms to be fixed during geometry optimization.

    All const_U routines in this package use ASE -- see the ASE website for
    more details on how to set up atoms and calculator objects.

    You should set the specific calculator parameters that you want in the
    calculator object that is passed to this routine (e.g. ENCUT, kpts, ...).

    Optional argument: z_cutoff specifices a z coordinate below which all
    atoms will be fixed during geometry optimization. Alternately, you can
    just fix the atoms you want and pass that atoms object into this routine.

    Optional arguments:
    optimizer     -- must be 'ase'. 'vtst' is rejected with an explanation:
                      the bond-length constraint only exists on the ASE side,
                      so an in-VASP optimizer would silently drop it.
    ase_optimizer -- 'bfgs' (default, current behavior) or 'fire'. FIRE is
                      the same algorithm as VTST IOPT=7 but runs on the ASE
                      side, so the bond constraint stays enforced.
    she_U         -- SHE reference potential, in V. Default (None) resolves
                      to the module-level common._she_U at call time.
                      Passed through to set_pot.
    tolerance_U   -- potential convergence tolerance, in V. Default (None)
                      resolves to the module-level common._tolerance_U at
                      call time. Passed through to set_pot.
    """
    import os,sys,pickle,math,time
    from ase.constraints import FixBondLength,FixAtoms
    from ase.optimize import QuasiNewton,BFGS,FIRE

    if she_U is None:
        she_U = _she_U
    if tolerance_U is None:
        tolerance_U = _tolerance_U

    if optimizer == 'vtst':
        raise ValueError(
            'const_U_FBL cannot use a VTST (in-VASP) optimizer: the fixed bond '
            'length is enforced by ase.constraints.FixBondLength on the Python '
            'side only. VTST optimizers (IBRION=3/POTIM=0/IOPT) move the ions '
            "inside VASP, which never sees this constraint, and VASP's ICONST "
            'constraint file is only active in molecular dynamics (IBRION=0), '
            'not in VTST optimizations, so the constraint would be silently '
            'dropped. The cost is real -- ASE-side stepping cannot use '
            "VASP's density/wavefunction extrapolation between steps -- but "
            'a native FBL algorithm would have to be added to VTST itself '
            "(Fortran). If BFGS is failing to converge, use optimizer='ase' "
            "with ase_optimizer='fire'.")
    if optimizer != 'ase':
        raise ValueError("optimizer must be 'ase' for const_U_FBL, got %r" % (optimizer,))
    if ase_optimizer not in ('bfgs','fire'):
        raise ValueError("ase_optimizer must be 'bfgs' or 'fire', got %r" % (ase_optimizer,))

    calc.float_params['ediffg'] = -1*ediffg
    atoms.set_calculator(calc)

    # set up constraints once, before the loop
    c = atoms.constraints
    fbl = FixBondLength(ind1,ind2)
    if fbl not in c:
        c.append(fbl)
    if z_cutoff is not None:
        fix_inds = [atom.index for atom in atoms if atom.z < z_cutoff]
        c.append(FixAtoms(indices=fix_inds))
    atoms.set_constraint(c)

    converged = 0
    i = 0
    while converged == 0:
        i += 1
        if i > 10:
            print('Stuck in a loop -- bug report?')
            exit()
        # first optimize NELECT
        set_pot(atoms,calc,desired_U,she_U=she_U,tolerance_U=tolerance_U)
        # calc.int_params['nsw'] = 300
        calc.bool_params['lwave'] = True
        nel_data = pickle.load(open('./nelect_data.pkl','rb'))

        # geometry optimize
        print('Starting geometry optimization, iteration %i'%i)
        sys.stdout.flush()
        atoms.set_calculator(calc)

        # Run optimization using an ASE optimizer, which enforces the
        # FixBondLength constraint each step. In my experience, linesearch
        # methods fail to converge, so the choices are regular old BFGS
        # (default) or FIRE (same algorithm as VTST IOPT=7, but the bond
        # constraint stays enforced on the ASE side).
        if ase_optimizer == 'fire':
            dyn = FIRE(atoms,trajectory='./qn.traj',logfile='./qn.log')
        else:
            dyn = BFGS(atoms,trajectory='./qn.traj',logfile='./qn.log')
        dyn.run(fmax=ediffg)

        # update NELECT history
        nel_data['potential'].append(get_wf_implicit('./')-she_U)
        nel_data['energy'].append(atoms.get_potential_energy())
        nel_out = float(greplines('grep NELECT OUTCAR')[0].split()[2])
        nel_data['nelect'].append(nel_out)
        pickle.dump(nel_data,open('nelect_data.pkl','wb'))

        # restart from CONTCAR
        os.system('cp CONTCAR POSCAR')
        atoms.write('iter%02d.traj'%i)

        # check convergence criteria: max forces, and current potential
        if fmax(atoms) < ediffg and abs(float(get_wf_implicit('./'))-she_U - desired_U) < tolerance_U:
            converged = 1
        else:
            print('Not yet converged')
            print('U = %.2f V vs SHE'%(float(get_wf_implicit('./'))-she_U))
            print('max force = %.2f eV/A'%fmax(atoms))
        sys.stdout.flush()

    print('\nFinished!\n')


def match_pbcs(fs_atoms,is_atoms,moving_atoms=[],tolerance=1.0):
    """Try to match atoms across a reaction coordinate so that interpolation
    can be used without minimum image convention (which sometimes breaks).

    Takes as input the FS and IS atoms, and a list of atoms which you expect
    to move across the reaction coordinate.

    Returns an updated FS atoms.
    """
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

def get_interp(atoms,end_inds,interp_ind,bl_1,bl_2,n_images=15):
    # linearly interpolates interp_ind between the two end_inds
    # such that desired bond lengths bl_1 and bl_2 at the end
    # points are satisfied. Useful for FBL calculations.
    # expected inputs:
    # atoms         -- an ASE atoms object containing the geometry to be used
    # end_inds      -- a list of the two end indices that interp_ind should 
    #                   be moved between
    # interp_ind    -- the atom to be moved between end_inds (e.g., H)
    # bl_1          -- desired bond length between end_inds[0] and interp_ind
    # bl_2          -- desired bond length between end_inds[1] and interp_ind
    # n_images      -- number of images to be created

	import os
	from scipy.optimize import fmin
	import numpy as np

	def f(t,a1,b1,a2,b2,a3,b3,ind,interp_ind,atoms,bl_des):
		# inner function used to find t parameter
		# corresponding to the end two bond lengths
		atoms[interp_ind].x = a1*t+b1
		atoms[interp_ind].y = a2*t+b2
		atoms[interp_ind].z = a3*t+b3
		dist = (bl_des-atoms.get_distance(ind,interp_ind))**2
		return dist

	# define three dimensional vector connecting
	# the two end points
	a1 = atoms[end_inds[1]].x - atoms[end_inds[0]].x
	a2 = atoms[end_inds[1]].y - atoms[end_inds[0]].y
	a3 = atoms[end_inds[1]].z - atoms[end_inds[0]].z
	b1 = atoms[end_inds[0]].x
	b2 = atoms[end_inds[0]].y
	b3 = atoms[end_inds[0]].z

	# get vector between the two end points
	t_end1 = fmin(f,0.5,args=(a1,b1,a2,b2,a3,b3,end_inds[0],interp_ind,atoms,bl_1),disp=False)[0]
	t_end2 = fmin(f,0.5,args=(a1,b1,a2,b2,a3,b3,end_inds[1],interp_ind,atoms,bl_2),disp=False)[0]

    # now interpolate between end states and write images
    # into subdirectories labeled with the VASP neb formalism
	ts = np.linspace(t_end1,t_end2,n_images)
	i = 0
	for t in ts:
		atoms[interp_ind].x = b1+t*a1
		atoms[interp_ind].y = b2+t*a2
		atoms[interp_ind].z = b3+t*a3
		dir = '%02d'%i
		if not os.path.exists(dir):
			os.mkdir(dir)
		atoms.write('%s/init.traj'%dir)
		i += 1

def handle_restart(dimer=False):
    """Short script that prepares a directory for a restarted job -- makes a
    backup directory and copies contents from the previous submission, etc.

    Use this at the top of a driver script to catch the case the in-loop
    backups cannot: a job killed by walltime part way through a VASP run,
    which never reaches the end-of-iteration snapshot in const_U_dimer.

    Optional arguments:
    dimer -- if True, also carry over the VTST dimer files (CENTCAR,
             DIMCAR, MODECAR, NEWMODECAR) into the backup, and restart
             from CENTCAR/NEWMODECAR rather than CONTCAR. Default False
             keeps the plain-relaxation behavior.
    """
    import os
    if not os.path.exists('./OUTCAR'):
        return False # job has not run here

    backup_path = make_backup(_DIMER_BACKUP_FILES if dimer else _BACKUP_FILES)

    if dimer and os.path.isfile('CENTCAR'):
        # CENTCAR is the dimer center, and NEWMODECAR the mode it ended on
        os.system('cp CENTCAR POSCAR')
        if os.path.isfile('NEWMODECAR'):
            os.system('cp MODECAR oldMODECAR')
            os.system('cp NEWMODECAR MODECAR')
    else:
        os.system('cp CONTCAR POSCAR')
    if os.path.isfile('WAVECAR'):
        os.remove('WAVECAR')
    return backup_path

def get_nearest_neighbors(atoms, scale=1.1):
    """Returns a dict mapping each atom index to a list of its nearest neighbor indices.

    Neighbor criterion: distance < scale * (r_i + r_j) where r_i and r_j are
    the covalent radii of each respective atom.
    """
    from ase.data import covalent_radii
    import numpy as np

    # Get per-atom radii as a 1D array, shape (N,)
    atomic_radii = np.array([covalent_radii[z] for z in atoms.get_atomic_numbers()])
    
    positions = atoms.positions  # shape (N, 3)
    N = len(atoms)

    # Pairwise displacement vectors, shape (N, N, 3)
    delta = positions[:, np.newaxis, :] - positions[np.newaxis, :, :]

    # Pairwise distances, shape (N, N)
    dists = np.sqrt(np.einsum('ijk,ijk->ij', delta, delta))

    # Pairwise cutoffs: scale * (r_i + r_j), shape (N, N)
    cutoffs = scale * (atomic_radii[:, np.newaxis] + atomic_radii[np.newaxis, :])

    # Boolean neighbor mask — exclude self (diagonal) explicitly
    neighbor_mask = (dists < cutoffs) & (dists > 0)

    # Build neighbor list dict
    neighbor_dict = {i: list(np.where(neighbor_mask[i])[0]) for i in range(N)}

    return neighbor_dict

