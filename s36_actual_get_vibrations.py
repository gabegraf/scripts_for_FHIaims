#!/usr/bin/env python3
"""
  modified by Simon Erker

  modified by Alaa Akkoush- Calculation of vibrational modes for
  periodic systems and vibrational Raman/IR intensities for periodic and non-
  periodic systems.

  modified by Connor Box - tidied up, removed sklearn requirement.
    I also removed the ability to submit calculations within the script, or call aims.
    Quite a lot of this functionality was significantly slowing down the processing of the results.

  This script is currently in the beta stage. 
  Please use it with caution. 
  It is designed for analyzing the vibrational modes of molecular adsorbates on large surfaces.
  Typically, only a small number of atoms need to be displaced, and the computational workload can be significant, potentially requiring the use of multiple CPUs.

  Plotting the IR and Raman-spectrum requires matplotlib.

  Example use:

      To generate finite difference displaced structures run:
          python get_vibrations.py name 0
      Run the calculations
      Then process the results with:
          python get_vibrations.py -Ipm --xmin 500 --xmax 3400 name 2

"""
from scipy.constants import value
# from numpy import array, reshape, zeros, ones, linalg, identity, float64, append, sqrt, arange, newaxis, delete, sum
import copy, os, shutil
# from sklearn import preprocessing
import numpy as np
import time
import sys

USAGE = """%prog [options] <name> <mode>

<name> will be used as prefix for all output files.
<mode> is an integer that determines the calculation mode:
  0: Generate finite-difference displaced structures
  1: Calculate vibrational modes
  2: Calculate vibrational modes and IR/Raman intensities
"""

def element(elem):
    """Element symbols"""
    data = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
            "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
            "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
            "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
            "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
            "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
            "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
            "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
            "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
            "Pa", "U"]
    try:
        return data[elem-1]
    except TypeError:
        return data.index(elem.capitalize())+1

def mass(elem):
    """Atomic masses"""
    data = [1.00794, 4.002602, 6.941, 9.012182, 10.811, 12.0107, 14.0067,
            15.9994, 18.9984032, 20.1797, 22.9897, 24.305, 26.9815, 28.0855,
            30.973762, 32.065, 35.453, 39.948, 39.0983, 40.078, 44.9559,
            47.867, 50.9415, 51.9961, 54.938, 55.845, 58.9332, 58.6934,
            63.546, 65.409, 69.723, 72.64, 74.9216, 78.96, 79.904, 83.8,
            85.4678, 87.62, 88.90585, 91.224, 92.9064, 95.94, 98, 101.07,
            102.9055, 106.42, 107.8682, 112.411, 114.818, 118.71, 121.76,
            127.6, 126.90447, 131.293, 132.9055, 137.327, 138.9055, 140.116,
            140.9077, 144.24, 145, 150.36, 151.964, 157.25, 158.9253, 162.5,
            164.9303, 167.259, 168.9342, 173.04, 174.967, 178.49, 180.9479,
            183.84, 186.207, 190.23, 192.217, 195.078, 196.9665, 200.59,
            204.3833, 207.2, 208.9804, 209, 210, 222, 223, 226,227, 232.0381,
            231.0359, 238.02891]
    try:
        return data[elem-1]
    except TypeError:
        return data[element(elem)-1]

class Atom:
    """Atom object"""
    def __init__(self, kind, coord, constraint=None, force=None):
        self.kind = kind
        self.coord = np.array(coord,np.float64)
        self.constraint = constraint if constraint else False
        self.force = np.array(force,np.float64) if force else np.zeros(3)
    def Z(self):
        return element(self.kind)
    def mass(self):
        return mass(self.Z())
    def copy(self):
        return copy.deepcopy(self)
    def to_str(self, elem_symb=True):
        elem = (element if elem_symb else str)(self.Z())
        if self.constraint:
          constraint = ' constrain_relaxation .true.\n'
        else:
          constraint = ''
        line='atom{0:14.8f}'.format(self.coord[0])
        line=line+'{0:14.8f}'.format(self.coord[1])
        line=line+'{0:14.8f}'.format(self.coord[2])
        line=line+'{0:>4}\n'.format(self.kind)
        line=line+constraint
        return line
class structure:
    """Structure of atoms object"""
    def __init__(self, atoms=None, vacuum_level=None, desc=None,\
                 periodic=None):
        self.atoms = atoms if atoms else np.array([])
        self.desc = desc if desc else ""
        self.vacuum_level = vacuum_level
        self.constrained=np.array([],bool)
        self.lattice_vector=np.zeros([0,3])
        self.periodic = periodic if periodic else False

    @classmethod
    def copy(self):
        return copy.deepcopy(self)
    def n(self):
        return len(self.atoms)
    def mass(self):
        return np.sum([mass(a.Z) for a in self.atoms])
    def join(self, add):
        self.atoms=np.append(self.atoms,add.copy())
        self.constrained=np.append(self.constrained,add.constraint)
    def to_str(self):
      line=''
      if self.periodic:
       # line=line+'set_vacuum_level{0:14.8f}\n'.format(self.vacuum_level)
        for i in range(3):
          line=line+'lattice_vector{0:14.8f}'.format(self.lattice_vector[i,0])
          line=line+'{0:14.8f}'.format(self.lattice_vector[i,1])
          line=line+'{0:14.8f}\n'.format(self.lattice_vector[i,2])
      for j in range(self.n()):
        line=line+self.atoms[j].to_str()
      return line

def split_line(lines):
  """Split input line"""
  line_array=np.array(lines.strip().split(' '))
  line_vals=line_array[line_array!='']
  return line_vals

def replace_submission(template_job, name, counter, filename):
    """
     Prepare submission script
     Only jobname and output file will be replaced
     Feel free to add more
    """
    template_job=template_job.replace('<jobname>',name+"_"+str(counter))
    template_job=template_job.replace('<outfile>',filename)
    job_out=open('job.sh','w')
    job_out.write(template_job)
    job_out.close()

def lorentz(g,x,x0):
    """
      Lorentzian function
      Parameters:
        - g  : Broadening value, gamma
        - x  : Frequency array upon which to calculate the lorentzian
        - x0 : Central frequency of lorentzian
    """
    lr = g/(np.pi*((x-x0)**2+g**2))

    return lr

def plot_spectra(x, intensity, freq, xlims, broadening=20, filename="spectrum.pdf"):
    """
    Plot IR and Raman spectra.

    Parameters:
        x (array-like): Array containing frequency values.
        infrared_intensity (array-like): Array containing IR intensity values.
        raman_intensity (array-like): Array containing Raman intensity values.
        freq (array-like): Array containing frequency values.
        filename (str): Name used for saving files.
        broadening (float): Width of the Lorentzian for convolution.

    Returns:
        None
    """
    # import matplotlib as mpl

    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.use('Agg')

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    spec_conv=np.zeros(len(x))

    for i in range(len(freq)):
        spec_conv += intensity[i] * lorentz(broadening, x, freq[i])

    ax.plot(x,spec_conv,'r',lw=2)

    ax.set_xlim(xlims)
    ax.set_ylim([-0.005,2000])
    ax.set_yticks([])
    ax.set_xlabel('Frequency [1/cm]',size=20)
    ax.set_ylabel('Intensity [a.u.]',size=20)
    fig.savefig(filename)
    
    np.savetxt(filename.replace("pdf","csv"), np.column_stack([x, spec_conv]), fmt=["%.15f",'%.15f'])

    return

def read_aims_output(folder, aimsoutname, struc, n_atoms):
    volume = -9999
    pol_jump = -9999
    dip_jump = -9999
    polr_jump = -9999
    forces_reached = False
    atom_count = 0
    
    with open(folder+'/'+aimsoutname, "r") as data:
        
        for line in data:
            if '| Unit cell volume ' in line:
                volume = np.float64(line.split()[-2])
                continue

            if 'Polarizability' in line:
                pol_jump = np.float64(line.split()[-6:])  # Periodic/cluster
                continue

            if '| Total dipole moment [eAng]' in line:
                dip_jump = np.float64(line.split()[-3:])  # Cluster
                continue

            if 'Cartesian Polarization' in line:
                polr_jump = np.float64(line.split()[-3:])  # Periodic
                continue

            if forces_reached and atom_count < n_atoms and 'Total atomic forces' not in line:
                struc.atoms[atom_count].force = np.float64(line.split()[2:])
                atom_count += 1
                if atom_count == n_atoms:
                    forces_reached = False
            if 'Total atomic forces' in line:
                forces_reached = True

    return volume, pol_jump, dip_jump, polr_jump

def get_phonopy_hessian():
    import yaml
    import numpy as np

    data = yaml.load(open("qpoints.yaml"),Loader=yaml.SafeLoader)
    dynmat = []
    dynmat_data = data['phonon'][0]['dynamical_matrix']
    for row in dynmat_data:
        vals = np.reshape(row, (-1, 2))
        dynmat.append(vals[:, 0] + vals[:, 1] * 1j)
    dynmat = np.array(dynmat)

    return dynmat

def main():
  import optparse
  # from numpy import sum

  # Parse command line
  parser = optparse.OptionParser(usage=USAGE)
  parser.add_option("-p", "--plot", action="store_true",
                    help="Generate pdf with IR-spectrum or Raman-spectrum, broadened with Lorentzian")
  parser.add_option("-I", "--IR", action="store_true",
                    help="Calculating infrared intensity")
  parser.add_option("-R", "--Raman", action="store_true",
                    help="Calculating Raman Intensity")
  parser.add_option("-M", "--IRRaman", action="store_true",
                    help="Calculating IR & Raman Intensity")
  parser.add_option("-i", "--info", action="store_true",
                      help="Set up/ Calculate vibrations & quit")
  # parser.add_option("-s", "--suffix", action="store",
  #                   help="Call suffix for binary e.g. 'mpirun -n 4 '",
  #                   default='')
  # parser.add_option("-r", "--run", action="store",
  #                   help="path to FHI-aims binary",default='')
  # parser.add_option("-x", "--relax", action="store_true",
  #                   help="Relax initial geometry")
  parser.add_option("-m", "--molden", action="store_true",
                    help="Output in molden format")
  parser.add_option("-w", "--distort", action="store_true",
                    help="Output geometry distorted along imaginary modes")
#   parser.add_option("-t", "--submit", action="store",
#                     help="""\
# Path to submission script, string <jobname>
# will be replaced by name + counter, string
#                             <outfile> will be replaced by filename""")
  parser.add_option("-d", "--delta", action="store", type="float",
                    help="Displacement", default=0.0025)
  parser.add_option("-b", "--broadening", action="store", type="float",
                    help="Broadening for IR and Raman spectrum in cm^{-1}", default=20)
  # parser.add_option("--kx", action="store", type="int", nargs=3, dest="gridx")
  # parser.add_option("--ky", action="store", type="int", nargs=3, dest="gridy")
  # parser.add_option("--kz", action="store", type="int", nargs=3, dest="gridz")

  parser.add_option("--xmin", action="store", type="float", default=-9999)
  parser.add_option("--xmax", action="store", type="float", default=-9999) 

  options, args = parser.parse_args()
  if options.info:
      print(__doc__)
      sys.exit(0)
  if len(args) != 2:
      parser.error("Need exactly two arguments")

  # AIMS_CALL=options.suffix+' '+options.run
  hessian_thresh = -1
  name=args[0]
  mode=args[1]
  # nx = options.gridx
  # ny = options.gridy
  # nz = options.gridz

  delta=options.delta
  broadening=options.broadening
  if mode=='2' and not options.Raman and not options.IR and not options.IRRaman:
      parser.error("Need an additional argument -I(only infrared), -R(only raman) or -M(raman and infrared) for vibrational calculations to run, check -h for more info")

  # run_aims=False
  # if options.run!='': run_aims=True

  # submit_script = options.submit is not None

  if options.plot or  mode=='2' or mode=='1':
      from pylab import transpose, eigh, argsort, sort, sign, dot



  # Constants
  bohr = value('Bohr radius') * 1e10  # Convert Bohr radius to Angstrom
  at_u = value('atomic mass unit-kilogram relationship')
  eV = value('electron volt-joule relationship')
  c = value('speed of light in vacuum')
  Ang = 1.0e-10
  hbar = value('reduced Planck constant')
  hessian_factor = eV / (at_u * Ang * Ang)
  grad_dipole_factor = (eV / (1. / (10 * c))) / Ang  # e -> D/Ang
  polr_factor = 0.062415091  # C/m2 to e/A2
  ir_factor = 1

  # File names
  atomicmasses = f'masses.{name}.dat'
  xyzfile = f'{name}.xyz'
  moldenname = f'{name}.molden'
  hessianname = f'hessian.{name}.dat'
  mass_vec = f'mass_vec.{name}.dat'
  cartesian_eig_vec = f'car_eig_vec.{name}.dat'
  normalized = f'massweighted_Hessian.{name}.dat'
  eigenvectors = f'eigen_vectors.{name}.dat'
  normalmodes = f'normalmodes.{name}.dat'
  graddipolename = f'grad_dipole.{name}.dat'
  irname = f'{name}.ir'
  Ramanname = f'{name}.Raman'
  aimsoutname = "aims.out"

  # Define deltas and coefficients
  deltas=np.array([-delta,delta])
  coeff=np.array([-1,1])
  c_zero = - 1. / (2. * delta)

  # Read templates
  with open('control.in', 'r') as f:
      template_control = f.read()

  # if submit_script:
  #     with open(options.submit, 'r') as f:
  #         template_job = f.read()


  # Central Point
  folder=''                                  # Dummy
  # if options.relax and (mode == '0' or mode == '2' or mode == '1'):
  #     filename = name + '.out'
  #     folder = name + '_relaxation'
  #     if not os.path.exists(folder):
  #         os.mkdir(folder)
  #     shutil.copy('geometry.in', folder + '/geometry.in')
  #     new_control = open(folder + '/control.in', 'w')
  #     new_control.write(template_control + 'relax_geometry trm 1E-3\n')
  #     new_control.close()
  #     os.chdir(folder)
  #     print('Central Point')
  #     # if run_aims:
  #     #     os.system(AIMS_CALL + ' > ' + filename)
  #     # if submit_script:
  #     #     replace_submission(template_job, name, 0, filename)
  #     os.chdir('..')



  # Check for relaxed geometry
  if os.path.isfile('geometry.in.next_step'):
      geometry = open('geometry.in.next_step', 'r')
  else:
      geometry = open('geometry.in', 'r')


  # Read input geometry
  n_line = 0
  struc = structure()
  lines = geometry.readlines()
  m = 0
  rec_lat = np.zeros(shape=(3, 3))
  for line in lines:
      n_line = n_line + 1
      if line.rfind('set_vacuum_level') != -1:
          struc.vacuum_level = float(split_line(line)[-1])
      if line.rfind('lattice_vector') != -1:
          lat = split_line(line)[1:]
          struc.lattice_vector = np.append(struc.lattice_vector, np.float64(np.array(lat))[np.newaxis, :], axis=0)
          struc.periodic = True
          rec_lat[m] = line.split()[1:4]
          m = m + 1
      if line.rfind('atom') != -1:
          line_vals = split_line(line)
          at = Atom(line_vals[-1], line_vals[1:-1])
          if n_line < len(lines):
              nextline = lines[n_line]
              if nextline.rfind('constrain_relaxation') != -1:
                  at = Atom(line_vals[-1], line_vals[1:-1], True)
              else:
                  at = Atom(line_vals[-1], line_vals[1:-1])
          struc.join(at)
      if line.rfind('atom_frac') != -1:
          parser.error("Use atom not atom_frac")
  geometry.close()

  n_atoms = struc.n()
  n_constrained = n_atoms - np.sum(struc.constrained)
  print('Total number of atoms found: ' + str(n_atoms))
  print('Number of atoms considered in frequency calculation: ' + str(n_constrained))

  if struc.periodic:
      rec_lat = np.linalg.pinv(rec_lat).transpose()
      # X_normalized = preprocessing.normalize(rec_lat, norm="l2")

      norms = np.linalg.norm(rec_lat, axis=1)
      X_normalized = rec_lat / norms[:, np.newaxis]

      x = np.linalg.solve(X_normalized, np.array([1, 0, 0]))
      y = np.linalg.solve(X_normalized, np.array([0, 1, 0]))
      z = np.linalg.solve(X_normalized, np.array([0, 0, 1]))
      # R = np.array([x, y, z])


  # Atomic mass file
  mass_file = open(atomicmasses, 'w')
  mass_vector = np.zeros([0])
  mass_vector_inv = np.zeros([0])
  for at_unconstrained in struc.atoms[struc.constrained == False]:
      mass_vector = np.append(mass_vector, np.ones(3) * 1 / np.sqrt(at_unconstrained.mass()))
      mass_vector_inv = np.append(mass_vector_inv, np.ones(3) * np.sqrt(at_unconstrained.mass()))
      line = '{0:10.5f}'.format(at_unconstrained.mass())
      for i in range(3):
          line = line + '{0:11.4f}'.format(at_unconstrained.coord[i])
      line = line + '{0:}\n'.format(at_unconstrained.kind)
      mass_file.writelines(line)
  mass_file.close()

  # Init
  polr = np.zeros([n_constrained * 3, 3])
  dip = np.zeros([n_constrained*3,3])
  pol = np.zeros([n_constrained*3,6])
  dip_jump=np.zeros([])
  pol_jump=np.zeros([])
  polr_jump=np.zeros([])
  hessian = np.zeros([n_constrained*3,n_constrained*3])
  index=0
  counter=1

  # if mode=='2' and  options.IR and struc.periodic and not (options.gridx or  options.gridy or  options.gridz):
  #     parser.error("Need to determine the polarization intergration grid")
  # if mode=='2' and options.IRRaman and struc.periodic and not (options.gridx or  options.gridy or  options.gridz):
  #     parser.error("Need to determine the polarization intergration grid")
  # Set up / Read folders for displaced atoms
  for atom in np.arange(n_atoms)[struc.constrained==False]:
    for coord in np.arange(3):
      for delta in deltas:
        #filename=name+'.i_atom_'+str(atom)+'.i_coord_'+str(coord)+'.displ_'+str(delta)+'.out'
        folder=name+'.i_atom_'+str(atom)+'.i_coord_'+str(coord)+'.displ_'+str(delta)

        struc_new = copy.deepcopy(struc)
        struc_new.atoms[atom].coord[coord] += delta
        geoname = 'geometry.i_atom_' + str(atom) + '.i_coord_' + str(coord) + '.displ_' + str(delta) + '.in'
        if not os.path.exists(folder):
            os.mkdir(folder)
        with open(folder + '/geometry.in', 'w') as new_geo:
            newline = '#\n# temporary structure-file for finite-difference calculation of forces\n'
            newline += '# displacement {0:8.4f} of atom '.format(delta) + '{0:5} direction {1:5}\n#\n'.format(atom, coord)
            new_geo.writelines(newline + struc_new.to_str())
        with open(folder + '/control.in', 'w') as new_control:
            template_control = template_control.replace('relax_geometry', '#relax_geometry')
            new_control.write(template_control)

        # if  mode=='2':   # Put new geometry and control.in into folder
        #   with open(folder + '/control.in', 'a') as new_control:
        #     if struc.periodic:
              
        #       if  options.IRRaman:
        #           new_control.write(template_control+'compute_forces .true. \n'+'final_forces_cleaned '+'.true. \n'+'DFPT dielectric\n'+  "KS_method serial \n"
        #                           + "output polarization    "
        #                           + str(1)
        #                           + " {} {} {}\n".format(nx[0], nx[1], nx[2])
        #                           + "output polarization    "
        #                           + str(2)
        #                           + " {} {} {}\n".format(ny[0], ny[1], ny[2])
        #                           + "output polarization    "
        #                           + str(3)
        #                           + " {} {} {}\n".format(nz[0], nz[1], nz[2]))
        #       elif options.IR:
        #           new_control.write(template_control+'compute_forces .true. \n'+'final_forces_cleaned '+'.true. \n'+ "KS_method serial \n"
        #                           + "output polarization    "
        #                           + str(1)
        #                           + " {} {} {}\n".format(nx[0], nx[1], nx[2])
        #                           + "output polarization    "
        #                           + str(2)
        #                           + " {} {} {}\n".format(ny[0], ny[1], ny[2])
        #                           + "output polarization    "
        #                           + str(3)
        #                           + " {} {} {}\n".format(nz[0], nz[1], nz[2]))
        #       elif options.Raman:
        #           new_control.write(template_control+'compute_forces .true. \n'+'final_forces_cleaned '+'.true. \n'+'DFPT dielectric\n')
        #     else:
        #       if  options.IRRaman:
        #           new_control.write(template_control+'compute_forces .true. \n'+'final_forces_cleaned '+'.true. \n'+'output dipole \n'+'DFPT polarizability\n')
        #       elif options.IR:
        #           new_control.write(template_control+'compute_forces .true. \n'+'final_forces_cleaned '+'.true. \n'+'output dipole \n')
        #       elif options.Raman:
        #           new_control.write(template_control+'compute_forces .true. \n'+'final_forces_cleaned '+'.true. \n'+'DFPT polarizability\n')

        # if mode=='1':   # Put new geometry and control.in into folder
        #   with open(folder + '/control.in', 'a') as new_control:
        #     template_control=template_control.replace('relax_geometry','#relax_geometry')
        #     new_control.write(template_control+'compute_forces .true. \n'+'final_forces_cleaned '+'.true. \n')

        # Run aims or write submit script
        # if mode in ['1','2']:
        #   os.chdir(folder)                                   # Change directoy
        #   print('Processing atom: '+str(atom+1)+'/'+str(n_atoms)+', coord.: '+str(coord+1)+'/'+str(3)+', delta: '+str(delta))
        #   if run_aims:
        #     os.system(AIMS_CALL+' > '+filename) # Run aims and pipe the output into a file named 'filename'
        #   if submit_script: replace_submission(template_job, name, counter, filename)
        #   #os.system('sbatch -N 4 job.sh') # Mind the environment variables <======================================== Modify according to you system
        #   counter=counter+1
        #   os.chdir('..')
        #   time.sleep(2)

        if mode in ['1','2']:
          volume, pol_jump, dip_jump, polr_jump = read_aims_output(folder, aimsoutname, struc, n_atoms)


          #forces=np.array([])
          #for at_unconstrained in struc.atoms[struc.constrained==False]:
           #   forces=np.append(forces,coeff[deltas==delta]*at_unconstrained.force)
          #hessian[index,:]=hessian[index,:]+forces*c_zero

        if  mode=='2':   # Read output
          if struc.periodic:
            pol[index,:]=pol[index,:]+pol_jump*coeff[deltas==delta]*c_zero #grad polar using finite difference
            polr[index, :] = (polr[index, :] + polr_jump * coeff[deltas == delta] * c_zero) #C/m2
          else:
            dip[index,:]=dip[index,:]+dip_jump*coeff[deltas==delta]*c_zero
            pol[index,:]=pol[index,:]+pol_jump*coeff[deltas==delta]*c_zero#grad polar using finite difference
      index=index+1


  if mode in ['1','2']:
    print('Entering hessian diagonalization')
    print('Number of atoms                = '+str(n_atoms))
    print('Name of Hessian input file     = '+hessianname)

    print('Name of Masses  input file     = '+atomicmasses)
    print('Name of XYZ output file        = '+xyzfile)
    print('Threshold for Matrix elements  = '+str(hessian_thresh))
    #if (hessian_thresh < 0.0): print('     All matrix elements are taken into account by default\n')
    #np.savetxt(hessianname,hessian)

    mass_mat=mass_vector[:,np.newaxis]*mass_vector[np.newaxis,:]  #1/sqrt(mimj)
    # mass_mat_inv=mass_vector_inv[:,np.newaxis]*mass_vector_inv[np.newaxis,:]  # sqrt(mimj)
    #hessian[abs(hessian)<hessian_thresh]=0.0

    # Mass weighted Hessian (Dynamical matrix)
    # Dij=Hij/sqrt(mimj)
    hessian=get_phonopy_hessian()
    hessian=hessian*hessian_factor #  unit was ev/amu*Ang*Ang changed to SI
    #hessian=(hessian+transpose(hessian))/2.

    # Diagonalize hessian (scipy)
    print('Solving eigenvalue system for Hessian Matrix')
    freq, eig_vec = eigh(hessian)

    print('Eigen vectors input file       = '+eigenvectors)
    np.savetxt(normalized,hessian)


    eig_vec=eig_vec[:,argsort(freq)]  # sort in ascending order in colums
    np.savetxt(eigenvectors,eig_vec)  # Qij
    # Finding Cartesian eigen vectors X=M^-1/2Q
    eig_vec_car = eig_vec*mass_vector[:,np.newaxis]*np.ones(len(mass_vector))[np.newaxis,:]
    np.savetxt(cartesian_eig_vec, transpose(eig_vec_car))
    # saving Diagonal mass matrix of 3N size
    M=mass_vector[:,np.newaxis]*np.identity(len(mass_vector))[np.newaxis,:]
    np.savetxt(mass_vec,M.reshape((n_constrained*3,-1)), fmt="%s")
    # From hessian_factor freq is in SI units
    freq=sort(sign(freq)*np.sqrt(abs(freq)))
    #ZPE (eV) half the vibrational frequency
    ZPE=hbar*(freq)/(2.0*eV)
    # Finding frrquency in cm-1:
    freq=freq/(200.*np.pi*c)


    # saving normal modes
    newline_eig=''
    # import numpy as np
    for x in np.arange(len(freq)):
          newline_eig=newline_eig+'Atoms: {0:2}  Mode #: {1:2} Freq: {2:11} cm-1,Mass weighted displacment dq is: \n {3:}\n'.format(n_constrained,x+1,freq[x],np.array2string(eig_vec[:, x], max_line_width=10000))
    eig=open(normalmodes,'w')
    eig.writelines(newline_eig)
    eig.close()

    eig_vec = eig_vec*mass_vector[:,np.newaxis]*np.ones(len(mass_vector))[np.newaxis,:]
    reduced_mass=np.sum(eig_vec**2,axis=0)




    if  mode=='2': # Calculate vibrations

      print('Name of grad dipole input file = '+graddipolename)
      np.savetxt(graddipolename,dip)

      grad_dipole = dip * grad_dipole_factor

      if options.IR or options.IRRaman:
        if struc.periodic:
      #Calculation of Infrared

            grad_polr= polr*volume *polr_factor* grad_dipole_factor # D/Ang
            infrared_intensity = np.sum(dot(transpose(grad_polr),eig_vec)**2,axis=0)*ir_factor  #D^2/A^2*amu
        else:
      # norm = sqrt(reduced_mass)
            infrared_intensity = np.sum(dot(transpose(grad_dipole),eig_vec)**2,axis=0)*ir_factor  #D^2/A^2*amu
      # eig_vec = eig_vec/norm
      if options.Raman or options.IRRaman:
        alphas=dot(transpose(pol),eig_vec)

      # mean polarizabilty derivative
        alphasxx=np.zeros(n_constrained*3)

        alphasyy=np.zeros(n_constrained*3)
        alphaszz=np.zeros(n_constrained*3)
        alphasxy=np.zeros(n_constrained*3)
        alphasxz=np.zeros(n_constrained*3)
        alphasyz=np.zeros(n_constrained*3)



        for ii in range(n_constrained*3):
            alphasxx[ii]=  alphas[0,ii]
            alphasyy[ii]=  alphas[1,ii]
            alphaszz[ii]=  alphas[2,ii]
            alphasxy[ii]=  alphas[3,ii]
            alphasxz[ii]=  alphas[4,ii]
            alphasyz[ii]=  alphas[5,ii]
        alpha= (alphasxx + alphasyy + alphaszz)*(1./3)
        # polarizability tensor derivative
        beta=(alphasxx-alphasyy)**2+(alphasxx-alphaszz)**2+(alphasyy-alphaszz)**2+6*((alphasxy)**2+(alphasxz)**2+(alphasyz)**2)
      #Raman Scattering Intensity:
        raman_intensity=45*(alpha**2)+(7./2)*(beta)
        raman_intensity=raman_intensity*0.02195865620442408 #bohr^6/ang^2 to ang^4


    # The rest is output, xyz, IR,...
    print('Results\n')
    print('List of all frequencies found:')
    if options.IR:
           print('Mode number      Frequency [cm^(-1)]   Zero point energy [eV]   IR-intensity [D^2/Ang^2amu]')
    elif options.Raman:
           print('Mode number      Frequency [cm^(-1)]   Zero point energy [eV]       Raman-intensity [Ang^4/amu] ')
    elif options.IRRaman:
           print('Mode number      Frequency [cm^(-1)]   Zero point energy [eV]   IR-intensity [D^2/Ang^2amu]     Raman-intensity [Ang^4/amu] ')

    for i in range(len(freq)):
        if options.IR:
           print('{0:11}{1:25.8f}{2:25.8f}{3:25.8f}'.format(i+1,freq[i],ZPE[i],infrared_intensity[i]))

        elif options.Raman:
           print('{0:11}{1:25.8f}{2:25.8f}{3:25.8f}'.format(i+1,freq[i],ZPE[i],raman_intensity[i]))

        elif options.IRRaman:
           print('{0:11}{1:25.8f}{2:25.8f}{3:25.8f}{4:25.8f}'.format(i+1,freq[i],ZPE[i],infrared_intensity[i],raman_intensity[i]))

    print('\n')
    print('Summary of zero point energy for entire system:')
    print('| Cumulative ZPE               = {0:15.8f} eV'.format(np.sum(ZPE)))
    if struc.periodic:
         print('| without first three eigenmodes = {0:15.8f} eV\n'.format(np.sum(ZPE)-np.sum(ZPE[:3])))
         print('Stability checking - eigenvalues should all be positive for a stable structure. ')
         print('The three smallest frequencies should be (almost) zero:')
         string=''
         for zz in freq[:3]: string=string+'{0:25.8f}'.format(zz)
         print(string)

    else:
         print('| without first six eigenmodes = {0:15.8f} eV\n'.format(np.sum(ZPE)-np.sum(ZPE[:6])))
         print('Stability checking - eigenvalues should all be positive for a stable structure. ')
         print('The six smallest frequencies should be (almost) zero:')
         string=''
         for zz in freq[:6]: string=string+'{0:25.8f}'.format(zz)
         print(string)

    print('Compare this with the largest eigenvalue, ')
    print('{0:25.8f}'.format(freq[-1]))

    # Molden stuff
    nums=np.arange(n_atoms)[struc.constrained==False]
    nums2=np.arange(n_atoms)[struc.constrained]
    newline=''
    newline_ir='[INT]\n'
    newline_Raman='[INT]\n'
    if options.molden:
      newline_molden='[Molden Format]\n[GEOMETRIES] XYZ\n'
      newline_molden=newline_molden+'{0:6}\n'.format(n_atoms)+'\n'
      for i_atoms in range(n_constrained):
        newline_molden=newline_molden+'{0:6}'.format(struc.atoms[nums[i_atoms]].kind)
        for i_coord in range(3):
          newline_molden=newline_molden+'{0:10.4f}'.format(struc.atoms[nums[i_atoms]].coord[i_coord])
        newline_molden=newline_molden+'\n'
      newline_molden=newline_molden+'[FREQ]\n'
      for i in range(len(freq)):
        newline_molden=newline_molden+'{0:10.3f}\n'.format(freq[i])
      newline_molden=newline_molden+'[INT]\n'
      if mode=='2':
        if options.IR or options.IRRaman:
          for i in range(len(freq)):
            newline_molden=newline_molden+'{0:17.6e}\n'.format(infrared_intensity[i])
          newline_molden=newline_molden+'[FR-COORD]\n'
          newline_molden=newline_molden+'{0:6}\n'.format(n_atoms)+'\n'
      for i_atoms in range(n_constrained):
        newline_molden=newline_molden+'{0:6}'.format(struc.atoms[nums[i_atoms]].kind)
        for i_coord in range(3):
          newline_molden=newline_molden+'{0:10.4f}'.format(struc.atoms[nums[i_atoms]].coord[i_coord]/bohr)
        newline_molden=newline_molden+'\n'
      newline_molden=newline_molden+'[FR-NORM-COORD]\n'



    for i in range(len(freq)):

      # Distoring the structure along the eigenmodes
      newline=newline+'{0:6}\n'.format(n_atoms)
      if freq[i]>0:
        newline=newline+'stable frequency at '
      elif freq[i]<0:
        newline=newline+'unstable frequency at '
        if options.distort and freq[i]<-50:
          struc_new=copy.deepcopy(struc)
          for i_atoms in range(n_constrained):
            for i_coord in range(3):
              struc_new.atoms[i_atoms].coord[i_coord]=\
              struc_new.atoms[i_atoms].coord[i_coord]+eig_vec[(i_atoms)*3+i_coord,i]
          geoname=name+'.distorted.vibration_'+str(i+1)+'.geometry.in'
          new_geo=open(geoname,'w')
          newline_geo='#\n# distorted structure-file for based on eigenmodes\n'
          newline_geo=newline_geo+'# vibration {0:5} :{1:10.3f} 1/cm\n#\n'.format(i+1,freq[i])
          new_geo.writelines(newline_geo+struc_new.to_str())
          new_geo.close()
      elif freq[i]==0:
        newline=newline+'translation or rotation '
      newline=newline+'{0:10.3f} 1/cm IR int. is '.format(freq[i])

      #XYZ file
      if mode=='2':
        if options.IR or options.IRRaman:
          newline=newline+'{0:10.4e} D^2/Ang^2; red. mass is '.format(infrared_intensity[i])
      newline=newline+'{0:5.3f} a.m.u.; force const. is '.format(1.0/reduced_mass[i])
      newline=newline+'{0:5.3f} mDyne/Ang.\n'.format(((freq[i]*(200*np.pi*c))**2)*(1.0/reduced_mass[i])*at_u*1.e-2)
      if options.molden: newline_molden=newline_molden+'vibration {0:6}\n'.format(i+1)
      for i_atoms in range(n_constrained):
        newline=newline+'{0:6}'.format(struc.atoms[nums[i_atoms]].kind)
        for i_coord in range(3):
          newline=newline+'{0:10.4f}'.format(struc.atoms[nums[i_atoms]].coord[i_coord])
        for i_coord in range(3):
          newline=newline+'{0:10.4f}'.format(eig_vec[(i_atoms)*3+i_coord,i])
          if options.molden: newline_molden=newline_molden+'{0:10.4f}'.format(eig_vec[(i_atoms)*3+i_coord,i]/bohr)
        newline=newline+'\n'
        if options.molden: newline_molden=newline_molden+'\n'
      for i_atoms in range(n_atoms-n_constrained):
        newline=newline+'{0:6}'.format(struc.atoms[nums2[i_atoms]].kind)
        for i_coord in range(3):
          newline=newline+'{0:10.4f}'.format(struc.atoms[nums2[i_atoms]].coord[i_coord])
        for i_coord in range(3):
          newline=newline+'{0:10.4f}'.format(0.0)
        newline=newline+'\n'

      if mode=='2':
        if options.IR or options.IRRaman:
          newline_ir=newline_ir+'{0:11}{1:25.8f}{2:25.8f}{3:25.8f}\n'.format(i+1,freq[i],ZPE[i],infrared_intensity[i])

        if options.Raman or options.IRRaman:
          newline_Raman=newline_Raman+'{0:11}{1:25.8f}{2:25.8f}{3:25.8f}\n'.format(i+1,freq[i],ZPE[i],raman_intensity[i])
          
    xyz=open(xyzfile,'w')
    xyz.writelines(newline)
    xyz.close()


    if options.molden:
      molden=open(moldenname,'w')
      molden.writelines(newline_molden)
      molden.close()


    if  mode=='2': # Calculate vibrations

      if options.IR or options.IRRaman:
        ir=open(irname,'w')
        ir.writelines('#Mode number      Frequency [cm^(-1)]   Zero point energy [eV]   IR-intensity [D^2/Ang^2amu]')
        ir.writelines(newline_ir)
        ir.close()

      if options.Raman or options.IRRaman:
        Raman=open(Ramanname,'w')
        Raman.writelines('#Mode number      Frequency [cm^(-1)]   Zero point energy [eV]   Raman-intensity  [Ang^4/amu]')
        Raman.writelines(newline_Raman)
        Raman.close()

      if options.plot:
        x=np.linspace(freq.min()-500,freq.max()+500,1000)

        xlims = [freq.min()+90,freq.max()+50]
        if options.xmin!=-9999:
            xlims[0] = options.xmin
        if options.xmax!=-9999:
            xlims[1]=options.xmax
          
        if options.IR or options.IRRaman:
          plot_spectra(x,infrared_intensity,freq,xlims,broadening,name+'_IR_spectrum.pdf')

        if options.Raman or options.IRRaman:
            plot_spectra(x,raman_intensity,freq,xlims,broadening,name+'_Raman_spectrum.pdf')

if __name__ == "__main__":
    main()
    print('\n Done. ')
