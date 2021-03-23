import sys, os, re
import numpy as np
import scipy as sp
from numpy.matlib import *
#from numpy.matlib import sqrt
from scipy.linalg import *
#from scipy.linalg import
from numba import jit


__author__ = 'Sean M. Ryno'
__copyright__ = 'Copyright 2016, Sean M. Ryno'
__credits__ = ['Sean M. Ryno', 'J. J. Goings', 'Brad Rose']
__license__= 'GPL v3.0'
__version__ = '1.13'
__maintainer__ = 'Sean M. Ryno'
__email__ = 'sean.m.ryno@gmail.com'
__status__ = 'Development'


def TI(argv):

    if (len(argv) != 2) and (len(argv) != 3) and (len(argv) != 7):
        print("Usage:", argv[0], "system_prefix [g09/orca/qchem] optional([homo1] [lumo1] [homo2] [lumo2])", file=sys.stderr)
        print("\"", argv[0], "-h\" will print a help file.", file=sys.stderr)
        print("Exiting...", file=sys.stderr)
        sys.exit(1)

    if argv[1] == "-h":
        TI_header()
        sys.exit(0)
    else:
        pass

    if len(argv) == 3:
        mon1_homo = 0
        mon1_lumo = 0
        mon2_homo = 0
        mon2_lumo = 0
    elif len(argv) == 7:
        mon1_homo = int(argv[3])
        mon1_lumo = int(argv[4])
        mon2_homo = int(argv[5])
        mon2_lumo = int(argv[6])
    else:
        print("There is an error reading the input homo and lumos", file=sys.stderr)
        print("Exiting...", file=sys.stderr)
        sys.exit(1)

    if argv[2] == 'g09':
        TI_g09(argv[1], mon1_homo, mon1_lumo, mon2_homo, mon2_lumo)
    elif argv[2] == 'orca':
        TI_orca(argv[1], mon1_homo, mon1_lumo, mon2_homo, mon2_lumo)
    elif argv[2] == 'adf':
        TI_adf(argv[1], mon1_homo, mon1_lumo, mon2_homo, mon2_lumo)
    elif argv[2] == 'qchem':
        TI_qchem(argv[1], mon1_homo, mon1_lumo, mon2_homo, mon2_lumo)

    else:
        print("Unrecognized program type.", file=sys.stderr)
        print("Exiting...", file=sys.stderr)
        sys.exit(1)


def TI_header():
    # Prints the header for this program
    print("################################################################################", file=sys.stdout)
    print("################################################################################", file=sys.stdout)
    print("##                                                                            ##", file=sys.stdout)
    print("##                        Electronic Coupling Program                         ##", file=sys.stdout)
    print("##                                                                            ##", file=sys.stdout)
    print("##                               Sean M. Ryno                                 ##", file=sys.stdout)
    print("##                                                                            ##", file=sys.stdout)
    print("##                             January 24, 2018                               ##", file=sys.stdout)
    print("##                                                                            ##", file=sys.stdout)
    print("##                                                                            ##", file=sys.stdout)
    print("## Usage: transfer_integral.py system_name [g09/orca/qchem]                   ##", file=sys.stdout)
    print("##                             optional([homo1] [lumo1] [homo2] [lumo2])      ##", file=sys.stdout)
    print("##                                                                            ##", file=sys.stdout)
    print("## Here, homo[1/2], lumo[1/2] are the upper and lower bounds of orbital of    ##", file=sys.stdout)
    print("##       interest for monomer [1/2] using the orbital numbering of            ##", file=sys.stdout)
    print("##       monomer [1/2]. All possible combinations of couplings will be        ##", file=sys.stdout)
    print("##       calculated.                                                          ##", file=sys.stdout)
    print("##                                                                            ##", file=sys.stdout)
    print("## Note: This program makes use of components of the gauss_parse.py           ##", file=sys.stdout)
    print("##       program of J. J. Goings                                              ##", file=sys.stdout)
    print("##                                                                            ##", file=sys.stdout)
    print("## Note: The output files from calculations of dimer and monomers must        ##", file=sys.stdout)
    print("##       include molecular orbital coefficients and the dimer overlap matrix  ##", file=sys.stdout)
    print("##       in addition to standard output                                       ##", file=sys.stdout)
    print("##                                                                            ##", file=sys.stdout)
    print("## Note: Program auto-detects restricted or unrestricted calculations         ##", file=sys.stdout)
    print("##       and provide respective couplings                                     ##", file=sys.stdout)
    print("##                                                                            ##", file=sys.stdout)
    print("## Note: For combined Gaussian jobs, the dimer calculation must come first.   ##", file=sys.stdout)
    print("##                                                                            ##", file=sys.stdout)
    print("## Note: Gaussian must adopt the following naming convention:                 ##", file=sys.stdout)
    print("##       dimer: *_A_B.fchk & *_A_B.log                                        ##", file=sys.stdout)
    print("##       monomer 1: *_A.fchk                                                  ##", file=sys.stdout)
    print("##       monomer 2: *_B.fchk                                                  ##", file=sys.stdout)
    print("##                                                                            ##", file=sys.stdout)
    print("## Gaussian keywords: # b3lyp/6-31g** nosymm pop=full iop(3/33=1,3/32=2)      ##", file=sys.stdout)
    print("##                                                                            ##", file=sys.stdout)
    print("## Additional Orca keywords: print[ P_mos ] 1 & print[ P_overlap ] 1          ##", file=sys.stdout)
    print("##                                                                            ##", file=sys.stdout)
    print("## Q-Chem $rem section must include:                                          ##", file=sys.stdout)
    print("##       iprint            200                                                ##", file=sys.stdout)
    print("##       print_orbitals    1000000                                            ##", file=sys.stdout)
    print("##                                                                            ##", file=sys.stdout)
    print("## Based upon: E. F. Valeev, V. Coropceanu, D. A. da Silva Filho, S. Salman,  ##", file=sys.stdout)
    print("##             J. L. Bredas, J. Am. Chem. Soc. 128, 9882-9886 (2006)          ##", file=sys.stdout)
    print("##                                                                            ##", file=sys.stdout)
    print("################################################################################", file=sys.stdout)
    print("################################################################################", file=sys.stdout)


def TI_outfile_header(fout):
    #Prints the header of the outfile
    fout.write("################################################################################\n")
    fout.write("################################################################################\n")
    fout.write("##                                                                            ##\n")
    fout.write("##                        Electronic Coupling Program                         ##\n")
    fout.write("##                                                                            ##\n")
    fout.write("##                               Sean M. Ryno                                 ##\n")
    fout.write("##                                                                            ##\n")
    fout.write("##                             January 24, 2018                               ##\n")
    fout.write("##                                                                            ##\n")
    fout.write("## Note: This program makes use of components of the gauss_parse.py           ##\n")
    fout.write("##       program of J. J. Goings                                              ##\n")
    fout.write("##                                                                            ##\n")
    fout.write("## Based upon: E. F. Valeev, V. Coropceanu, D. A. da Silva Filho, S. Salman,  ##\n")
    fout.write("##             J. L. Bredas, J. Am. Chem. Soc. 128, 9882-9886 (2006)          ##\n")
    fout.write("##                                                                            ##\n")
    fout.write("################################################################################\n")
    fout.write("################################################################################\n")
    fout.write("\n")
    fout.write("\n")


def TI_adf(system_name, mon1_homo, mon1_lumo, mon2_homo, mon2_lumo):

    d_file = open(system_name+'_A_B.out', 'r')
    d_file_line = []
    for line in d_file:
        d_file_line.append(line)

    d_file.close()

    # Get dimer values
    nbf_d = adf_get_basis_functions(d_file_line)
    print("nbf_d", file=sys.stdout)
    e_d = adf_get_alpha_orbital_energies(d_file_line, nbf_d)
    c_d = adf_get_MO_coefficients(d_file_line, nbf_d)
    s_d = adf_get_overlap_matrix(d_file_line, nbf_d)
    f_d = adf_calculate_fock(c_d, s_d, nbf_d, e_d)

    d_file_line = None

    # Open monomer A out file
    mon1_file = open(system_name+'_A.out', 'r')
    mon1_file_line = []
    for line in mon1_file:
        mon1_file_line.append(line)

    mon1_file.close()

    # Get monomer A values
    nbf_mon1 = adf_get_basis_functions(mon1_file_line)
    if mon1_homo == 0:
        mon1_homo, mon1_lumo = adf_get_homo_lumo(mon1_file_line)
    else:
        pass
    c_mon1 = adf_get_MO_coefficients(mon1_file_line, nbf_mon1)

    mon1_file_line = None

    # Open monomer B out file
    mon2_file = open(system_name+'_B.out', 'r')
    mon2_file_line = []
    for line in mon2_file:
        mon2_file_line.append(line)

    mon2_file.close()

    # Get monomer B values
    nbf_mon2 = adf_get_basis_functions(mon2_file_line)
    if mon2_homo == 0:
        mon2_homo, mon2_lumo = adf_get_homo_lumo(mon2_file_line)
    else:
        pass
    c_mon2 = adf_get_MO_coefficients(mon2_file_line, nbf_mon2)

    mon2_file_line = None

    cfrag = adf_build_cfrag(nbf_d, nbf_mon1, nbf_mon2, c_mon1, c_mon2)

    adf_calculate_coupling(system_name, f_d, s_d, mon1_homo, mon1_lumo, mon2_homo, mon2_lumo, cfrag, nbf_mon1)


def adf_get_basis_functions(adf_output):
        # Get the number of basis functions
        for line in adf_output:
            if "Nr. of SFOs" in line:
                nbf = int(line.split()[4])
                break
            else:
                pass

        return nbf


def adf_get_alpha_orbital_energies(adf_output, nbf):
        # Get alpha orbital energies
        e = np.matlib.zeros(nbf)
        for line in adf_output:
            if "Orbital Energies, per Irrep and Spin" in line:
                line_index = int(adf_output.index(line)) + 5
                break
            else:
                pass
        count = 0
        while count < nbf:
            mo_line = adf_output[line_index + count].split()
            e[0, count] = float(mo_line[2])
            count += 1

        return e


def adf_get_MO_coefficients_old(qchem_output, nbf):
        # Get alpha MO coefficients
        c_temp = np.matlib.zeros(nbf * nbf)
        for line in qchem_output:
            if "RESTRICTED (RHF) MOLECULAR ORBITAL COEFFICIENTS" in line:
                line_index = int(qchem_output.index(line)) + 2
                break
            else:
                pass
        count = 0
        icount = 0
        while count < nbf * nbf:
            if ("eigenvalues:" in qchem_output[line_index + icount]) or (len(qchem_output[line_index + icount].split()) < 8):
                icount += 1
                pass
            else:
                mo_line = qchem_output[line_index + icount].split()
                icount += 1
                for i in range(len(mo_line[4:])):
                    c_temp[0, count] = float(mo_line[i+4])
                    count += 1
        c = c_temp.reshape((nbf, nbf))

        return c


def adf_get_MO_coefficients(qchem_output, nbf):
    c = np.matlib.zeros((nbf, nbf))
    for line in qchem_output:
        if "RESTRICTED (RHF) MOLECULAR ORBITAL COEFFICIENTS" in line:
            line_index = int(qchem_output.index(line)) + 2
            break
        else:
            pass
    acounter = 0
    bcounter = 0
    ccounter = 0
    count = 0
    icount = 0
    while count < nbf * nbf:
        if ("eigenvalues" in qchem_output[line_index + icount]):
            icount += 1
            acounter = 0
            pass
        elif (len(qchem_output[line_index + icount].split()) < 8):
            icount += 1
            bcounter += 6
            acounter = 0
            pass
        else:
            mo_line = qchem_output[line_index + icount].split()
            ccounter = bcounter
            for i in range(len(mo_line[4:])):
                c[acounter, ccounter] = float(mo_line[i+4])
                ccounter += 1
                count += 1
            acounter += 1
            icount += 1
    return c.T


def adf_get_overlap_matrix(qchem_output, nbf):
    s = np.matlib.zeros((nbf, nbf))
    for line in qchem_output:
        if "Overlap Matrix" in line:
            line_index = int(qchem_output.index(line)) + 2
            break
        else:
            pass
    acounter = 0
    bcounter = 0
    ccounter = 0
    count = 0
    icount = 0
    while count < nbf * nbf:
        if (len(qchem_output[line_index + icount].split()) < 6):
            icount += 1
            bcounter += 5
            acounter = 0
            pass
        else:
            mo_line = qchem_output[line_index + icount].split()
            ccounter = bcounter
            for i in range(len(mo_line[1:])):
                s[acounter, ccounter] = float(mo_line[i+1])
                ccounter += 1
                count += 1
            acounter += 1
            icount += 1
    return s


@jit
def adf_calculate_fock(c, s, nbf, e):
        # Calculate the Fock matrix
        sc_temp = sp.linalg.blas.dgemm(1.0, np.array(c), np.array(s), 1.0)
        sc = np.asmatrix(sc_temp)
        sce = np.matlib.zeros((nbf, nbf))
        for i in range(nbf):
            sce[i, :] = e[0, i] * sc[i, :]
        c_lu, ipiv, info = sp.linalg.lapack.dgetrf(c)
        ic_lu, info = sp.linalg.lapack.dgetri(c_lu, ipiv)
        f_temp = sp.linalg.blas.dgemm(1.0, np.array(ic_lu), np.array(sce), 1.0)
        f = np.asmatrix(f_temp)

        return f


def adf_get_homo_lumo(qchem_output):
        # Get HOMO and LUMO levels
        for line in qchem_output:
            if "There are" and "alpha" and "beta electrons" in line:
                homo = int(line.split()[2])
                lumo = homo + 1
                break
            else:
                pass
        return homo, lumo


@jit
def adf_build_cfrag(nbf_d, nbf_mon1, nbf_mon2, c_mon1, c_mon2):
        # Build the block matrix of monomer MO coefficients
        cfrag = np.matlib.zeros((nbf_d, nbf_d))
        for i in range(nbf_mon1):
            for j in range(nbf_mon1):
                cfrag[j, i] = c_mon1[j, i]

        for i in range(nbf_mon2):
            for j in range(nbf_mon2):
                cfrag[j + nbf_mon1, i + nbf_mon1] = c_mon2[j, i]

        return cfrag


@jit
def adf_calculate_coupling(system_name, f, s, homo_mon1, lumo_mon1, homo_mon2, lumo_mon2, cfrag, nbf_mon1):
        # Calculate intermolecular couplings
        fout = open(system_name + '.TI.out', 'w')
        TI_outfile_header(fout)
        fout.write("%-8s %-8s %-12s %-12s %-14s %-12s %-14s %-14s %-17s %-17s %-14s" \
                   % ("MO 1", "MO 2", "e 1, eV", "e 2, eV", "J_12, meV", "S_12", "e^eff_1, eV", "e^eff_2, eV",
                      "J^eff_12, meV",
                      "J^eff_12, cm-1", "dE_12, meV"))
        fout.write('\n')
        for i in range(homo_mon1 - 1, lumo_mon1):
            for j in range(homo_mon2 - 1 + nbf_mon1, lumo_mon2 + nbf_mon1):
                ctmp = cfrag[i, :] * f
                e1 = sp.linalg.blas.ddot(cfrag[i, :], ctmp) * 27.2114
                ctmp = cfrag[j, :] * f
                e2 = sp.linalg.blas.ddot(cfrag[j, :], ctmp) * 27.2114
                ctmp = cfrag[j, :] * f
                J12 = sp.linalg.blas.ddot(cfrag[i, :], ctmp) * 27.2114
                ctmp = cfrag[j, :] * s
                S12 = sp.linalg.blas.ddot(cfrag[i, :], ctmp)
                ee1 = 0.5 * ((e1 + e2) - 2.0 * J12 * S12 + (e1 - e2) * sqrt(1.0 - S12 * S12)) / (1.0 - S12 * S12)
                ee2 = 0.5 * ((e1 + e2) - 2.0 * J12 * S12 - (e1 - e2) * sqrt(1.0 - S12 * S12)) / (1.0 - S12 * S12)
                Je12 = (J12 - 0.5 * (e1 + e2) * S12) / (1.0 - S12 * S12)
                dE = sqrt((ee1 - ee2) ** 2 + (2 * Je12) ** 2)
                fout.write('%-8s %-8s %8.5f %12.5f %12.5f %14.5f %12.5f %14.5f %14.5f %17.5f %18.5f'
                           % (i + 1, j + 1 - nbf_mon1, e1, e2, J12 * 1000, S12, ee1, ee2, Je12 * 1000, Je12 * 8065.544,
                              dE * 1000))
                fout.write('\n')


def TI_g09(system_name, mon1_homo, mon1_lumo, mon2_homo, mon2_lumo):

    # Open fchk and log files
    dimer_fchk = system_name + '_A_B.fchk'
    mon1_fchk = system_name + '_A.fchk'
    mon2_fchk = system_name + '_B.fchk'

    try:
        dimer_log = system_name + '_A_B.log'
        log_check = open(dimer_log, 'r')
        log_check.close()
    except IOError:
        dimer_log = system_name + '_A_B.out'
        log_check = open(dimer_log, 'r')
        log_check.close()

    d_fchk = open(dimer_fchk, 'r')
    d_fchk_line = []
    for line in d_fchk:
        d_fchk_line.append(line)
    d_fchk.close()

    mon1_fchk = open(mon1_fchk, 'r')
    mon1_fchk_line = []
    for line in mon1_fchk:
        mon1_fchk_line.append(line)
    mon1_fchk.close()

    mon2_fchk = open(mon2_fchk, 'r')
    mon2_fchk_line = []
    for line in mon2_fchk:
        mon2_fchk_line.append(line)
    mon2_fchk.close()

    # Determine if restricted or unrestricted
    restrict_line = d_fchk_line[1].split()
    if restrict_line[1][0] == "U":
        TI_g09_unrestricted(system_name, dimer_log, d_fchk_line, mon1_fchk_line, mon2_fchk_line, mon1_homo,
                              mon1_lumo, mon2_homo, mon2_lumo)
    else:
        TI_g09_restricted(system_name, dimer_log, d_fchk_line, mon1_fchk_line, mon2_fchk_line, mon1_homo,
                            mon1_lumo, mon2_homo, mon2_lumo)


def TI_g09_restricted(system_name, dimer_log, d_fchk_line, mon1_fchk_line, mon2_fchk_line, homo_mon1, lumo_mon1, homo_mon2, lumo_mon2):
    # Get dimer values
    nbf_d = g09_get_basis_functions(d_fchk_line)
    e_d = g09_get_alpha_orbital_energies(d_fchk_line, nbf_d)
    c_d = g09_get_MO_coefficients(d_fchk_line, nbf_d)
    s_d = get_overlap(dimer_log, nbf_d)
    f_d = g09_calculate_fock(c_d, s_d, nbf_d, e_d)

    # Get monomer values
    nbf_mon1 = g09_get_basis_functions(mon1_fchk_line)
    if homo_mon1 == 0:
        homo_mon1, lumo_mon1 = g09_get_homo_lumo(mon1_fchk_line)
    else:
        pass
    c_mon1 = g09_get_MO_coefficients(mon1_fchk_line, nbf_mon1)
    nbf_mon2 = g09_get_basis_functions(mon2_fchk_line)
    if homo_mon2 == 0:
        homo_mon2, lumo_mon2 = g09_get_homo_lumo(mon2_fchk_line)
    else:
        pass
    c_mon2 = g09_get_MO_coefficients(mon2_fchk_line, nbf_mon2)
    cfrag = g09_build_cfrag(nbf_d, nbf_mon1, nbf_mon2, c_mon1, c_mon2)

    g09_calculate_coupling(system_name, f_d, s_d, homo_mon1, lumo_mon1, homo_mon2, lumo_mon2, cfrag, nbf_mon1)


def TI_g09_unrestricted(system_name, dimer_log, d_fchk_line, mon1_fchk_line, mon2_fchk_line, homo_mon1, lumo_mon1, homo_mon2, lumo_mon2):
    # Get dimer values
    nbf_d = g09_get_basis_functions(d_fchk_line)
    e_a_d = g09_get_alpha_orbital_energies(d_fchk_line, nbf_d)
    e_b_d = g09_get_beta_orbital_energies(d_fchk_line, nbf_d)
    c_a_d = g09_get_MO_coefficients(d_fchk_line, nbf_d)
    c_b_d = g09_get_beta_MO_coefficients(d_fchk_line, nbf_d)
    s_d = get_overlap(dimer_log, nbf_d)
    f_a_d = g09_calculate_fock(c_a_d, s_d, nbf_d, e_a_d)
    f_b_d = g09_calculate_fock(c_b_d, s_d, nbf_d, e_b_d)

    # Get monomer values
    nbf_mon1 = g09_get_basis_functions(mon1_fchk_line)
    if homo_mon1 == 0:
        homo_a_mon1, lumo_a_mon1 = g09_get_homo_lumo(mon1_fchk_line)
        homo_b_mon1, lumo_b_mon1 = g09_get_beta_homo_lumo(mon1_fchk_line)
    else:
        homo_a_mon1 = homo_mon1
        homo_b_mon1 = homo_mon1
        lumo_a_mon1 = lumo_mon1
        lumo_b_mon1 = lumo_mon1
    c_a_mon1 = g09_get_MO_coefficients(mon1_fchk_line, nbf_mon1)
    c_b_mon1 = g09_get_beta_MO_coefficients(mon1_fchk_line, nbf_mon1)
    nbf_mon2 = g09_get_basis_functions(mon2_fchk_line)
    if homo_mon1 == 0:
        homo_a_mon2, lumo_a_mon2 = g09_get_homo_lumo(mon2_fchk_line)
        homo_b_mon2, lumo_b_mon2 = g09_get_beta_homo_lumo(mon2_fchk_line)
    else:
        homo_a_mon2 = homo_mon2
        homo_b_mon2 = homo_mon2
        lumo_a_mon2 = lumo_mon2
        lumo_b_mon2 = lumo_mon2
    c_a_mon2 = g09_get_MO_coefficients(mon2_fchk_line, nbf_mon2)
    c_b_mon2 = g09_get_beta_MO_coefficients(mon2_fchk_line, nbf_mon2)
    cfrag_a = g09_build_cfrag(nbf_d, nbf_mon1, nbf_mon2, c_a_mon1, c_a_mon2)
    cfrag_b = g09_build_cfrag(nbf_d, nbf_mon1, nbf_mon2, c_b_mon1, c_b_mon2)

    g09_calculate_unrestricted_coupling(system_name, f_a_d, f_b_d, s_d, homo_a_mon1, homo_b_mon1, lumo_a_mon1, lumo_b_mon1,
                           homo_a_mon2, homo_b_mon2, lumo_a_mon2, lumo_b_mon2, cfrag_a, cfrag_b, nbf_mon1)


def g09_get_basis_functions(g09_output):
    # Get the number of basis functions
    for line in g09_output:
        if "Number of basis functions" in line:
            nbf = int(line.split()[5])
            break
        else: pass

    return nbf


def g09_get_homo_lumo(g09_output):
    # Get HOMO and LUMO levels
    for line in g09_output:
        if "Number of alpha electrons" in line:
            homo = int(line.split()[5])
            lumo = homo + 1
            break
        else: pass

    return homo, lumo


def g09_get_beta_homo_lumo(g09_output):
    # Get HOMO and LUMO levels
    for line in g09_output:
        if "Number of beta electrons" in line:
            homo = int(line.split()[5])
            lumo = homo + 1
            break
        else: pass

    return homo, lumo


def g09_get_alpha_orbital_energies(g09_output, nbf):
    # Get alpha orbital energies
    e = np.matlib.zeros(nbf)
    for line in g09_output:
        if "Alpha Orbital Energies" in line:
            line_index = int(g09_output.index(line))+1
            break
        else: pass
    count = 0
    icount = 0
    while count < nbf:
        mo_line = g09_output[line_index+icount].split()
        icount += 1
        for i in range(len(mo_line)):
            e[0, count] = float(mo_line[i])
            count += 1

    return e


def g09_get_beta_orbital_energies(g09_output, nbf):
    # Get alpha orbital energies
    e = np.matlib.zeros(nbf)
    for line in g09_output:
        if "Beta Orbital Energies" in line:
            line_index = int(g09_output.index(line))+1
            break
        else: pass
    count = 0
    icount = 0
    while count < nbf:
        mo_line = g09_output[line_index+icount].split()
        icount += 1
        for i in range(len(mo_line)):
            e[0, count] = float(mo_line[i])
            count += 1

    return e


def g09_get_MO_coefficients(g09_output, nbf):
    # Get alpha MO coefficients
    c_temp = np.matlib.zeros(nbf*nbf)
    for line in g09_output:
        if "Alpha MO coefficients" in line:
            line_index = int(g09_output.index(line))+1
            break
        else: pass
    count = 0
    icount = 0
    while count < nbf*nbf:
        mo_line = g09_output[line_index+icount].split()
        icount += 1
        for i in range(len(mo_line)):
            c_temp[0, count] = float(mo_line[i])
            count += 1
    c = c_temp.reshape((nbf, nbf))

    return c


def g09_get_beta_MO_coefficients(g09_output, nbf):
    # Get alpha MO coefficients
    c_temp = np.matlib.zeros(nbf*nbf)
    for line in g09_output:
        if "Beta MO coefficients" in line:
            line_index = int(g09_output.index(line))+1
            break
        else: pass
    count = 0
    icount = 0
    while count < nbf*nbf:
        mo_line = g09_output[line_index+icount].split()
        icount += 1
        for i in range(len(mo_line)):
            c_temp[0, count] = float(mo_line[i])
            count += 1
    c = c_temp.reshape((nbf, nbf))

    return c


@jit
def g09_calculate_fock(c, s, nbf, e):
    # Calculate the Fock matrix
    sc_temp = sp.linalg.blas.dgemm(1.0, np.array(c), np.array(s), 1.0)
    sc = np.asmatrix(sc_temp)
    sce = np.matlib.zeros((nbf, nbf))
    for i in range(nbf):
        sce[i,:] = e[0, i]*sc[i, :]
    c_lu, ipiv, info = sp.linalg.lapack.dgetrf(c)
    ic_lu, info = sp.linalg.lapack.dgetri(c_lu, ipiv)
    f_temp = sp.linalg.blas.dgemm(1.0, np.array(ic_lu), np.array(sce), 1.0)
    f = np.asmatrix(f_temp)

    return f


@jit
def g09_build_cfrag(nbf_d, nbf_mon1, nbf_mon2, c_mon1, c_mon2):
    # Build the block matrix of monomer MO coefficients
    cfrag = np.matlib.zeros((nbf_d, nbf_d))
    for i in range(nbf_mon1):
        for j in range(nbf_mon1):
            cfrag[j, i] = c_mon1[j, i]

    for i in range(nbf_mon2):
        for j in range(nbf_mon2):
            cfrag[j+nbf_mon1, i+nbf_mon1] = c_mon2[j, i]

    return cfrag


@jit
def g09_calculate_coupling(system_name, f, s, homo_mon1, lumo_mon1, homo_mon2, lumo_mon2, cfrag, nbf_mon1):
    # Calculate intermolecular couplings
    fout = open(system_name+'.TI.out','w')
    TI_outfile_header(fout)
    fout.write("%-8s %-8s %-12s %-12s %-14s %-12s %-14s %-14s %-17s %-17s %-14s" \
          % ("MO 1", "MO 2", "e 1, eV", "e 2, eV", "J_12, meV", "S_12", "e^eff_1, eV", "e^eff_2, eV", "J^eff_12, meV",
             "J^eff_12, cm-1", "dE_12, meV"))
    fout.write('\n')
    for i in range(homo_mon1-1, lumo_mon1):
        for j in range(homo_mon2-1+nbf_mon1, lumo_mon2+nbf_mon1):
            ctmp = cfrag[i, :] * f
            e1 = sp.linalg.blas.ddot(cfrag[i, :], ctmp)*27.2116
            ctmp = cfrag[j, :] * f
            e2 = sp.linalg.blas.ddot(cfrag[j, :], ctmp)*27.2116
            ctmp = cfrag[j, :] * f
            J12 = sp.linalg.blas.ddot(cfrag[i, :], ctmp)*27.2116
            ctmp = cfrag[j, :] * s
            S12 = sp.linalg.blas.ddot(cfrag[i, :], ctmp)
            ee1 = 0.5*((e1+e2)-2.0*J12*S12+(e1-e2)*sqrt(1.0-S12*S12))/(1.0-S12*S12)
            ee2 = 0.5*((e1+e2)-2.0*J12*S12-(e1-e2)*sqrt(1.0-S12*S12))/(1.0-S12*S12)
            Je12 = (J12-0.5*(e1+e2)*S12)/(1.0-S12*S12)
            dE = sqrt((ee1-ee2)**2 + (2*Je12)**2)
            fout.write('%-8s %-8s %8.5f %12.5f %12.5f %14.5f %12.5f %14.5f %14.5f %17.5f %18.5f'
                % (i+1, j+1-nbf_mon1, e1, e2, J12*1000, S12, ee1, ee2, Je12*1000, Je12*8065.544, dE*1000))
            fout.write('\n')


@jit
def g09_calculate_unrestricted_coupling(system_name, f_a, f_b, s, homo_a_mon1, homo_b_mon1, lumo_a_mon1, lumo_b_mon1, homo_a_mon2, homo_b_mon2, lumo_a_mon2, lumo_b_mon2, cfrag_a, cfrag_b, nbf_mon1):
    # Calculate intermolecular couplings
    fout_a = open(system_name+'_alpha.TI.out','w')
    fout_b = open(system_name+'_beta.TI.out','w')
    TI_outfile_header(fout_a)
    TI_outfile_header(fout_b)
    fout_a.write("%-8s %-8s %-12s %-12s %-14s %-12s %-14s %-14s %-17s %-17s %-14s" \
          % ("MO 1", "MO 2", "e 1, eV", "e 2, eV", "J_12, meV", "S_12", "e^eff_1, eV", "e^eff_2, eV", "J^eff_12, meV",
             "J^eff_12, cm-1", "dE_12, meV"))
    fout_b.write("%-8s %-8s %-12s %-12s %-14s %-12s %-14s %-14s %-17s %-17s %-14s" \
          % ("MO 1", "MO 2", "e 1, eV", "e 2, eV", "J_12, meV", "S_12", "e^eff_1, eV", "e^eff_2, eV", "J^eff_12, meV",
             "J^eff_12, cm-1", "dE_12, meV"))

    fout_a.write('\n')
    fout_b.write('\n')
    for i in range(homo_a_mon1-1, lumo_a_mon1):
        for j in range(homo_a_mon2-1+nbf_mon1, lumo_a_mon2+nbf_mon1):
            ctmp = cfrag_a[i, :] * f_a
            e1 = sp.linalg.blas.ddot(cfrag_a[i, :], ctmp)*27.2116
            ctmp = cfrag_a[j, :] * f_a
            e2 = sp.linalg.blas.ddot(cfrag_a[j, :], ctmp)*27.2116
            ctmp = cfrag_a[j, :] * f_a
            J12 = sp.linalg.blas.ddot(cfrag_a[i, :], ctmp)*27.2116
            ctmp = cfrag_a[j, :] * s
            S12 = sp.linalg.blas.ddot(cfrag_a[i, :], ctmp)
            ee1 = 0.5*((e1+e2)-2.0*J12*S12+(e1-e2)*sqrt(1.0-S12*S12))/(1.0-S12*S12)
            ee2 = 0.5*((e1+e2)-2.0*J12*S12-(e1-e2)*sqrt(1.0-S12*S12))/(1.0-S12*S12)
            Je12 = (J12-0.5*(e1+e2)*S12)/(1.0-S12*S12)
            dE = sqrt((ee1-ee2)**2 + (2*Je12)**2)
            fout_a.write('%-8s %-8s %8.5f %12.5f %12.5f %14.5f %12.5f %14.5f %14.5f %17.5f %18.5f'
                % (i+1, j+1-nbf_mon1, e1, e2, J12*1000, S12, ee1, ee2, Je12*1000, Je12*8065.544, dE*1000))
            fout_a.write('\n')

    for i in range(homo_b_mon1-1, lumo_b_mon1):
        for j in range(homo_b_mon2-1+nbf_mon1, lumo_b_mon2+nbf_mon1):
            ctmp = cfrag_b[i, :] * f_b
            e1 = sp.linalg.blas.ddot(cfrag_b[i, :], ctmp)*27.2116
            ctmp = cfrag_b[j, :] * f_b
            e2 = sp.linalg.blas.ddot(cfrag_b[j, :], ctmp)*27.2116
            ctmp = cfrag_b[j, :] * f_b
            J12 = sp.linalg.blas.ddot(cfrag_b[i, :], ctmp)*27.2116
            ctmp = cfrag_b[j, :] * s
            S12 = sp.linalg.blas.ddot(cfrag_b[i, :], ctmp)
            ee1 = 0.5*((e1+e2)-2.0*J12*S12+(e1-e2)*sqrt(1.0-S12*S12))/(1.0-S12*S12)
            ee2 = 0.5*((e1+e2)-2.0*J12*S12-(e1-e2)*sqrt(1.0-S12*S12))/(1.0-S12*S12)
            Je12 = (J12-0.5*(e1+e2)*S12)/(1.0-S12*S12)
            dE = sqrt((ee1-ee2)**2 + (2*Je12)**2)
            fout_b.write('%-8s %-8s %8.5f %12.5f %12.5f %14.5f %12.5f %14.5f %14.5f %17.5f %18.5f'
                % (i+1, j+1-nbf_mon1, e1, e2, J12*1000, S12, ee1, ee2, Je12*1000, Je12*8065.544, dE*1000))
            fout_b.write('\n')


def TI_orca(system_name, homo_mon1, lumo_mon1, homo_mon2, lumo_mon2):
    # Open dimer out file
    d_file = open(system_name+'.out', 'r')
    d_file_line = []
    for line in d_file:
        d_file_line.append(line)

    d_file.close()

    nbf_mon1, nbf_mon2, nbf_d = orca_get_basis_functions(d_file_line)
    if homo_mon1 == 0:
        homo_mon1, lumo_mon1, homo_mon2, lumo_mon2 = orca_get_homo_lumo(d_file_line)
    else:
        pass
    e_d = orca_get_alpha_orbital_energies(d_file_line, nbf_d)
    c_temp_d, c_temp_mon1, c_temp_mon2, cfrag = orca_get_MO_coefficients(d_file_line, nbf_d, nbf_mon1, nbf_mon2)
    s = orca_get_overlap_matrix(d_file_line, nbf_d)
    f = orca_calculate_fock(c_temp_d, e_d, s, nbf_d)
    orca_calculate_coupling(system_name, homo_mon1, lumo_mon1, homo_mon2, lumo_mon2, cfrag, f, s, nbf_mon1)


def orca_get_basis_functions(orca_output_lines):
    # Get number of basis functions
    counter = 0
    for line in orca_output_lines:
        if counter == 0:
            if "# of contracted basis functions" in line:
                nbf_mon1 = int(line.split()[6])
                counter = 1
        elif counter == 1:
            if "# of contracted basis functions" in line:
                nbf_mon2 = int(line.split()[6])
                counter = 2
        elif counter == 2:
            if "# of contracted basis functions" in line:
                nbf_d = int(line.split()[6])
                counter = 3
                break
        else: pass

    return nbf_mon1, nbf_mon2, nbf_d


def orca_get_homo_lumo(orca_output_lines):
    # Get HOMOs and LUMOs
    counter = 0
    for line in orca_output_lines:
        if counter == 0:
            if "N(Alpha)" in line:
                homo_mon1 = int(round(float(line.split()[2])))
                lumo_mon1 = homo_mon1+1
                counter = 1
        elif counter == 1:
            if "N(Alpha)" in line:
                homo_mon2 = int(round(float(line.split()[2])))
                lumo_mon2 = homo_mon2+1
                counter = 2
                break
        else: pass

    return homo_mon1, lumo_mon1, homo_mon2, lumo_mon2


def orca_get_alpha_orbital_energies(orca_output_lines, nbf_d):
    # Get alpha orbital energies
    e_d = np.matlib.zeros(nbf_d)

    counter = 0
    for i in range(len(orca_output_lines)):
        line = orca_output_lines[i]
        if counter == 0:
            if "ORBITAL ENERGIES" in line:
                mon1_index = i+4
                counter = 1
        elif counter == 1:
            if "ORBITAL ENERGIES" in line:
                mon2_index = i+4
                counter = 2
        elif counter == 2:
            if "ORBITAL ENERGIES" in line:
                d_index = i+4
                counter = 3
                break
        else: pass

    count = 0
    while count < nbf_d:
        e_d[0, count] = float(orca_output_lines[d_index+count].split()[2])
        count += 1

    return e_d


def orca_get_MO_coefficients(orca_output_lines, nbf_d, nbf_mon1, nbf_mon2):
    # Get orbital coefficients
    c_temp_mon1 = np.matlib.zeros((nbf_mon1, nbf_mon1))
    c_temp_mon2 = np.matlib.zeros((nbf_mon2, nbf_mon2))
    c_temp_d = np.matlib.zeros((nbf_d, nbf_d))

    counter = 0
    for i in range(len(orca_output_lines)):
        line = orca_output_lines[i]
        if counter == 0:
            if "MOLECULAR ORBITALS" in line:
                mon1_index = i+6
                counter = 1
        elif counter == 1:
            if "MOLECULAR ORBITALS" in line:
                mon2_index = i+6
                counter = 2
        elif counter == 2:
            if "MOLECULAR ORBITALS" in line:
                d_index = i+6
                counter = 3
                break
        else: pass

    acounter = 0
    bcounter = 0
    ccounter = 0
    looping = 1
    while looping:
        if acounter < nbf_mon1:
            temp_line = orca_output_lines[mon1_index+acounter][16:]
            n = 10
            temp_line = [temp_line[i:i+n] for i in range(0, len(temp_line), n)]
            for i in range(len(temp_line)):
                temp_line[i] = temp_line[i].strip()
            temp_line = [_f for _f in temp_line if _f]
            for i in range(len(temp_line)):
                c_temp_mon1[acounter, bcounter] = float(temp_line[i])
                bcounter += 1
            bcounter = ccounter
            acounter += 1
            print(temp_line, file=sys.stdout)
            for i in range(len(temp_line)):
                c_temp_mon1[acounter, bcounter] = float(temp_line[i])
                bcounter += 1
            bcounter = ccounter
            acounter += 1
        elif acounter == nbf_mon1:
            if ccounter+6 >= nbf_mon1:
                looping = 0
            else:
                mon1_index += nbf_mon1+4
                ccounter += 6
                bcounter = ccounter
                acounter = 0
    acounter = 0
    bcounter = 0
    ccounter = 0
    looping = 1
    while looping:
        if acounter < nbf_mon2:
            temp_line = orca_output_lines[mon2_index+acounter][16:]
            n = 10
            temp_line = [temp_line[i:i+n] for i in range(0, len(temp_line), n)]
            for i in range(len(temp_line)):
                temp_line[i] = temp_line[i].strip()
            temp_line = [_f for _f in temp_line if _f]
            for i in range(len(temp_line)):
                c_temp_mon2[acounter, bcounter] = float(temp_line[i])
                bcounter += 1
            bcounter = ccounter
            acounter += 1
        elif acounter == nbf_mon2:
            if ccounter+6 >= nbf_mon2:
                looping = 0
            else:
                mon2_index += nbf_mon2+4
                ccounter += 6
                bcounter = ccounter
                acounter = 0
    acounter = 0
    bcounter = 0
    ccounter = 0
    looping = 1
    while looping:
        if acounter < nbf_d:
            temp_line = orca_output_lines[d_index+acounter][16:]
            n = 10
            temp_line = [temp_line[i:i+n] for i in range(0, len(temp_line), n)]
            for i in range(len(temp_line)):
                temp_line[i] = temp_line[i].strip()
            temp_line = [_f for _f in temp_line if _f]
            for i in range(len(temp_line)):
                c_temp_d[acounter, bcounter] = float(temp_line[i])
                bcounter += 1
            bcounter = ccounter
            acounter += 1
        elif acounter == nbf_d:
            if ccounter+6 >= nbf_d:
                looping = 0
            else:
                d_index += nbf_d+4
                ccounter += 6
                bcounter = ccounter
                acounter = 0

    cfrag = np.matlib.zeros((nbf_d, nbf_d))
    for i in range(nbf_mon1):
        for j in range(nbf_mon1):
            cfrag[j, i] = c_temp_mon1[j, i]
    for i in range(nbf_mon2):
        for j in range(nbf_mon2):
            cfrag[j+nbf_mon1, i+nbf_mon1] = c_temp_mon2[j, i]

    return c_temp_d.T, c_temp_mon1, c_temp_mon2, cfrag.T


def orca_get_overlap_matrix(orca_output_lines, nbf_d):
    # Get overlap matrix
    s = np.matlib.zeros((nbf_d, nbf_d))
    for i in range(len(orca_output_lines)):
        line = orca_output_lines[i]
        if "OVERLAP MATRIX" in line:
            d_index = i+3
            break
        else: pass

    acounter = 0
    bcounter = 0
    ccounter = 0
    counter = 0
    looping = 1
    nbf_flag = 1
    #print >> sys.stdout, "nbf_d=", nbf_d
    while looping:
            if acounter < nbf_d:
                temp_line = orca_output_lines[d_index+acounter].split()
                temp_line.pop(0)
                for i in range(len(temp_line)):
                    temp_line[i] = temp_line[i].strip()
                for i in range(len(temp_line)):
                    s[acounter,bcounter] = float(temp_line[i])
                    bcounter += 1
                bcounter = ccounter
                acounter += 1
            elif acounter == nbf_d:
                if ccounter+6 >= nbf_d:
                    looping = 0
                else:
                    d_index += nbf_d+1
                    ccounter += 6
                    bcounter = ccounter
                    acounter = 0

    return s


@jit
def orca_calculate_fock(c_temp_d, e_d, s, nbf_d):
    # Calculate the dimer Fock matrix
    c = c_temp_d
    e = e_d
    sc_temp = sp.linalg.blas.dgemm(1.0, np.array(c), np.array(s), 1.0)
    sc = np.asmatrix(sc_temp)
    sce = np.matlib.zeros((nbf_d, nbf_d))
    for i in range(nbf_d):
        sce[i,:] = e[0, i]*sc[i, :]
    c_lu, ipiv, info = sp.linalg.lapack.dgetrf(c)
    ic_lu, info = sp.linalg.lapack.dgetri(c_lu, ipiv)
    f_temp = sp.linalg.blas.dgemm(1.0, np.array(ic_lu), np.array(sce), 1.0)
    f = np.asmatrix(f_temp)

    return f


@jit
def orca_calculate_coupling(system_name, homo_mon1, lumo_mon1, homo_mon2, lumo_mon2, cfrag, f, s, nbf_mon1):
    # Calculate intermolecular couplings
    fout = open(system_name+'.TI.out','w')
    TI_outfile_header(fout)
    fout.write("%-8s %-8s %-12s %-12s %-14s %-12s %-14s %-14s %-17s %-17s %-14s" \
          % ("MO 1", "MO 2", "e 1, eV", "e 2, eV", "J_12, meV", "S_12", "e^eff_1, eV", "e^eff_2, eV", "J^eff_12, meV",
             "J^eff_12, cm-1", "dE_12, meV"))
    fout.write('\n')
    for i in range(homo_mon1-1, lumo_mon1):
        for j in range(homo_mon2-1+nbf_mon1, lumo_mon2+nbf_mon1):
            ctmp = cfrag[i, :] * f
            e1 = sp.linalg.blas.ddot(cfrag[i, :], ctmp)*27.2116
            ctmp = cfrag[j, :] * f
            e2 = sp.linalg.blas.ddot(cfrag[j, :], ctmp)*27.2116
            ctmp = cfrag[j, :] * f
            J12 = sp.linalg.blas.ddot(cfrag[i, :], ctmp)*27.2116
            ctmp = cfrag[j, :] * s
            S12 = sp.linalg.blas.ddot(cfrag[i, :], ctmp)
            ee1 = 0.5*((e1+e2)-2.0*J12*S12+(e1-e2)*sqrt(1.0-S12*S12))/(1.0-S12*S12)
            ee2 = 0.5*((e1+e2)-2.0*J12*S12-(e1-e2)*sqrt(1.0-S12*S12))/(1.0-S12*S12)
            Je12 = (J12-0.5*(e1+e2)*S12)/(1.0-S12*S12)
            dE = sqrt((ee1-ee2)**2 + (2*Je12)**2)
            fout.write('%-8s %-8s %8.5f %12.5f %12.5f %14.5f %12.5f %14.5f %14.5f %17.5f %18.5f'
                % (i+1, j+1-nbf_mon1, e1, e2, J12*1000, S12, ee1, ee2, Je12*1000, Je12*8065.544, dE*1000))
            fout.write('\n')


def TI_qchem(system_name, homo_mon1, lumo_mon1, homo_mon2, lumo_mon2):
    # Open dimer out file
    d_file = open(system_name+'_A_B.out', 'r')
    d_file_line = []
    for line in d_file:
        d_file_line.append(line)

    d_file.close()

    # Get dimer values
    nbf_d = qchem_get_basis_functions(d_file_line)
    e_d = qchem_get_alpha_orbital_energies(d_file_line, nbf_d)
    c_d = qchem_get_MO_coefficients(d_file_line, nbf_d)
    s_d = qchem_get_overlap_matrix(d_file_line, nbf_d)
    f_d = qchem_calculate_fock(c_d, s_d, nbf_d, e_d)

    d_file_line = None

    # Open monomer A out file
    mon1_file = open(system_name+'_A.out', 'r')
    mon1_file_line = []
    for line in mon1_file:
        mon1_file_line.append(line)

    mon1_file.close()

    # Get monomer A values
    nbf_mon1 = qchem_get_basis_functions(mon1_file_line)
    if homo_mon1 == 0:
        homo_mon1, lumo_mon1 = qchem_get_homo_lumo(mon1_file_line)
    else:
        pass
    c_mon1 = qchem_get_MO_coefficients(mon1_file_line, nbf_mon1)

    mon1_file_line = None

    # Open monomer B out file
    mon2_file = open(system_name+'_B.out', 'r')
    mon2_file_line = []
    for line in mon2_file:
        mon2_file_line.append(line)

    mon2_file.close()

    # Get monomer B values
    nbf_mon2 = qchem_get_basis_functions(mon2_file_line)
    if homo_mon2 == 0:
        homo_mon2, lumo_mon2 = qchem_get_homo_lumo(mon2_file_line)
    else:
        pass
    c_mon2 = qchem_get_MO_coefficients(mon2_file_line, nbf_mon2)

    mon2_file_line = None

    cfrag = qchem_build_cfrag(nbf_d, nbf_mon1, nbf_mon2, c_mon1, c_mon2)

    qchem_calculate_coupling(system_name, f_d, s_d, homo_mon1, lumo_mon1, homo_mon2, lumo_mon2, cfrag, nbf_mon1)


def qchem_get_basis_functions(qchem_output):
        # Get the number of basis functions
        for line in qchem_output:
            if "There are" and "shells and" and "basis functions" in line:
            #if "Number of Cartesian basis functions" in line:          //Not guaranteed to use Cartesian basis functions
                nbf = int(line.split()[5])
                break
            else:
                pass

        return nbf


def qchem_get_alpha_orbital_energies(qchem_output, nbf):
        # Get alpha orbital energies
        e = np.matlib.zeros(nbf)
        for line in qchem_output:
            if "Orbital Energies (a.u.)" in line:
                line_index = int(qchem_output.index(line)) + 5
                break
            else:
                pass
        count = 0
        icount = 0
        while count < nbf:
            if "Virtual" in qchem_output[line_index + icount]:
                icount += 1
                pass
            else:
                mo_line = qchem_output[line_index + icount].split()
                icount += 1
                for i in range(len(mo_line)):
                    e[0, count] = float(mo_line[i])
                    count += 1

        return e


def qchem_get_MO_coefficients_old(qchem_output, nbf):
        # Get alpha MO coefficients
        c_temp = np.matlib.zeros(nbf * nbf)
        for line in qchem_output:
            if "RESTRICTED (RHF) MOLECULAR ORBITAL COEFFICIENTS" in line:
                line_index = int(qchem_output.index(line)) + 2
                break
            else:
                pass
        count = 0
        icount = 0
        while count < nbf * nbf:
            if ("eigenvalues:" in qchem_output[line_index + icount]) or (len(qchem_output[line_index + icount].split()) < 8):
                icount += 1
                pass
            else:
                mo_line = qchem_output[line_index + icount].split()
                icount += 1
                for i in range(len(mo_line[4:])):
                    c_temp[0, count] = float(mo_line[i+4])
                    count += 1
        c = c_temp.reshape((nbf, nbf))

        return c


def qchem_get_MO_coefficients(qchem_output, nbf):
    c = np.matlib.zeros((nbf, nbf))
    for line in qchem_output:
        if "RESTRICTED (RHF) MOLECULAR ORBITAL COEFFICIENTS" in line:
            line_index = int(qchem_output.index(line)) + 2
            break
        else:
            pass
    acounter = 0
    bcounter = 0
    ccounter = 0
    count = 0
    icount = 0
    element_list = ['H',
                    'B', 'C', 'N', 'O', 'F',
                    'Al', 'Si', 'P', 'S', 'Cl',
                    'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br',
                    'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',
                    'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi']
    while count < nbf * nbf:
        # print(qchem_output[line_index + icount])
        if ("eigenvalues" in qchem_output[line_index + icount]):
            icount += 1
            acounter = 0
            pass
        elif qchem_output[line_index + icount].strip() == "":
            break
        elif len(qchem_output[line_index + icount].split()) == 1:
            icount += 1
            bcounter += 6
            acounter = 0
            pass
        elif (qchem_output[line_index + icount].split()[1] in element_list):
            mo_line = qchem_output[line_index + icount]
            ccounter = bcounter
            dcounter = 0
            while dcounter < 6:
                try:
                    c[acounter, ccounter] = float(mo_line[18:28].strip())
                    ccounter += 1
                    count += 1
                    dcounter += 1
                except (ValueError, IndexError):
                    dcounter += 1
                    pass
                try:
                    c[acounter, ccounter] = float(mo_line[28:38].strip())
                    ccounter += 1
                    count += 1
                    dcounter += 1
                except ValueError:
                    dcounter += 1
                    pass
                try:
                    c[acounter, ccounter] = float(mo_line[38:48].strip())
                    ccounter += 1
                    count += 1
                    dcounter += 1
                except ValueError:
                    dcounter += 1
                    pass
                try:
                    c[acounter, ccounter] = float(mo_line[48:58].strip())
                    ccounter += 1
                    count += 1
                    dcounter += 1
                except ValueError:
                    dcounter += 1
                    pass
                try:
                    c[acounter, ccounter] = float(mo_line[58:68].strip())
                    ccounter += 1
                    count += 1
                    dcounter += 1
                except ValueError:
                    dcounter += 1
                    pass
                try:
                    c[acounter, ccounter] = float(mo_line[68:78].strip())
                    ccounter += 1
                    count += 1
                    dcounter += 1
                except ValueError:
                    dcounter += 1
                    pass

            acounter += 1
            icount += 1
        else:
            icount += 1
            bcounter += 6
            acounter = 0
            pass
    return c.T


def qchem_get_overlap_matrix(qchem_output, nbf):
    s = np.matlib.zeros((nbf, nbf))
    for line in qchem_output:
        if "Overlap Matrix" in line:
            line_index = int(qchem_output.index(line)) + 2
            break
        else:
            pass
    acounter = 0
    bcounter = 0
    ccounter = 0
    count = 0
    icount = 0
    while count < nbf * nbf:
        #if (len(qchem_output[line_index + icount].split()) < 6):       // Used in a previous test where all lines were filled
        if ('.' not in qchem_output[line_index + icount].split()[0]) and \
                ('.' not in qchem_output[line_index + icount].split()[1]):
            icount += 1
            bcounter += 5
            acounter = 0
            pass
        else:
            temp_line = qchem_output[line_index + icount]
            mo_line = qchem_output[line_index + icount].split()
            ccounter = bcounter
            for i in range(len(mo_line[1:])):
                s[acounter, ccounter] = float(mo_line[i+1])
                ccounter += 1
                count += 1
            acounter += 1
            icount += 1
    return s


@jit
def qchem_calculate_fock(c, s, nbf, e):
        # Calculate the Fock matrix
        sc_temp = sp.linalg.blas.dgemm(1.0, np.array(c), np.array(s), 1.0)
        sc = np.asmatrix(sc_temp)
        sce = np.matlib.zeros((nbf, nbf))
        for i in range(nbf):
            sce[i, :] = e[0, i] * sc[i, :]
        c_lu, ipiv, info = sp.linalg.lapack.dgetrf(c)
        ic_lu, info = sp.linalg.lapack.dgetri(c_lu, ipiv)
        f_temp = sp.linalg.blas.dgemm(1.0, np.array(ic_lu), np.array(sce), 1.0)
        f = np.asmatrix(f_temp)

        return f


def qchem_get_homo_lumo(qchem_output):
        # Get HOMO and LUMO levels
        for line in qchem_output:
            if "There are" and "alpha" and "beta electrons" in line:
                homo = int(line.split()[2])
                lumo = homo + 1
                break
            else:
                pass
        return homo, lumo


@jit
def qchem_build_cfrag(nbf_d, nbf_mon1, nbf_mon2, c_mon1, c_mon2):
        # Build the block matrix of monomer MO coefficients
        cfrag = np.matlib.zeros((nbf_d, nbf_d))
        for i in range(nbf_mon1):
            for j in range(nbf_mon1):
                cfrag[j, i] = c_mon1[j, i]

        for i in range(nbf_mon2):
            for j in range(nbf_mon2):
                cfrag[j + nbf_mon1, i + nbf_mon1] = c_mon2[j, i]

        return cfrag


@jit
def qchem_calculate_coupling(system_name, f, s, homo_mon1, lumo_mon1, homo_mon2, lumo_mon2, cfrag, nbf_mon1):
        # Calculate intermolecular couplings
        fout = open(system_name + '.TI.out', 'w')
        TI_outfile_header(fout)
        fout.write("%-8s %-8s %-12s %-12s %-14s %-12s %-14s %-14s %-17s %-17s %-14s" \
                   % ("MO 1", "MO 2", "e 1, eV", "e 2, eV", "J_12, meV", "S_12", "e^eff_1, eV", "e^eff_2, eV",
                      "J^eff_12, meV",
                      "J^eff_12, cm-1", "dE_12, meV"))
        fout.write('\n')
        for i in range(homo_mon1 - 1, lumo_mon1):
            for j in range(homo_mon2 - 1 + nbf_mon1, lumo_mon2 + nbf_mon1):
                ctmp = cfrag[i, :] * f
                e1 = sp.linalg.blas.ddot(cfrag[i, :], ctmp) * 27.2114
                ctmp = cfrag[j, :] * f
                e2 = sp.linalg.blas.ddot(cfrag[j, :], ctmp) * 27.2114
                ctmp = cfrag[j, :] * f
                J12 = sp.linalg.blas.ddot(cfrag[i, :], ctmp) * 27.2114
                ctmp = cfrag[j, :] * s
                S12 = sp.linalg.blas.ddot(cfrag[i, :], ctmp)
                ee1 = 0.5 * ((e1 + e2) - 2.0 * J12 * S12 + (e1 - e2) * sqrt(1.0 - S12 * S12)) / (1.0 - S12 * S12)
                ee2 = 0.5 * ((e1 + e2) - 2.0 * J12 * S12 - (e1 - e2) * sqrt(1.0 - S12 * S12)) / (1.0 - S12 * S12)
                Je12 = (J12 - 0.5 * (e1 + e2) * S12) / (1.0 - S12 * S12)
                dE = sqrt((ee1 - ee2) ** 2 + (2 * Je12) ** 2)
                fout.write('%-8s %-8s %8.5f %12.5f %12.5f %14.5f %12.5f %14.5f %14.5f %17.5f %18.5f'
                           % (i + 1, j + 1 - nbf_mon1, e1, e2, J12 * 1000, S12, ee1, ee2, Je12 * 1000, Je12 * 8065.544,
                              dE * 1000))
                fout.write('\n')


def get_overlap(g09file, nbf):
    """
    Extracts the overlap matrix from a Gaussian logfile.
    Returns a numpy matrix.
    """
    logfile = open(g09file, 'r')
    data = logfile.read()
    overlap_matrix = np.zeros((nbf, nbf))
    # grab all text between  "Overlap ***" and "*** Kinetic"
    raw_overlap_string = re.findall(r'Overlap \*\*\*(.*?)\*\*\* Kinetic', data, re.DOTALL)
    raw_overlap_string = raw_overlap_string[0].replace('D', 'E')
    raw_overlap_elements = raw_overlap_string.split()
    matrix_elements = []
    for overlap_value in raw_overlap_elements:
        if 'E' in overlap_value:
            matrix_elements.append(overlap_value)
    overlap = create_matrix(overlap_matrix, matrix_elements, nbf)
    overlap = make_symmetric(overlap)
    logfile.close()
    overlap = np.matrix(overlap)

    return overlap


def create_matrix(matrix, elements, nbf):
    """ create lower triangular matrix from list of matrix elements
    indexed like so:
        [[0,0,0,...,0],
        [1,2,0,...,0],
        [3,4,5,...,0]]
    nbf is number of basis functions
    elements is a list of matrix elements indexed like above, e.g.
        [0,1,2,3,...]
    Gaussian prints every 5 columns, so the mod5 accounts for this
    """
    count = 0   # count is our index
    # fill the main block, leaving remainder for triangle fill
    for i in range(0, nbf-nbf%5, 5):
        matrix, count = triangle_fill(matrix, i, i+5, i, count, elements)
        matrix, count = block_fill(matrix, i+5, nbf, i, i+5, count, elements)
    # finish filling the last triangle bit
    matrix, count = triangle_fill(matrix, nbf-nbf%5, nbf, nbf-nbf%5, count, elements)

    return matrix

@jit
def make_symmetric(matrix):
    return matrix+matrix.T - np.diag(np.diag(matrix))

@jit
def triangle_fill(matrix, istrt, iend, jstrt, count, element):
    for i in range(istrt, iend):
        for j in range(jstrt, i+1):
            matrix[i, j] = element[count]
            count += 1
    return matrix, count

@jit
def block_fill(matrix, istrt, nbf, jstrt, jend, count, element):
    for i in range(istrt, nbf):
        for j in range(jstrt, jend):
            matrix[i, j] = element[count]
            count += 1
    return matrix, count


if __name__ == '__main__':
    TI(sys.argv[:])
