import math
from sympy import *
import os.path


# Basic cross product formula, a and b should be three element tuples
def crossProduct(a, b):
    return a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]


# Basic 3 element dot product. a and b should be three element tuples
def dotProduct(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


# This function scales the output values f1, f2, f3 (contained in tuple r) such that they output a vector within the
# Brillouin zone. NOTE: this does not actually compute the brillouin zone, it merely uses an approximation by
# checking in 26 distinct directions
def scale(t, r, b1, b2, b3):
    s = 1
    n = (-1, 0, 1)
    for x in range(3):
        for y in range(3):
            for z in range(3):

                # This if check is necessary so the component in the direction (0,0,0) is not calculated
                if x == 1 and y == 1 and z == 1:
                    continue
                dir = (n[x] * b1[0] + n[y] * b2[0] + n[z] * b3[0],
                       n[x] * b1[1] + n[y] * b2[1] + n[z] * b3[1],
                       n[x] * b1[2] + n[y] * b2[2] + n[z] * b3[2])
                comp = component(dir, t)
                if comp > 0.5:
                    if s > (1 / (2 * comp)):
                        # print("max dir: ", dir)
                        s = 1 / (2 * comp)
    # print("s: ", s)
    return tuple(s * x for x in r)


# This function returns the amount of a that b's projection on a covers (unit-less)
def component(a, b):
    return dotProduct(a, b) / dotProduct(a, a)


# This function returns the magnitude of the vector a
def magnitude(a):
    return math.sqrt((a[0]) ** 2 + (a[1]) ** 2 + (a[2]) ** 2)


# Determines if the lattice vectors are linearly independent
def latticeValid(a, b, c):
    M = Matrix([[a[0], a[1], a[2]], [b[0], b[1], b[2]], [c[0], c[1], c[2]]])
    return M.det() != 0


# This function returns the lattice vectors given the file exists within the same directory and the lattice vectors
# are linearly independent
def form(filename):
    # Check to ensure geometry.in file exists
    if not os.path.exists(filename):
        print("ERROR: Given geometry.in file does not exist")
        return

    # Reads in lattice vector information from geometry.in file
    with open(filename, "r") as f:
        lv = []
        for ln in f:
            if ln.startswith("lattice_vector"):
                lv.append(ln[15:])

    # Check to ensure exactly 3 lattice vectors were found
    if len(lv) != 3:
        print("ERROR: Wrong number of lattice vectors in input geometry.in file")
        return

    # Formats lattice vectors into separate variables
    aS = lv[0].split()
    bS = lv[1].split()
    cS = lv[2].split()
    a = float(aS[0]), float(aS[1]), float(aS[2])
    b = float(bS[0]), float(bS[1]), float(bS[2])
    c = float(cS[0]), float(cS[1]), float(cS[2])

    if not latticeValid(a, b, c):
        print("ERROR: Lattice vectors in geometry.in file linearly dependent")
        return
    return a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]


# Main function, input filepath of geometry.in file as well as the vector (t1,t2,t3) expressed in lattice real-space
# of the desired band direction
# inputs:
# filename: path to geometry.in file
# target: 3 element tuple values either cartesian length values or units of lattice vectors based on lattice_defined var
# lattice_defined: true or false, if true target is interpreted as multiples of lattice vectors
# if false target is interpreted as cartesian values (in angstroms)
# debug determines whether debug information is printed
# mag: if not input or 0 the result is scaled to be in an approximation of the brillouin zone, otherwise the result will be scaled
# to have a magnitude as determined by this parameter
def generate(filename, target, lattice_defined=True, debug=False, mag=0):
    # The form function returns a 9-element tuple containing the lattice vectors from the given geometry.in file
    lattice = form(filename)
    a = (lattice[0], lattice[1], lattice[2])
    b = (lattice[3], lattice[4], lattice[5])
    c = (lattice[6], lattice[7], lattice[8])
    if debug and lattice_defined:
        print("Real space lattice vectors")
        print(a)
        print(b)
        print(c)

    # The vector (A1,A2,A3) is the target that the result will be parallel to
    if lattice_defined:
        A1 = (target[0] * a[0] + target[1] * b[0] + target[2] * c[0])
        A2 = (target[0] * a[1] + target[1] * b[1] + target[2] * c[1])
        A3 = (target[0] * a[2] + target[1] * b[2] + target[2] * c[2])
    else:
        A1 = target[0]
        A2 = target[1]
        A3 = target[2]
    target_vector = (A1, A2, A3)

    # B1, B2, and B3 define the three reciprocal unit cell vectors. V is the volume of the unit cell.
    V = abs(dotProduct(a, crossProduct(b, c)))
    B1 = tuple((2 * math.pi / V) * x for x in crossProduct(b, c))
    B2 = tuple((2 * math.pi / V) * x for x in crossProduct(c, a))
    B3 = tuple((2 * math.pi / V) * x for x in crossProduct(a, b))
    if debug:
        print("target: ", target_vector)
        print("reciprocal lattice vectors:")
        print(B1)
        print(B2)
        print(B3)

    # M defines the augmented matrix consisting of B1, B2, and B3 as columns with (A1,A2,A3) as the augmented column.
    # We know this is full rank by definition. Thus, row reducing produces the desired result
    M = Matrix([[B1[0], B2[0], B3[0], A1], [B1[1], B2[1], B3[1], A2], [B1[2], B2[2], B3[2], A3]])
    M_rref = M.rref()
    result = (M_rref[0][3], M_rref[0][7], M_rref[0][11])

    final = scale(target_vector, result, B1, B2, B3)
    if debug:
        print("initial f vals: ", result)
        r_x = B1[0] * final[0] + B2[0] * final[1] + B3[0] * final[2]
        r_y = B1[1] * final[0] + B2[1] * final[1] + B3[1] * final[2]
        r_z = B1[2] * final[0] + B2[2] * final[1] + B3[2] * final[2]
        print(A1)
        print(A2)
        print(A3)
        print(r_x)
        print(r_y)
        print(r_z)
    if mag != 0:
        leng = math.sqrt(final[0] ** 2 + final[1] ** 2 + final[2] ** 2)
        final = (mag * final[0] / leng, mag * final[1] / leng, mag * final[2] / leng)
    if debug:
        print("Result:")
    print(str(final[0]) + " " + str(final[1]) + " " + str(final[2]))
    return final


base = "."
f = base + "/geometry.in"
generate(f, [0, 1, 0], lattice_defined=True, debug=False)
