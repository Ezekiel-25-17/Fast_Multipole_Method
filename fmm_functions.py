import numpy as np
from quadtree import QuadTree, build_tree
from scipy.special import binom

# Functions Explained:
## potential_fmm: fmm evaluation of all-to-all potential
## potential_fmm_tree: fmm evaluation of all-to-all potential with a prebuilt tree
## potential_ds: direct sum calculation of all-to-all potential
## potential_dds: direct sum calculation of all-to-all potential from different sources


def potential_fmm(particles, bbox=None, tree_threshold=None, nterms=5, boundary='wall'):
    tree = build_tree(particles, tree_threshold, bbox=bbox, boundary=boundary)
    outer_multipole_exp(tree.root, nterms)
    tree.root.inner = np.zeros((nterms + 1), dtype=complex)
    any(inner_exp(child) for child in tree.root)


def potential_fmm_tree(tree, nterms=5):
    outer_multipole_exp(tree.root, nterms)
    tree.root.inner = np.zeros((nterms + 1), dtype=complex)
    any(inner_exp(child) for child in tree.root)


def potential_ds(particles):
    phi = np.zeros((len(particles),))

    for i, particle in enumerate(particles):
        for source in (particles[:i] + particles[i+1:]):
            r = distance(particle.pos, source.pos)
            particle.phi -= particle.q*np.log(r)
        phi[i] = particle.phi

    return phi


def potential_dds(particles, sources):
    for i, particle in enumerate(particles):
        for source in sources:
            r = distance(particle.pos, source.pos)
            particle.phi -= particle.q*np.log(r)


def distance(point_1, point_2):
    distance = np.sqrt((point_1[0] - point_2[0])**2 + (point_1[1] - point_2[1])**2)
    return distance


def multipole(particles, center=(0,0), nterms=5):
    coeffs = np.empty(nterms + 1, dtype=complex)
    coeffs[0] = sum(p.q for p in particles)
    coeffs[1:] = [sum([-p.q*complex(p.x - center[0], p.y - center[1])**k/k
                  for p in particles]) for k in range(1, nterms+1)]

    return coeffs


def shift_multipole_exp(coeffs, z0):
    shift = np.empty_like(coeffs)
    shift[0] = coeffs[0]
    shift[1:] = [sum([coeffs[k]*z0**(l - k)*binom(l-1, k-1) - (coeffs[0]*z0**l)/l
                  for k in range(1, l)]) for l in range(1, len(coeffs))]

    return shift


def shift_taylor_exp(coeffs, z0):
    shift = np.empty_like(coeffs)
    shift = [sum([coeffs[k]*binom(k,l)*(-z0)**(k-l)
              for k in range(l,len(coeffs))])
              for l in range(len(coeffs))]
    return shift


def outer_multipole_exp(tnode, nterms):
    if tnode.is_leaf():
        tnode.outer = multipole(tnode.get_points(), center=tnode.center, nterms=nterms)
    else:
        tnode.outer = np.zeros((nterms + 1), dtype=complex)
        for child in tnode:
            outer_multipole_exp(child, nterms)
            z0 = complex(*child.center) - complex(*tnode.center)
            tnode.outer += shift_multipole_exp(child.outer, z0)


def convert_outer_to_inner(coeffs, z0):
    inner = np.empty_like(coeffs)
    inner[0] = (sum([(coeffs[k]/z0**k)*(-1)**k for k in range(1, len(coeffs))]) + coeffs[0]*np.log(-z0))
    inner[1:] = [(1/z0**l)*sum([(coeffs[k]/z0**k)*binom(l+k-1, k-1)*(-1)**k
                 for k in range(1, len(coeffs))]) - coeffs[0]/((z0**l)*l)
                 for l in range(1, len(coeffs))]
    return inner


def inner_exp(tnode):
    z0 = complex(*tnode.parents.center) - complex(*tnode.center)
    tnode.inner = shift_taylor_exp(tnode.parents.inner, z0)
    for tin in tnode.interaction_set():
        z0 = complex(*tin.center) - complex(*tnode.center)
        tnode.inner += convert_outer_to_inner(tin.outer, z0)

    if tnode.is_leaf():
        z0, coeffs = complex(*tnode.center), tnode.inner

        for p in tnode.get_points():
            z = complex(*p.pos)
            p.phi -= np.real(np.polyval(coeffs[::-1], z-z0))

        for nn in tnode.nearest_neighbors:
            potential_dds(tnode.get_points(), nn.get_points())
        _ = potential_ds(tnode.get_points())
    else:
        for child in tnode:
            inner_exp(child)
