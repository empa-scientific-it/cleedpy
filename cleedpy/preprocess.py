import numpy as np

from .interface.cleed import Crystal
from .physics.constants import ANGSTROM_TO_BOHR, BOHR_TO_ANGSTROM

GEO_TOLERANCE = 0.0001


def transform_cells(parameters):
    """Converto to right handed coordinate system."""
    a1 = np.array(parameters.unit_cell.a1)
    a2 = np.array(parameters.unit_cell.a2)
    m_transformation = np.array([[1.0, 0], [0, 1.0]])

    # Check the angle between a1 and a2 (1: must be >= pi/2, i.e. a1*a2 < 0).
    # if not: take a1 and -a2 as basic vectors and modify transformation matrix.
    if np.dot(a1, a2) > 0.0:
        a2 *= -1.0
        m_transformation[1] *= -1.0

    # Check the angle between a1 and a2 (2: must be < pi, i.e. a1*a2 > 0)
    # if not: exchange a1 and a2 and modify transformation matrix.
    if np.cross(a1, a2)[2] < 0.0:
        a1, a2 = a2, a1
        m_transformation[0, 1] = m_transformation[1, 0]

    parameters.unit_cell.a1, parameters.unit_cell.a2 = tuple(a1), tuple(a2)

    return m_transformation


def extract_bulk_parameters(parameters, transformation_matrix):
    bulk = Crystal()

    bulk.vr, bulk.vi = parameters.optical_potential
    bulk.temp = parameters.sample_temperature

    def to_matrix(v1, v2):
        """Converts two vectors to a matrix."""
        return np.array([v1, v2])

    # Unit cell.
    a = (
        to_matrix(parameters.unit_cell.a1[:2], parameters.unit_cell.a2[:2])
        * ANGSTROM_TO_BOHR
    )
    bulk.a = (0.0, a[0, 0], a[1, 0], a[0, 1], a[1, 1])

    # Reciprocal lattice vectors.
    a_recip = np.linalg.inv(a).T * 2.0 * np.pi
    print(a_recip / BOHR_TO_ANGSTROM)
    bulk.a_1 = (0.0, a_recip[0, 0], a_recip[0, 1], a_recip[1, 0], a_recip[1, 1])

    # Superstructure matrix
    superstructure_matrix = to_matrix(
        parameters.superstructure_matrix.m1, parameters.superstructure_matrix.m2
    )
    superstructure_matrix @= transformation_matrix
    bulk.m_super = (
        0.0,
        superstructure_matrix[0, 0],
        superstructure_matrix[1, 0],
        superstructure_matrix[0, 1],
        superstructure_matrix[1, 1],
    )
    bulk.m_trans = (
        0.0,
        transformation_matrix[0, 0],
        transformation_matrix[1, 0],
        transformation_matrix[0, 1],
        transformation_matrix[1, 1],
    )

    # Reciprocal superstructure matrix.
    reciprocal_matrix = np.linalg.inv(superstructure_matrix).T
    bulk.m_recip = (
        0.0,
        reciprocal_matrix[0, 0],
        reciprocal_matrix[0, 1],
        reciprocal_matrix[1, 0],
        reciprocal_matrix[1, 1],
    )

    return bulk
