from ctypes import POINTER, Structure, c_char_p, c_double, c_int, cdll
from dataclasses import dataclass
from pathlib import Path

import cleedpy
from cleedpy.interface.Matrix import MatPtr


class Atom(Structure):
    _fields_ = [
        ("layer", c_int),
        ("type", c_int),
        ("t_type", c_int),
        ("pos", c_double * 4),
        ("dwf", c_double),
    ]


class Layer(Structure):
    _fields_ = [
        ("no_of_layer", c_int),
        ("bulk_over", c_int),
        ("periodic", c_int),
        ("natoms", c_int),
        ("a_lat", c_double * 5),
        ("rel_area", c_double),
        ("reg_shift", c_double * 4),
        ("vec_from_last", c_double * 4),
        ("vec_to_next", c_double * 4),
        ("atoms", POINTER(Atom)),
    ]


class Crystal(Structure):
    _fields_ = [
        ("vr", c_double),
        ("vi", c_double),
        ("temp", c_double),
        ("n_rot", c_int),
        ("rot_axis", c_double * 4),
        ("n_mir", c_int),
        ("m_plane", POINTER(c_double)),
        ("alpha", POINTER(c_double)),
        ("symmetry", c_int),
        ("a", c_double * 5),
        ("a_1", c_double * 5),
        ("area", c_double),
        ("m_trans", c_double * 5),
        ("m_super", c_double * 5),
        ("m_recip", c_double * 5),
        ("b", c_double * 5),
        ("b_1", c_double * 5),
        ("rel_area_sup", c_double),
        ("nlayers", c_int),
        ("layers", POINTER(Layer)),
        ("dmin", c_double),
        ("natoms", c_int),
        ("ntypes", c_int),
        ("comments", POINTER(c_char_p)),
    ]


class PhaseShifts(Structure):
    _fields_ = [
        ("lmax", c_int),
        ("neng", c_int),
        ("t_type", c_int),
        ("eng_max", c_double),
        ("eng_min", c_double),
        ("energy", POINTER(c_double)),
        ("pshift", POINTER(c_double)),
        ("dr", c_double * 4),
        ("input_file", c_char_p),
    ]


class LeedResult(Structure):
    _fields_ = [("x", c_double), ("y", c_double), ("result", c_double)]


class EnergyLoopVariables(Structure):
    _fields_ = [
        ("eng_r", c_double),
        ("eng_i", c_double),
        ("vr", c_double),
        ("vi_prev", c_double),
        ("vi_exp", c_double),
        ("theta", c_double),
        ("phi", c_double),
        ("k_in", c_double * 4),
        ("epsilon", c_double),
        ("l_max", c_int),
        ("p_t1", POINTER(MatPtr)),
    ]


@dataclass
class CleedInputs:
    bulk: Crystal
    overlayers: Crystal
    energies: list[float]
    phase_shifts: PhaseShifts
    params: EnergyLoopVariables


def prepare_cleed_input() -> CleedInputs:
    pass


def call_cleed():
    path = Path(cleedpy.__file__).parent / "cleed" \
    / "build" / "lib" / "libcleed.so"
    lib = cdll.LoadLibrary(str(path))
    lib.leed


if __name__ == "__main__":
    call_cleed()
