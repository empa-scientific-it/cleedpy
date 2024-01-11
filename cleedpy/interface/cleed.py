from ctypes import POINTER, Structure, c_char_p, c_double, c_int, cdll
from dataclasses import dataclass
from pathlib import Path
from typing import Type, Union
import math

import cleedpy
from cleedpy.interface.Matrix import MatPtr
from cleedpy.config import (SearchParameters, 
                            AtomParametersStructured, 
                            AtomParametersVariants,
                            Position, 
                            PositionOptimizationParameters, 
                            VibrationalDisplacementParametersVariant,
                            DebyeTemperatureMass,
                            MeanSquareDisplacements,
                            MeanSquareDisplacementsND,
                            RadialMeanSquareDisplacement,
                            )

from cleedpy.physics import constants


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
        ("eng_v", c_double),
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


def generate_energies(input: list[float]) -> list[c_double]:
    return [c_double(i) for i in input]



def convert_energy_loop_variables(input: SearchParameters) -> EnergyLoopVariables:
    """
    This corresponds to the following c function: inp_rdpar
    """
    return EnergyLoopVariables(
        eng_r=c_double(input.energy_range.initial / constants.HART),
        eng_i=c_double(input.energy_range.final / constants.HART),
        eng_v=c_double(input.energy_range.step / constants.HART),
        vr=c_double(input.optical_potential[0]),
        vi_prev=c_double(input.optical_potential[1]),
        vi_exp=c_double(input.optical_potential[1]),
        theta=c_double(math.radians(input.polar_incidence_angle)),
        phi=c_double(math.radians(input.azimuthal_incidence_angle)),
        k_in=(c_double * 4)(0, 0, 0, 0),
        epsilon=c_double(input.epsilon),
        l_max=c_int(input.maximum_angular_momentum),
        p_t1=POINTER(MatPtr),
    )



def generate_single_atom_structure(input: AtomParametersVariants) -> Atom:
    """
    Generate a single atom for cleed from the input atom parameters
    """
    match input:
        case AtomParametersStructured(phase_file, Position(x, y, z), vibrational_displacement):
            a = Atom()
            a.layer = 0
            a.type = 0
            a.t_type = 0
            a.pos = (x, y, z, 0)
            a.dwf

        case AtomParametersStructured(phase_file, PositionOptimizationParameters(), vibrational_displacement):
            raise ValueError("Cleed needs a position for each atom, not an optimization range")
        case _:
            raise ValueError("Not implemented")

def generate_atom_structure(input: list[AtomParametersVariants]) -> list[Atom]:
    """
    Generate a list of atoms for cleed from the input list of atom parameters
    It also sorts the atoms by the z-coordinate
    """
    return [Atom() for i in input]

def generate_crystal_structure(input: AtomParametersVariants) -> Crystal:
    match input:
        case AtomParametersStructured(phase_file, Position(x, y, z), vibrational_displacement):
            pass
        case AtomParametersStructured(phase_file, PositionOptimizationParameters(), vibrational_displacement):
            raise ValueError("Cleed needs a position for each atom, not an optimization range")
        case _:
            raise ValueError("Not implemented")

def call_cleed():
    path = Path(cleedpy.__file__).parent / "cleed" /  "build" / "lib" / "libcleed.so"
    lib = cdll.LoadLibrary(str(path))
    lib.leed


if __name__ == "__main__":
    call_cleed()
