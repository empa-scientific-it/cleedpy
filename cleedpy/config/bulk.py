from pathlib import Path
from typing import Any, Literal, NamedTuple, Tuple

from numpy.typing import ArrayLike
from pydantic import BaseModel


class UnitCellParameters(BaseModel):
    """Unit cell parameters for bulk calculations"""

    a1: tuple[float, float, float]
    a2: tuple[float, float, float]
    a3: tuple[float, float, float]


class SuperstructureMatrix(BaseModel):
    """Superstructure matrix for bulk calculations"""

    m1: tuple[int, int]
    m2: tuple[int, int]


class VibrationalDisplacementParameters(BaseModel):
    type: Literal["dmt"] | Literal["dr1"] | Literal["dr3"] | Literal["nd3"]


class DebyeTemperatureMass(NamedTuple):
    """Debye temperature and atomic mass for bulk calculations"""

    type: Literal["dmt"]
    debye_temperature: float
    atomic_mass: float
    temperature: float


class RadialMeanSquareDisplacement(NamedTuple):
    """Radial mean square displacement for bulk calculations"""

    type: Literal["dr1"]
    radial_mean_square_displacement: float


class MeanSquareDisplacements(NamedTuple):
    """Mean square displacements for bulk calculations- along the coordinate axes"""

    type: Literal["dr3"]
    x_mean_square_displacement: float
    y_mean_square_displacement: float
    z_mean_square_displacement: float


class MeanSquareDisplacementsND(NamedTuple):
    """Mean square displacements for bulk calculations along the coordinates with a non diagonal matrix"""

    type: Literal["nd3"]
    x_mean_square_displacement: float
    y_mean_square_displacement: float
    z_mean_square_displacement: float


VibrationalDisplacementParametersVariant = (
    DebyeTemperatureMass
    | RadialMeanSquareDisplacement
    | MeanSquareDisplacements
    | MeanSquareDisplacementsND
)


class Position(NamedTuple):
    x: float
    y: float
    z: float


class OptimiziationRangeParameters(NamedTuple):
    initial: float
    min: float
    max: float


class PositionOptimizationParameters(BaseModel):
    x: OptimiziationRangeParameters
    y: OptimiziationRangeParameters
    z: OptimiziationRangeParameters


class AtomParameters(BaseModel):
    """Atom parameters for overlayers"""

    phase_file: str | Path
    position: PositionOptimizationParameters | Position
    vibrational_displacement: VibrationalDisplacementParametersVariant


class EnergyRangeParameters(BaseModel):
    """Energy range parameters for bulk calculations"""

    initial: float
    final: float
    step: float


class NonGeometricalParameters(BaseModel):
    """Non geometrical parameters for bulk calculations"""

    unit_cell: UnitCellParameters
    superstructure_matrix: SuperstructureMatrix
    overlayers: list[AtomParameters]
    bulk_layers: list[AtomParameters]
    optical_potential: tuple[float, float]
    energy_range: EnergyRangeParameters
    polar_incidence_angle: float
    azimuthal_incidence_angle: float
    epsilon: float
    maximum_angular_momentum: int


class SearchRadiusParameters(BaseModel):
    """Search radius parameters for bulk calculations"""

    phase: str | Path
    radius: float


class RotationalSymmetryParameters(BaseModel):
    """Rotational symmetry parameters for bulk calculations"""

    degree: int
    x: float
    y: float


class SearchParameters(BaseModel):
    unit_cell: UnitCellParameters
    superstructure_matrix: SuperstructureMatrix
    overlayers: list[AtomParameters]
    search_radius: list[SearchRadiusParameters]
    z_range: tuple[float, float]
    parameter_variation: Literal["vertical"] | Literal["both"]
    rotational_symmetry: RotationalSymmetryParameters
    angle_search: bool
    polar_incidence_angle: float
    azimuthal_incidence_angle: float
    experimental_data: str | Path

    def model_post_init(self, __context: Any) -> None:
        return super().model_post_init(__context)

    def get_iv_curve(self) -> ArrayLike:
        """Get the IV curve from the experimental data"""
