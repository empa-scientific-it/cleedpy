import numpy as np
from scipy import optimize

from . import config, physics, rfactor
from .interface.cleed import call_cleed


def cleed_result_to_iv(result) -> np.ndarray:
    """Convert CLEED result to IV array."""
    iv = []
    n_beams = result.n_beams
    n_energies = result.n_energies
    for i in range(n_energies):
        for j in range(n_beams):
            iv.append(
                [
                    result.beam_index1[j],
                    result.beam_index2[j],
                    0,
                    result.energies[i] * physics.constants.HART,
                    result.iv_curves[i * n_beams + j],
                ]
            )
    return np.array(iv)


class CleedSearchCoordinator:
    def __init__(
        self,
        config: config.InputParameters,
        phase_path: str,
        experimental_iv_file: str,
    ) -> None:
        self.config = config
        self.phase_path = phase_path
        self.iteration = 0
        self.current_rfactor = 0.0
        self.experimental_iv = np.loadtxt(experimental_iv_file)
        self.x = []
        self.correspondence = []
        self.optimize_shift = True
        self.optimal_shift = 0.0
        self.largest_rfactor = 1.0
        self.result = None

    def start_optimization(self, method: str = "Nelder-Mead") -> None:
        """Start the optimization process."""
        x_init = np.array(self.x)
        self.result = optimize.minimize(
            self.function_to_minimize,
            x_init,
            method=method,
            options={"disp": True},
            callback=self.print_things_to_file,
            tol=5e-4,
        )

    def print_things_to_file(self, intermediate_result):
        """Print the current parameters and R-factor to a file."""

        with open("optimization_history.log", "a") as fobj:
            fobj.write(f"Iteration {self.iteration} ")
            fobj.write(f"Shift: {self.optimal_shift} ")
            fobj.write(f"R-factor: {self.current_rfactor}\n")

    def set_search_parameters(
        self, overlayer_atoms: list[str] | str = "xyz", optimize_shift: bool = True
    ) -> None:
        """Set the search parameters."""

        if isinstance(overlayer_atoms, str):
            overlayer_atoms = [overlayer_atoms] * len(self.config.overlayers)

        for i, atom in enumerate(self.config.overlayers):
            if "x" in overlayer_atoms[i]:
                self.x.append(atom.position.x)
                self.correspondence.append(f"overlayers.{i}.position.x")
            if "y" in overlayer_atoms[i]:
                self.x.append(atom.position.y)
                self.correspondence.append(f"overlayers.{i}.position.y")
            if "z" in overlayer_atoms[i]:
                self.x.append(atom.position.z)
                self.correspondence.append(f"overlayers.{i}.position.z")

        self.optimize_shift = optimize_shift

    def set_params(self, x: np.typing.ArrayLike) -> None:
        """Set the parameters in the configuration object."""
        if self.config is None:
            msg = "Configuration object is not connected."
            raise ValueError(msg)

        for i, value in enumerate(x):
            path = self.correspondence[i].split(".")
            obj = self.config
            for p in path[:-1]:
                if p.isdigit():
                    obj = obj[int(p)]
                else:
                    obj = getattr(obj, p)
            setattr(obj, path[-1], value)

    def function_to_minimize(self, x: np.typing.ArrayLike) -> float:
        self.iteration += 1

        # Update parameters in the config object
        self.set_params(x)

        # First compute the geometrical R-factor, as it is fast.
        geometrical_r = rfactor.compute_geometrical_rfactor(self.config)

        # No need to continue if the geometrical R-factor is too large.
        if geometrical_r > 1.0:
            return geometrical_r + self.largest_rfactor

        # Write the current parameters to a temporary file.
        old_format = config.OLD_FORMAT_TEMPLATE.render(**self.config.model_dump())
        with open("current_parameters.inp", "w") as fobj:
            fobj.write(old_format)

        # Call CLEED with the current parameters.
        result = call_cleed(
            "current_parameters.inp", "current_parameters.inp", self.phase_path
        )

        self.theoretical_iv = cleed_result_to_iv(result)

        # Optimize the shift if requested.
        if self.optimize_shift:
            shift_opt_result = optimize.minimize_scalar(
                self.function_to_minimize_shift,
                bounds=(-10.0, 10.0),
                method="bounded",
            )
            self.optimal_shift = shift_opt_result.x

        iv_r = rfactor.compute_rfactor(
            theoretical_iv=self.theoretical_iv,
            experimental_iv=self.experimental_iv,
            shift=self.optimal_shift,
            rfactor_type="pendry",
        )

        self.current_rfactor = iv_r + geometrical_r

        return self.current_rfactor

    def function_to_minimize_shift(self, shift: float) -> float:
        """Function to minimize the shift between theoretical and experimental IV curves."""
        r = rfactor.compute_rfactor(
            theoretical_iv=self.theoretical_iv,
            experimental_iv=self.experimental_iv,
            shift=shift,
            rfactor_type="pendry",
        )
        return r
