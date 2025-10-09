import numpy as np
import pytest

from cleedpy import search


@pytest.fixture
def dummy_config():
    # Create a mock config object with necessary attributes
    class Atom:
        def __init__(self, x, y, z):
            self.position = type("pos", (), {"x": x, "y": y, "z": z})

    class Config:
        overlayers = [Atom(1.0, 2.0, 3.0), Atom(4.0, 5.0, 6.0)]

        def model_dump(self):
            return {"key": "value"}

    return Config()


def test_init(monkeypatch):
    """Test the initialization of CleedSearchCoordinator."""

    # Mock the numpy.loadtxt function to avoid file I/O during tests
    monkeypatch.setattr("numpy.loadtxt", lambda x: np.array([[1, 2, 3, 4, 5]]))
    csc = search.CleedSearchCoordinator(
        "dummy_config", "dummy_phase_path", "dummy_experimental_iv_file"
    )

    assert csc.config == "dummy_config"
    assert csc.phase_path == "dummy_phase_path"
    assert csc.iteration == 0
    assert csc.current_rfactor == 0.0
    assert csc.x == []
    assert csc.correspondence == []
    assert csc.optimize_shift is True
    assert csc.optimal_shift == 0.0
    assert csc.largest_rfactor == 1.0
    assert csc.result is None


def test_set_search_parameters(monkeypatch, dummy_config):
    """Test the set_search_parameters method of CleedSearchCoordinator."""

    # Mock the numpy.loadtxt function to avoid file I/O during tests
    monkeypatch.setattr("numpy.loadtxt", lambda x: np.array([[1, 2, 3, 4, 5]]))

    csc = search.CleedSearchCoordinator(
        config=dummy_config,
        phase_path="dummy_phase_path",
        experimental_iv_file="dummy_experimental_iv_file",
    )
    csc.set_search_parameters(overlayer_atoms="xyz")
    assert len(csc.x) == 6  # 2 overlayers * 3 coordinates each
    assert len(csc.correspondence) == 6
    assert csc.correspondence == [
        "overlayers.0.position.x",
        "overlayers.0.position.y",
        "overlayers.0.position.z",
        "overlayers.1.position.x",
        "overlayers.1.position.y",
        "overlayers.1.position.z",
    ]
    assert csc.optimize_shift is True

    csc = search.CleedSearchCoordinator(
        config=dummy_config,
        phase_path="dummy_phase_path",
        experimental_iv_file="dummy_experimental_iv_file",
    )
    csc.set_search_parameters(overlayer_atoms=["x", "y"], optimize_shift=False)
    assert len(csc.x) == 2  # 2 overlayers * 1 coordinate each
    assert len(csc.correspondence) == 2
    assert csc.correspondence == ["overlayers.0.position.x", "overlayers.1.position.y"]
    assert csc.optimize_shift is False

    csc = search.CleedSearchCoordinator(
        config=dummy_config,
        phase_path="dummy_phase_path",
        experimental_iv_file="dummy_experimental_iv_file",
    )
    csc.set_search_parameters(overlayer_atoms=["xy", "xz"])
    assert len(csc.x) == 4  # 2 overlayers * 2 coordinates each
    assert len(csc.correspondence) == 4
    assert csc.correspondence == [
        "overlayers.0.position.x",
        "overlayers.0.position.y",
        "overlayers.1.position.x",
        "overlayers.1.position.z",
    ]


def test_set_params(monkeypatch, dummy_config):
    """Test the set_params method of CleedSearchCoordinator."""

    # Mock the numpy.loadtxt function to avoid file I/O during tests
    monkeypatch.setattr("numpy.loadtxt", lambda x: np.array([[1, 2, 3, 4, 5]]))

    csc = search.CleedSearchCoordinator(
        config=dummy_config,
        phase_path="dummy_phase_path",
        experimental_iv_file="dummy_experimental_iv_file",
    )
    csc.set_search_parameters(overlayer_atoms="xyz")
    new_values = np.array([10.0, 20.0, 30.0, 40.0, 50.0, 60.0])
    csc.set_params(new_values)

    assert dummy_config.overlayers[0].position.x == 10.0
    assert dummy_config.overlayers[0].position.y == 20.0
    assert dummy_config.overlayers[0].position.z == 30.0
    assert dummy_config.overlayers[1].position.x == 40.0
    assert dummy_config.overlayers[1].position.y == 50.0
    assert dummy_config.overlayers[1].position.z == 60.0
