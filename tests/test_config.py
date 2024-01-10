from cleedpy.config.bulk import AtomParametersWrapper, SearchParameters
import yaml
import pkg_resources
import pytest
import cleedpy

@pytest.fixture
def yaml_data() -> str:
    path =  pkg_resources.resource_filename("tests", "resources/test_config.yaml")
    with open(path, "r") as file:
        return file.read()

def test_overlayer_parameters(yaml_data: str):
    atom = yaml.load(yaml_data, Loader=yaml.FullLoader)
    atom_parameters = SearchParameters.model_validate(atom)
    assert atom_parameters.overlayers[0].phase_file == "a"