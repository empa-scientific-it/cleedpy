import os
import pathlib as pl
import platform
import subprocess
import sys

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name, source_dir=""):
        Extension.__init__(self, name, sources=[])
        self.source_dir = pl.Path(source_dir).absolute()


class CMakeBuild(build_ext):
    def run(self):
        try:
            subprocess.check_output(["cmake", "--version"])
        except subprocess.CalledProcessError as err:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            ) from err

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext: CMakeExtension):
        ext_dir = pl.Path(self.get_ext_fullpath(ext.name)).parent.absolute()

        cmake_args = ["-DPYTHON_EXECUTABLE=" + sys.executable]

        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]

        if platform.system() == "Windows":
            cmake_args += [
                f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={ext_dir}"
            ]
            if sys.maxsize > 2**32:
                cmake_args += ["-A", "x64"]
            build_args += ["--", "/m"]
        else:
            build_args += ["--", "-j2"]

        env = os.environ.copy()
        env["CXXFLAGS"] = '{} -DVERSION_INFO="{}"'.format(
            env.get("CXXFLAGS", ""), self.distribution.get_version()
        )

        build_dir = pl.Path(self.build_temp)
        build_dir.mkdir(parents=True, exist_ok=True)

        try:
            subprocess.check_call(
                ["cmake", "-S", ext.source_dir.as_posix(), "-B", "."] + cmake_args,
                cwd=build_dir,
                env=env,
            )
            subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=build_dir)
        except subprocess.CalledProcessError as err:
            raise RuntimeError(f"CMake error: {err}") from err


setup(
    name="cleedpy",
    ext_modules=[
        CMakeExtension(
            "cleed", (pl.Path(__file__).parent / "cleedpy" / "cleed").as_posix()
        )
    ],
    cmdclass=dict(build_ext=CMakeBuild),
)
