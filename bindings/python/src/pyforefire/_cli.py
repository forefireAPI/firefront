"""Entry point for the bundled ``forefire`` command line interpreter.

The interpreter is a native binary shipped inside this package rather than as
a wheel script. Wheel repair tools rewrite a binary's library references
relative to the binary's own location, and pip installs ``.data/scripts/``
and the package into different directories, so a binary placed in the scripts
directory cannot reach the vendored libraries once installed. Keeping it in
the package and exec'ing it from here avoids the problem on every platform.
"""

import os
import sys
from pathlib import Path

#: Name of the native executable installed next to this module by CMake.
_EXECUTABLE = "forefire.exe" if sys.platform == "win32" else "forefire"


def executable_path() -> Path:
    """Return the path to the bundled interpreter."""
    return Path(__file__).resolve().parent / _EXECUTABLE


def main() -> None:
    """Replace this process with the interpreter, forwarding all arguments."""
    binary = executable_path()
    if not binary.exists():
        raise SystemExit(
            f"the forefire executable is missing from {binary.parent}. "
            "This build of the package does not ship the command line "
            "interpreter; build from source to get it."
        )
    os.execv(str(binary), ["forefire", *sys.argv[1:]])
