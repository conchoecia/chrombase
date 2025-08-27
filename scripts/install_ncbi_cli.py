"""Install the latest NCBI datasets and dataformat CLI tools.

This script fetches the most recent release of the `ncbi/datasets` project
from GitHub, downloads the command line package appropriate for the user's
platform, and extracts the ``datasets`` and ``dataformat`` binaries into the
repository's ``bin/`` directory.

The GitHub release provides a single zip file containing the binaries for
each supported operating system/architecture pair.  We select the correct
asset based on ``platform.system()`` and ``platform.machine()``.  After
downloading the archive we extract the two executables and mark them
executable.

Running this script multiple times is safe – existing binaries will simply be
overwritten.
"""

from __future__ import annotations

import json
import os
import platform
import stat
import tempfile
import urllib.request
import zipfile


GITHUB_API = "https://api.github.com/repos/ncbi/datasets/releases/latest"


def _detect_platform() -> tuple[str, str]:
    """Return the GitHub asset prefix for the current OS and architecture."""

    system = platform.system().lower()
    machine = platform.machine().lower()

    if system.startswith("linux"):
        os_name = "linux"
    elif system.startswith("darwin"):
        os_name = "darwin"
    elif system.startswith("win"):
        os_name = "windows"
    else:
        raise RuntimeError(f"Unsupported operating system: {system}")

    arch_map = {
        "x86_64": "amd64",
        "amd64": "amd64",
        "aarch64": "arm64",
        "arm64": "arm64",
        "armv7l": "arm",
        "armv6l": "arm",
    }
    arch = arch_map.get(machine, machine)
    return os_name, arch


def _get_download_url(os_name: str, arch: str) -> str:
    """Fetch the download URL for the zip asset matching ``os_name``/``arch``."""

    asset_name = f"{os_name}-{arch}.cli.package.zip"
    with urllib.request.urlopen(GITHUB_API) as resp:
        release = json.load(resp)

    for asset in release.get("assets", []):
        if asset.get("name") == asset_name:
            return asset["browser_download_url"]

    raise RuntimeError(f"No release asset named {asset_name} found")


def install_tools() -> None:
    """Download and install the datasets and dataformat binaries."""

    os_name, arch = _detect_platform()
    url = _get_download_url(os_name, arch)

    repo_root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    bin_dir = os.path.join(repo_root, "bin")
    os.makedirs(bin_dir, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        zip_path = os.path.join(tmpdir, "cli.zip")
        urllib.request.urlretrieve(url, zip_path)
        with zipfile.ZipFile(zip_path) as zf:
            for member in ("datasets", "dataformat"):
                zf.extract(member, bin_dir)
                exe_path = os.path.join(bin_dir, member)
                # Ensure the files are executable
                st = os.stat(exe_path)
                os.chmod(
                    exe_path,
                    st.st_mode
                    | stat.S_IXUSR
                    | stat.S_IXGRP
                    | stat.S_IXOTH,
                )


if __name__ == "__main__":
    install_tools()

