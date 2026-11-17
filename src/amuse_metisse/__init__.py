"""
Interface for METISSE
"""
import os
from .interface import MetisseInterface
from .interface import Metisse
from .download import download_sample_metisse_tracks

if not os.path.exists(os.path.join(os.path.dirname(__file__), "data")):
    print("Downloading sample MESA tracks for METISSE")
    download_sample_metisse_tracks()

__all__ = ["MetisseInterface", "Metisse"]
