import sys
import os
import subprocess
import xml.etree.ElementTree as ET
from load_data import load_constant


class Outcar:
    """
    This class reads data from calculations, such as total energy and magnetization.

    Attributes:
        toten (float): total energy in eV.
        mag (float): magnetization in mu_B.
    """

    def __init__(self):
        self.toten = 0.0
        self.mag = 0.0

    def __str__(self):
        str_out = "OUTCAR:\n"
        str_out += f" toten = {self.toten} eV\n"
        str_out += f" mag = {self.mag} mu_B\n"
        return str_out

#######################################################################

    def read_vasp(self, file_name="OUTCAR"):
        cmd = "grep 'energy  without entropy=' " + file_name + " | tail -n 1 | awk '{print $NF}'"
        self.toten = float(subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.strip())

        cmd = "grep 'number of electron' " + file_name + " | tail -n 1 | awk '{print $NF}'"
        self.mag = float(subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.strip())

    def read_xml(self, filename=""):

        if filename == "":
            # find a .xml file
            files = os.listdir()
            for f in files:
                if f.endswith('.xml'):
                    filename = f
                    break
        if filename == "":
            print("Error: .xml file is not found")
            sys.exit()

        Ha = load_constant("rydberg")*2.0

        tree = ET.parse(filename)
        output = tree.getroot().find("output")

        toten = float(output.find("total_energy").find("etot").text)
        demet = float(output.find("total_energy").find("demet").text)
        self.toten = (toten - demet) * Ha

        self.mag = float(output.find("magnetization").find("total").text)
