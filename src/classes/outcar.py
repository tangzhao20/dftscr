import subprocess
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
        str_out += " toten = " + str(self.toten) + "\n"
        str_out += " mag = " + str(self.mag) + "\n"
        return str_out

#######################################################################

    def read_vasp(self, file_name="OUTCAR"):
        cmd = "grep 'energy  without entropy=' " + file_name + " | tail -n 1 | awk '{print $NF}'"
        self.toten = float(subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.strip())

        cmd = "grep 'number of electron' " + file_name + " | tail -n 1 | awk '{print $NF}'"
        self.mag = float(subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.strip())
