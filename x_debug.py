import subprocess
import aerosandbox as asb

_orig = subprocess.Popen.communicate
def _spy(self, input=None, timeout=None):
    print("=== STDIN TO XFOIL ===")
    print(input)
    print("======================")
    return _orig(self, input, timeout)
subprocess.Popen.communicate = _spy

xf = asb.XFoil(airfoil=asb.Airfoil("naca2412"), Re=1e6,
               working_directory="xfoil", verbose=True)
xf.alpha(0)