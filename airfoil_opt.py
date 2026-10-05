import neuralfoil as nf
import aerosandbox as asb
import os

PATH_TO_AIRFOILS = "airfoil"

##### CONSTANTS ####
Re = 1_250_000 # TODO: make the actual number at some point
alpha = 0
n_crit = 6
mach = 0.05

def analyze_airfoil(airfoil_name : str):
    coordinates = asb.geometry.airfoil.airfoil_families.get_file_coordinates(airfoil_name)
    airfoil = asb.geometry.Airfoil(coordinates=coordinates)

    aero = nf.get_aero_from_dat_file(
        airfoil_name,
        alpha=alpha,
        Re=Re,
    )
    # print(aero)
    Cd = aero["CD"]
    Cm = aero["CM"]
    Cl = aero["CL"]
    area = airfoil.area()
    return score(Cd, Cm, Cl, area) if aero["analysis_confidence"] > 0.8 else -float("inf")

    # TODO: xfoil borken
    # xf = XFoil(
    #     airfoil=Airfoil(airfoil_name),
    #     Re=1e6,
    #     # xfoil_command="/Users/alex/Desktop/DBF/Xfoil-for-Mac/bin/xfoil",
    #     working_directory="xfoil out",
    #     verbose=True,
    #     xfoil_repanel=True,
    # )
    # results = xf.alpha(alpha)
    # print(results)
    # return(results)

    # return score(Cd, Cm, Cl, area)


def main():
    airfoils = [os.path.join(PATH_TO_AIRFOILS, dat_file) for dat_file in os.listdir(PATH_TO_AIRFOILS)]
    print(airfoils)

if __name__ == "__main__":
    main()