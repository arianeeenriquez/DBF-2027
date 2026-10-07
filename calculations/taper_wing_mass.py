import aerosandbox as asb

airfoil = asb.Airfoil("sd7032")
c_r = 0.48
c_t = 0.48
b = 1.767


def const_quarter_chord(rc, tc):
    return rc / 4 - tc / 4

wing = asb.Wing(
    name="wing",
    symmetric=True,
    xsecs=[
        asb.WingXSec(
            xyz_le=[0, 0, 0],  # UPDATED per Taras/AVL documentation
            chord=c_r,
            airfoil=airfoil,
          ) 
        ,
        asb.WingXSec(
            xyz_le=[
                const_quarter_chord(c_r, c_t),
                b / 2,
                0
            ],
            chord=c_t,
            airfoil=airfoil
        )
    ]
)

print(wing.volume())
