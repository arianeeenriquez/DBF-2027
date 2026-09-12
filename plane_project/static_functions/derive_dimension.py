import aerosandbox as asb
import aerosandbox.numpy as np

def run_dimensions(self):
	to_rad(self)
	wing_dimension(self)
	tail_dimension(self)
	
def wing_dimension(self):	
	#b: span, c_r root chord, c_t tip chord, MAC mean aero chord
	self.b = np.sqrt(self.S * self.AR)
	self.mean_geo_chord = self.S / self.b
	self.c_r = 2 / (self.taper + 1) * self.mean_geo_chord
	self.c_t = self.taper * self.c_r
	self.MAC = 2/3 * self.c_r * (1 + self.taper + self.taper**2)/ (1 + self.taper)
	self.yMAC = self.b * (1 + 2 * self.taper) / (6 + 6 * self.taper)
	self.proj_span = self.b/np.cos(self.dihedral_rad)

def tail_dimension(self):
	self.b_h = np.sqrt(self.S_h * self.AR_h)
	self.b_v = np.sqrt(self.S_v * self.AR_v)
	self.c_h = self.S_h / self.b_h
	self.c_v = self.S_v / self.b_v
	#rectangular tails

def to_rad(self):
	self.dihedral_rad = self.dihedral * np.pi / 180
	self.washout_rad = self.washout * np.pi / 180
	self.decalage_rad = self.decalage * np.pi / 180
