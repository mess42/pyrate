import ctypes
'''
typedef struct
{
 
double x, y, z;     /* the coordinates */
double l, m, n;     /* the ray direction cosines */
double ln, mn, nn;  /* the surface normals */
   double path;        /* the optical path change */
   double sag1, sag2;  /* the sag and alternate hyperhemispheric sag */
double index, dndx, dndy, dndz; /* for GRIN surfaces only */
   double rel_surf_tran; /* for relative surface transmission data, if any */
   double udreserved1, udreserved2, udreserved3, udreserved4; /* for future expansion */
   char string[20];    /* for returning string data */
 
}USER_DATA;
'''
class USER_DATA(ctypes.Structure):
	_fields_ = [("x", ctypes.c_double), ("y", ctypes.c_double), ("z", ctypes.c_double),
			("l", ctypes.c_double), ("m", ctypes.c_double), ("n", ctypes.c_double),
                        ("ln", ctypes.c_double), ("mn", ctypes.c_double), ("nn", ctypes.c_double),
			("path", ctypes.c_double), ("sag1", ctypes.c_double), ("sag2", ctypes.c_double),
			("index", ctypes.c_double), ("dndx", ctypes.c_double), ("dndy", ctypes.c_double), ("dndz", ctypes.c_double),
                ("rel_surf_tran", ctypes.c_double), 
                ("udreserved1", ctypes.c_double),
                ("udreserved2", ctypes.c_double),
                ("udreserved3", ctypes.c_double),
                ("udreserved4", ctypes.c_double),
                 ("string", 20*ctypes.c_byte)]
'''
typedef struct
{

   int type, numb;     /* the requested data type and number */
   int surf, wave;     /* the surface number and wavelength number */
   double wavelength, pwavelength;      /* the wavelength and primary wavelength */
   double n1, n2;      /* the index before and after */
   double cv, thic, sdia, k; /* the curvature, thickness, semi-diameter, and conic */
   double param[9];    /* the parameters 1-8 */
   double fdreserved1, fdreserved2, fdreserved3, fdreserved4; /* for future expansion */
   double xdata[201];  /* the extra data 1-200 */
   char glass[21];     /* the glass name on the surface */

}FIXED_DATA;
'''
class FIXED_DATA(ctypes.Structure):
    _fields_ = [
            ("type", ctypes.c_int),
            ("numb", ctypes.c_int),
            ("surf", ctypes.c_int),
            ("wave", ctypes.c_int),
            ("wavelength", ctypes.c_double), 
            ("n1", ctypes.c_double),
            ("n2", ctypes.c_double),
            ("cv", ctypes.c_double),
            ("thic", ctypes.c_double),
            ("sdia", ctypes.c_double),
            ("k", ctypes.c_double),
            ("param", 9*ctypes.c_double),
            ("fdreserved1", ctypes.c_double),
            ("xdata", 201*ctypes.c_double),
            ("glass", 21*ctypes.c_byte) 
        ]

# Load DLL into memory.


mydll = ctypes.CDLL ("./us_stand_gcc.so")

userdefinedsurface = mydll.UserDefinedSurface


u = USER_DATA()
f = FIXED_DATA()

f.k = -1.0
f.cv = 0.01
u.x = 1.0
u.y = 1.0
f.type = 3 # ask for sag

retval = userdefinedsurface(ctypes.byref(u), ctypes.byref(f))

print("\n", retval)
print(u.x, " ", u.y)
print(u.sag1)
print(u.sag2)

'''
# Set up prototype and parameters for the desired function call.
# HLLAPI

hllApiProto = ctypes.WINFUNCTYPE (
    ctypes.c_int,      # Return type.
    ctypes.c_void_p,   # Parameters 1 ...
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_void_p)   # ... thru 4.
hllApiParams = (1, "p1", 0), (1, "p2", 0), (1, "p3",0), (1, "p4",0),

# Actually map the call ("HLLAPI(...)") to a Python name.

hllApi = hllApiProto (("HLLAPI", hllDll), hllApiParams)

# This is how you can actually call the DLL function.
# Set up the variables and call the Python name with them.

p1 = ctypes.c_int (1)
p2 = ctypes.c_char_p (sessionVar)
p3 = ctypes.c_int (1)
p4 = ctypes.c_int (0)
hllApi (ctypes.byref (p1), p2, ctypes.byref (p3), ctypes.byref (p4))
'''


