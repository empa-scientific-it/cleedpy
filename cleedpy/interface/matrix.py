from ctypes import POINTER, Structure, c_double, c_int


class Matrix(Structure):
    _fields_ = [
        ("mag_no", c_int),
        ("blk_type", c_int),
        ("mat_type", c_int),
        ("num_type", c_int),
        ("rows", c_int),
        ("cols", c_int),
        ("rel", POINTER(c_double)),
        ("img", POINTER(c_double)),
    ]


MatPtr = POINTER(Matrix)
