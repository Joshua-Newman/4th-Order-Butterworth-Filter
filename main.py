import control
import numpy as np

num = [0, 0, 0, 0, 9.3502 * 10 ** 15]
den = [1, 17182, 147610562, 742878907342, 1869341405846161]

sys = control.TransferFunction(num, den)

a, b, c, d = control.ssdata(sys)
a = np.flip(a, (0, 1))
b = np.flip(b, (0, 1))
c = np.flip(c, (0, 1))
d = np.flip(d, (0, 1))
print(a)
print(b)
print(c)
print(d)


#Z^4 + 17182 Z^3 + 147610562 Z^2 + 742878907342 Z + 1869341405846161