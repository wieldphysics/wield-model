# -*- coding: utf-8 -*-
"""
"""
import numpy as np
from .utils import str_m, str_D

class ComplexBeamParam(object):
    """
    All distances should be in the same units as the wavelength_mgth.
    """
    pi = np.pi
    I = 1j

    # make this complex number dispatch smarter
    @staticmethod
    def complex(R, I = None):
        if I is not None:
            return R + 1j*I
        else:
            return R

    gouy_phasor = 1

    def __init__(
            self,
            value,
            wavelength_m,
            gouy_phasor = None,
    ):
        self.value = self.complex(value)
        if gouy_phasor is not None:
            self.gouy_phasor = gouy_phasor
        self.wavelength_m = wavelength_m

    @classmethod
    def from_Z_W0(
            cls,
            Z,
            W0,
            wavelength_m,
            gouy_phasor = None,
    ):
        ZR = cls.pi*W0**2/wavelength_m
        return cls(
            Z + cls.I*ZR,
            wavelength_m = wavelength_m,
            gouy_phasor = gouy_phasor,
        )

    @classmethod
    def from_W_R(
            cls,
            W,
            R,
            wavelength_m,
            gouy_phasor = None,
    ):
        if np.any(R != 0):
            iq = 1/R - 1j * wavelength_m / (np.pi * W**2)
            return cls(
                1/iq,
                wavelength_m = wavelength_m,
                gouy_phasor = gouy_phasor,
            )
        else:
            return cls(
                1j * (np.pi * W**2) / wavelength_m,
                wavelength_m = wavelength_m,
                gouy_phasor = gouy_phasor,
            )

    @classmethod
    def from_Z_ZR(
            cls,
            Z,
            ZR,
            wavelength_m = None,
            gouy_phasor = None,
    ):
        return cls(
            Z + cls.I*ZR,
            wavelength_m = wavelength_m,
            gouy_phasor = gouy_phasor,
        )

    @classmethod
    def from_dRad_Z(
            cls,
            theta,
            Z,
            wavelength_m,
            gouy_phasor = None,
    ):
        W0 = 2 * wavelength_m / (cls.pi * theta)
        return cls.from_Z_W0(
            Z,
            W0,
            wavelength_m = wavelength_m,
            gouy_phasor = gouy_phasor,
        )

    @property
    def W0(self):
        ZR = (self.value).imag
        return np.sqrt(self.wavelength_m * ZR / self.pi)

    @property
    def Z(self):
        return self.value.real

    @property
    def ZR(self):
        return (self.value).imag

    @property
    def cplg02(self):
        #the order matters right now because the Complex class is stupid
        return (self.value / self.value.conjugate()) / self.ZR * (self.gouy_phasor / self.gouy_phasor.conjugate())

    @property
    def W(self):
        return np.sqrt(-self.wavelength_m/((1/self.value).imag * self.pi))

    @property
    def R(self):
        return 1/((1/self.value).real)

    @property
    def R_inv(self):
        return ((1/self.value).real)

    @property
    def divergence_rad(self):
        return 2 * self.wavelength_m / (self.pi * self.W0)

    @property
    def k(self):
        return 2 * self.pi / self.wavelength_m

    @property
    def sensitivity_matrix(self):
        ZR = self.ZR
        Z = -self.Z
        return -self.k / (4*ZR) * np.matrix([[1, Z], [Z, Z**2 + ZR**2]])

    @property
    def sensitivity_matrix_sqrt(self):
        SM = self.sensitivity_matrix
        a = -SM[0, 0]
        b = -SM[0, 1]
        c = -SM[1, 0]
        d = -SM[1, 1]
        det = a*d - b*c
        trace = a + d
        s = det**.5
        t = (trace + 2 * s)**.5
        return np.matrix([[a + s, b], [c, d + s]]) / t

    def propagate_matrix(self, abcd_mat):
        a = abcd_mat[..., 0, 0]
        b = abcd_mat[..., 0, 1]
        c = abcd_mat[..., 1, 0]
        d = abcd_mat[..., 1, 1]
        determinant = a*d - b*c

        tanPhiA = self.wavelength_m * b
        tanPhiB = (a + b*self.R_inv)*self.pi*self.W**2

        #gouy_rad = np.arctan2(tanPhiA, tanPhiB)
        return self.__class__(
            (self.value * a + b)/(self.value * c + d),
            wavelength_m = self.wavelength_m * determinant,
            gouy_phasor = self.gouy_phasor * (tanPhiB + tanPhiA*1j),
        )

    def propagate_distance(self, L_m):
        a = 1
        b = L_m
        c = 0
        d = 1
        determinant = a*d - b*c

        tanPhiA = self.wavelength_m * b
        tanPhiB = (a + b*self.R_inv)*self.pi*self.W**2

        #gouy_rad = np.arctan2(tanPhiA, tanPhiB)
        return self.__class__(
            (self.value * a + b)/(self.value * c + d),
            wavelength_m = self.wavelength_m * determinant,
            gouy_phasor = self.gouy_phasor * (tanPhiB + tanPhiA*1j),
        )

    def __complex__(self):
        return self.value
        
    def overlap(self, q2):
        """
        Computes the projection from one beam parameter to another to give a measure of the
        overlap between the two beam parameters.
        
        This function was provided by Paul Fulda and Antonio Perreca, which came originally
        from Chris Mueller.
        
        Added on 20/4/2015

        From PyKat for testing
        """
        return abs(4*self.value.imag * q2.value.imag)/abs(self.value.conjugate()-q2.value)**2

    def overlap_HG(self, other):
        """
        This is the HG overlap
        """
        return ((
            (self.value * other.value.conjugate()) / (other.value.conjugate() - self.value)
            * (2 * self.wavelength_m / (other.W * self.W * self.pi))
        ) * -self.I)**0.5

    def overlap_HG(self, other):
        """
        This is the HG overlap
        """
        return ((
            (self.value * other.value.conjugate()) / (other.value.conjugate() - self.value)
            * (2 * self.wavelength_m / (other.W * self.W * self.pi))
        ) * -self.I)**0.5

    def overlap_HG_2mode(self, other):
        """
        This is the HG overlap
        """
        w1 = self.W
        w2 = other.W
        a = -self.I * (
            (self.value * other.value.conjugate()) / (other.value.conjugate() - self.value)
            * (2 * self.wavelength_m / (w2 * w1 * self.pi))
        )
        o = a**0.5
        #the abs isn't in the original formula, it is here to give the differential
        #phase between HG0 and HG2, rather than include the common phase of both from MM
        return o, (a*w1/w2 - 1) * abs(o)/2**0.5

    def overlap_LG(self, other):
        """
        This is the LG1 overlap
        """
        return -self.I * (
            (self.value * other.value.conjugate()) / (other.value.conjugate() - self.value)
            * (2 * self.wavelength_m / (other.W * self.W * self.pi))
        )

    def overlap_LG_2mode(self, other):
        w1 = self.W
        w2 = other.W
        a = -self.I * (
            (self.value * other.value.conjugate()) / (other.value.conjugate() - self.value)
            * (2 * self.wavelength_m / (w2 * w1 * self.pi))
        )
        o = a
        #the abs isn't in the original formula, it is here to give the differential
        #phase between LG0 and LG1, rather than include the common phase of both from MM
        return o, (a*(w1/w2) - 1) * abs(o)

    def diff_self_overlap_m(self, other):
        return -self.I * (
            (self.value * other.value.conjugate()) / (other.value.conjugate() - self.value)
            * (2 * self.wavelength_m / (other.W * self.W * self.pi))
        )

    def __str__(self):
        if self.R_inv != 0:
            R_str = u"{0}".format(str_m(self.R))
        else:
            R_str = u"1/0"
        return (u"Q(ZR={ZR}, Z={Z}, W0={W0}, W={W}, R={R}, λ={wavelength_m:.0f}nm)".format(
            ZR      = str_m(self.ZR, space = False),
            Z       = str_m(self.Z),
            W0      = str_m(self.W0, space = False),
            W       = str_m(self.W, space = False),
            wavelength_m = self.wavelength_m*1e9,
            R       = R_str,
        ))

    def str_kw(self):
        return dict(
            W  = str_m(self.W, space = False),
            diam  = str_m(2*self.W, space = False),
            Z  = str_m(self.Z, space = False),
            Ri = str_D(self.R_inv, space = False),
            W0 = str_m(self.W0, space = False),
            ZR = str_m(self.ZR, space = False),
            q_short = self.str_short(),
        )

    def str_short(self):
        return (u"Z={Z}, ZR={ZR}, W={W}, 1/R={Ri}".format(
            W  = str_m(self.W, space = False),
            Z  = str_m(self.Z, space = False),
            Ri = str_D(self.R_inv, space = False),
            W0 = str_m(self.W0, space = False),
            ZR = str_m(self.ZR, space = False),
        ))

    def string(self):
        return (u"W={W}, Z={Z}, 1/R={Ri}, W0={W0}, ZR={ZR}".format(
            W  = str_m(self.W, space = False),
            Z  = str_m(self.Z),
            Ri = str_D(self.R_inv),
            W0 = str_m(self.W0, space = False),
            ZR = str_m(self.ZR, space = False),
        ))

    def reversed(self):
        return self.__class__(
            value = -self.value.conjugate(),
            wavelength_m = self.wavelength_m,
            gouy_phasor = self.gouy_phasor
        )


def beam_shape_1D(x, x_cbp, trans, tilt):
    """
    .. todo:: check the that the tilt and shift sign conventions make sense
    """
    cbp = complex(x_cbp)
    k = 2*np.pi / x_cbp.wavelength_m
    return (2/np.pi)**.25 * (1/x_cbp.W)**.5 * np.exp(-1j * k * (x + trans) ** 2 / (2*cbp) + 2*np.pi * 1j * x/x_cbp.wavelength_m * np.sin(tilt))


def beam_transverse(x, y, x_cbp, y_cbp = None, x_trans = 0, y_trans = 0, x_tilt = 0, y_tilt = 0):
    if y_cbp is None:
        y_cbp = x_cbp

    return np.outer(
        beam_shape_1D(y, y_cbp, trans = y_trans, tilt = y_tilt),
        beam_shape_1D(x, x_cbp, trans = x_trans, tilt = x_tilt),
    )


