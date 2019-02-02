#![allow(unused_imports)]
extern crate num_complex as nc;
use self::nc::Complex;

extern crate num_traits as nt;
use self::nt::{cast::FromPrimitive, Num};

use std::f64::consts;
use std::ops::Neg;

pub trait Inte where Self: FromPrimitive + Num + Copy + Neg<Output=Self> {
    fn legendre(z: Self, i: i8) -> Self {
        let j = Self::from_i8(i).unwrap();
        let o = Self::from_i8(1).unwrap();
        let t = Self::from_i8(2).unwrap();
        match i {
            0 => return o,
            1 => return z,
            _ => return ((t*j-o)*z*Inte::legendre(z, i-1)-(j-o)*Inte::legendre(z, i-2))/j
        }
    }

    fn dlegendre(z: Self, i: i8) -> Self {
        let j = Self::from_i8(i).unwrap();
        let f = j / (z * z - Self::from_i8(1).unwrap());
        (z * Inte::legendre(z, i) - Inte::legendre(z, i-1)) * f
    }

    /// Integrate a function with real or complex bounds.
    /// 
    /// # Example
    /// ```rust
    /// Inte::numint(1.5, 0.5, &my_func_name)
    /// ```
    fn numint(u: Self, l: Self, f: &Fn(Self) -> Self) -> Self {
        let o = Self::from_i8(1).unwrap();
        let t = Self::from_i8(2).unwrap();

        let a: [f64; 10] = [
            0.1488743389816312_f64,
            -0.1488743389816312_f64,
            0.43339539412924719_f64,
            -0.43339539412924719_f64,
            0.6794095682990244_f64,
            -0.6794095682990244_f64,
            0.8650633666889845_f64,
            -0.8650633666889845_f64,
            0.9739065285171717_f64,
            -0.9739065285171717_f64,
        ];

        let mut r = Self::from_i8(0).unwrap();

        for i in 0..10 {
            let b = Self::from_f64(a[i]).unwrap();
            let w = t / ((o - Inte::quad(b)) * Inte::quad(Inte::dlegendre(b, 10)));
            r = r + w * f((u - l) / t * b + (l + u) / t)
        }
        (u - l) / t * r
    }

    fn quad(z: Self) -> Self {
        z * z
    }
}

macro_rules! impl_Inte {
    ($($t:ty),+) => {$(
        impl Inte for $t {}
    )+};
}

impl_Inte!(f32, f64, Complex<f32>, Complex<f64>);
