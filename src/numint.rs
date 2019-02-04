extern crate num_complex as nc;
use self::nc::Complex;

extern crate num_traits as nt;
use self::nt::{cast::FromPrimitive, Num};

use std::{ops::Neg, cmp::PartialOrd};

use functions::Funcs;

pub trait Intg<I>
where
    Self: FromPrimitive + Num + Copy + Neg<Output = Self> + Funcs<I>,
    I: PartialOrd,
{
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
            let w = t / ((o - Intg::quad(b)) * Intg::quad(Funcs::dlegendre(b, 10)));
            r = r + w * f((u - l) / t * b + (l + u) / t)
        }
        (u - l) / t * r
    }

    fn quad(z: Self) -> Self {
        z * z
    }
}

macro_rules! impl_Intg {
    ($($t:ty),+) => {$(
        impl<I> Intg<I> for $t
        where
            I: Neg<Output=I> + Copy + Num + FromPrimitive + Funcs<I> + PartialOrd,
            Complex<I>: Funcs<I>,
        {}
    )+};
}

impl_Intg!(I, Complex<I>);
