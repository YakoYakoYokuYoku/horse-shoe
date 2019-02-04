#![allow(unused_imports)]
extern crate num_complex as nc;
use self::nc::Complex;

extern crate num_traits as nt;
use self::nt::{cast::FromPrimitive, Num, NumCast, Float, float::{CommonFloat, FloatCore}};

use std::{f64, ops::Neg, cmp::PartialOrd};

pub trait Funcs<T>
where
    Self: FromPrimitive + Num + Copy + CommonFloat<Typo=T>,
    T: PartialOrd,
{
    fn legendre(z: Self, i: i8) -> Self {
        let j = Self::from_i8(i).unwrap();
        let o = Self::from_i8(1).unwrap();
        let t = Self::from_i8(2).unwrap();
        match i {
            0 => return o,
            1 => return z,
            _ => {
                return ((t * j - o) * z * Funcs::legendre(z, i - 1)
                    - (j - o) * Funcs::legendre(z, i - 2))
                    / j;
            }
        }
    }

    fn dlegendre(z: Self, i: i8) -> Self {
        let j = Self::from_i8(i).unwrap();
        let f = j / (z * z - Self::from_i8(1).unwrap());
        (z * Funcs::legendre(z, i) - Funcs::legendre(z, i - 1)) * f
    }

    fn gamma(mut self) -> Self {
        let ha = Self::from_f32(0.5).unwrap();
        let on = Self::from_i8(1).unwrap();
        let fh = Self::from_f32(7.5).unwrap();
        let lan = Self::from_f64(2.50662827463100050241576528481104525_f64).unwrap();
        let arq = Self::from_f64(f64::consts::PI).unwrap();

        let arr: [f64; 9] = [
            0.99999999999980993,
            676.5203681218851,
            -1259.1392167224028,
            771.32342877765313,
            -176.61502916214059,
            12.507343278686905,
            -0.13857109526572012,
            9.9843695780195716e-6,
            1.5056327351493116e-7,
        ];

        if self.real() >= ha.real() {
            self = self - on;
            let mut x = Self::from_f64(arr[0]).unwrap();
            for q in 1..9 {
                let d = Self::from_i8(q).unwrap();
                let v = q as usize;
                let lec = Self::from_f64(arr[v]).unwrap();
                let xop = on / (self + d);
                x = x + xop * lec;
            }
            let fa = self + fh;
            let fb = fa.pown(self + ha);
            let ep = (-fa).exp();
            let rt = x * fb * ep * lan;

            rt
        } else {
            let efm = (self * arq).sin();
            let gam = Funcs::gamma(-self + on);
            let inv = on / (efm * gam);
            let ret = inv * arq;

            ret
        }
    }
}

impl<T> Funcs<T> for T 
where
    T: PartialOrd + CommonFloat<Typo = T> + FromPrimitive
{}

impl<T> Funcs<T> for Complex<T>
where
    T: FloatCore + Float + FromPrimitive
{}
