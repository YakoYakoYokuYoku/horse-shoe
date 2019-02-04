use numint as ni;

extern crate num_complex as nc;
use self::nc::Complex;

extern crate num_traits as nt;
use self::nt::{cast::FromPrimitive, Num, float::CommonFloat};

use functions::Funcs;

use std::cmp::PartialOrd;

pub trait Diff<D>
where
    Self: FromPrimitive + Num + Copy + ni::Intg<D> + CommonFloat<Typo=D>,
    D: PartialOrd,
{
    fn numder(z: Self, f: &Fn(Self) -> Self) -> Self {
        let h = Self::from_f64(1e-10).unwrap();
        (f(z + h) - f(z)) / h
    }
}

impl<D> Diff<D> for D
where
    D: PartialOrd + Funcs<D>,
    Complex<D>: Funcs<D>
{}

impl<D> Diff<D> for Complex<D> 
where
    D: PartialOrd + Funcs<D>,
    Complex<D>: Funcs<D>,
{}
