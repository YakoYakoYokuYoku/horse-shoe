extern crate num_complex as nc;
use self::nc::Complex;

extern crate num_traits as nt;
use self::nt::{cast::FromPrimitive, Num};

pub trait Diff where Self: FromPrimitive + Num + Copy {
    /// Take the derivative of a function at a certain real or complex number.
    /// 
    /// # Example
    /// ```rust
    /// Diff::numder(0.5, &my_func_name)
    /// ```
    fn numder(z: Self, f: &Fn(Self) -> Self) -> Self {
        let h = Self::from_f64(1e-10).unwrap();
        (f(z + h) - f(z)) / h
    }
}

macro_rules! impl_Diff {
    ($($t:ty),+) => {$(
        impl Diff for $t {}
    )+};
}

impl_Diff!(f32, f64, Complex<f32>, Complex<f64>);
