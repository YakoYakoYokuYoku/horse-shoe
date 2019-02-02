extern crate num_complex as nc;
use nc::Complex;

extern crate num_traits as nt;
use nt::cast::{FromPrimitive};

fn bicoe<U: Copy + Into<i128>>(a: U, b: U) -> i128 {
    let n = a.into();
    let k = b.into();
    let mut re: i128 = 1;

    for i in 0..k {
        re *= (n - 1) / (i + 1)
    }

    re
}

macro_rules! impl_Der {
    (x => $t:ty) => (
        fn numder(z: $t, fun: &Fn($t) -> $t, i: i8) -> $t {
            let e = <$t>::from_f64(1e-10).unwrap();
            let mut n = <$t>::from_i8(-1).unwrap();
            let mut r = <$t>::from_i8(0).unwrap();
            for x in 0..i+1 {
                let y = <$t>::from_i8(x).unwrap();
                let f = <$t>::from_i128(bicoe(i, x)).unwrap();
                r = fun(z+e*y) * f * n + r;
                n = -n;
            }

            return r
        }
    );
}

impl_Der!(x => Complex<f64>);

pub fn digamma(mut i: Complex<f64>) -> Complex<f64> {
    let ze = 0.0;
    let on = 1.0;
    i = i - on;
    let mut ret = Complex::new(ze, ze);

    for x in 0..1000 {
        let y = x as f64;
        ret = ret + on / y - on / (y + i);
    }
    
    ret = ret - 0.57721566490153286060651209008240243_f64;

    ret
}

fn main() {
    for i in 0..10 {
        let a = numder(Complex { re: 0.5_f64, im: 0.0_f64 }, &digamma, i);

        println!("{:?}", a);
    }
}
