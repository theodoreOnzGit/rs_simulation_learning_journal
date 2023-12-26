use std::f64::consts::PI;

use uom::si::{f64::*, angle::{radian, degree}};

pub fn find_cubic_roots(b: f64,
c: f64,
d: f64) -> Vec<f64>{

    // r = (p/3)^3 + (q/2)^2
    let r: f64;
    let p: f64;
    let q: f64;

    // p = (3c - b^2)/3.0
    p = (3.0 * c - b.powf(2.0))/3.0;

    // q = (27d - 9bc + 2b^3)/27

    q = (27.0 * d - 9.0*b*c + 2_f64*b.powf(3.0))/27.0;

    r = (p/3.0).powf(3.0) + (q/2.0).powf(2.0);

    // what we do next depends on r 

    let y_0: f64;
    // note that y_1 and y_2 may or may not be imaginary
    let y_1: Option<f64>;
    let y_2: Option<f64>;

    
    if r < 0.0 {

        // expect 3 real roots

        fn get_y_k(k: usize, p: f64,
        phi_degrees: f64,
        q:f64) -> f64 {
            let y_k;

            if q > 0.0 {
                
                y_k = -2.0 * (-p/3.0).sqrt()*
                ((phi_degrees/3.0 + 120.0 * k as f64) * PI/180.0).cos();


            } else if q < 0.0 {

                y_k = 2.0 * (-p/3.0).sqrt()*
                ((phi_degrees/3.0 + 120.0 * k as f64) * PI/180.0).cos();

            } else {
                // q = 0.0 case
                // q = 0.0, means that phi is Acos(0)
                // this implies phi is 90 degrees or 
                // pi/2

                y_k = 2.0 * (-p/3.0).sqrt()*
                ((phi_degrees/3.0 + 120.0 * k as f64)*PI/180.0).cos();


            }

            return y_k;
        }



        let phi: f64;

        let q_sq_by_4: f64 = q.powf(2.0)*0.25;
        let minus_p_cube_by_27: f64 = -p.powf(3.0)/27.0;

        // note that phi is in radians
        phi = (q_sq_by_4/minus_p_cube_by_27).acos();

        let phi_dimensioned = Angle::new::<radian>(phi);
        let phi_degrees: f64 = phi_dimensioned.get::<degree>();

        y_0 = get_y_k(0, p, phi_degrees, q);
        y_1 = Some(get_y_k(1, p, phi_degrees, q));
        y_2 = Some(get_y_k(2, p, phi_degrees, q));


    } else if r > 0.0 {


        // r > 0.0
        // expect 1 real root, 2 complex roots
        // which are conjugates of each other

        let A: f64;
        let B: f64;

        let minus_q_by_2: f64 = -q/2.0;

        A = (minus_q_by_2 + r.sqrt()).cbrt();
        B = (minus_q_by_2 - r.sqrt()).cbrt();

        // real root
        y_0 = A + B;

        // other roots, y_1 and y_2 are imaginary
        // in this case, we are not interested in imaginary 
        // roots, so y_1 and y_2 are undefined
        // the trouble with y_1 and y_2 is that sometimes,
        // y_1 and y_2 are undefined (imaginary)
        // we use the None variant of the Option Enum 
        // to say that y_1 and y_2 do not exist
        y_1 = None;
        y_2 = None;


        


    } else {

        // expect 3 real roots, two of which are 
        // equal 

        let A: f64;
        let B: f64;

        let minus_q_by_2: f64 = -q/2.0;

        A = (minus_q_by_2).cbrt();
        B = (minus_q_by_2).cbrt();

        // real root
        y_0 = A + B;

        // two other equal real roots
        // we use the Some variant for the 
        // Option Enum to say that, hey y_1 and y_2 exist
        y_1 = Some(-0.5*(y_0));
        y_2 = Some(-0.5*(y_0));

    }

    fn get_x_k(y_k: f64, b: f64) -> f64 {
        let x_k = y_k - b/3.0;

        return x_k;
    }
    
    let x_0 = get_x_k(y_0, b);

    let mut solution_vec = vec![x_0];
    
    match y_1 {
        Some(y_1_value) => {
            solution_vec.push(get_x_k(y_1_value, b));
        },
        None => (),
    }

    match y_2 {
        Some(y_2_value) => {
            solution_vec.push(get_x_k(y_2_value, b));
        },
        None => (),
    }


    return solution_vec;

}

#[test]
fn test_three_real_roots(){
    
    // x^3 - 7x + 7 = 0
    let b = 0.0;
    let c = -7.0;
    let d = 7.0;
    let solution = find_cubic_roots(b, c, d);

    panic!("{:?}",solution);

}