use uom::si::{f64::*, angle::{radian, degree}};

pub fn find_cubic_roots(b: f64,
c: f64,
d: f64){

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
    
    if r < 0.0 {

        // expect 3 real roots

        fn get_y_k(k: usize, p: f64,
        phi: f64,
        q:f64) -> f64 {
            let y_k;

            if q > 0.0 {
                
                y_k = -2.0 * (-p/3.0).sqrt()*
                (phi/3.0 + 120.0 * k as f64);


            } else if q < 0.0 {

                y_k = 2.0 * (-p/3.0).sqrt()*
                (phi/3.0 + 120.0 * k as f64);

            } else {
                // q = 0.0 case

                todo!("q=0.0 case for cubic root finding not done yet")


            }

            return y_k;
        }

        fn get_x_k(y_k: f64, b: f64) -> f64 {
            let x_k = y_k - b/3.0;

            return x_k;
        }


        let phi: f64;

        let q_sq_by_4: f64 = q.powf(2.0)*0.25;
        let minus_p_cube_by_27: f64 = -p.powf(3.0)/27.0;

        // note that phi is in radians
        phi = (q_sq_by_4/minus_p_cube_by_27).acos();

        let phi_dimensioned = Angle::new::<radian>(phi);
        let phi_degrees: f64 = phi_dimensioned.get::<degree>();


    } else if r == 0.0 {

        // expect 3 real roots, two of which are 
        // equal 

        let y_k: f64;
        let x_k: f64;

        let phi: f64;

        let q_sq_by_4: f64 = q.powf(2.0)*0.25;
        let minus_p_cube_by_27: f64 = -p.powf(3.0)/27.0;

        // note that phi is in radians
        phi = (q_sq_by_4/minus_p_cube_by_27).acos();

        let phi_dimensioned = Angle::new::<radian>(phi);
        let phi_degrees: f64 = phi_dimensioned.get::<degree>();

        

        if q > 0.0 {


        } else if q < 0.0 {


        } else {
            // q = 0.0 case


        }


    } else {

        // r > 0.0
        // expect 1 real root, 2 complex roots
        // which are conjugates of each other


    }


}