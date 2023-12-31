use std::{ops::Deref, process::id};

pub fn dot_product(vec1: Vec<f64>,vec2: Vec<f64>) -> f64 {


    // make sure vector lengths are identical, otherwise, we 
    // have a problem (todo)

    // suppose the vectors are of equal length/dimension,
    // we can start multiplying

    let mut sum: f64 = 0.0;

    let vec_length = vec1.len();
    //dbg!(vec_length);

    // this is the simplest way

    for idx in 0..vec_length {

        let a = vec1[idx];
        let b = vec2[idx];

        let ab = a*b;

        sum = sum + ab;
    }


    // this is a more concise way using iter()

    let mut sum: f64 = 0.0;
    for (idx, a_ref) in vec1.iter().enumerate() {

        let a: f64 = *a_ref;
        let b = vec2[idx];
        let ab = a*b;

        sum = sum + ab;


    }

    // this is a more concise way using iter()


    //// how to add 2 to every element in vec1
    //let vec3: Vec<f64> = vec1.iter().map(|x| 2.0 + x).collect();
    //dbg!(vec3);


    let vec3: Vec<f64> = vec1.iter().enumerate().map(
        |(idx, a_ref)| {
            let a = *a_ref;
            let b = vec2[idx];
            let ab = a*b;

            ab

        }
    ).collect();

    let sum = vec3.iter().sum();

    return sum;


}

#[test]
fn dot_product_sandbox() {

    let vec1 = vec![0.1,0.2,0.5,124.3];
    let vec2 = vec![0.1,3.2,0.5,12.3];

    let scalar = dot_product(vec1, vec2);

    dbg!(scalar);

    assert_relative_eq!(1529.79, scalar);

}