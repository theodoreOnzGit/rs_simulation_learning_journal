use std::{ops::Deref, process::id};

pub fn dot_product(vec1: Vec<f64>,vec2: Vec<f64>) -> Option<f64> {


    // make sure vector lengths are identical, otherwise, we 
    // have a problem (todo)
    
    if vec1.len() != vec2.len() {
        return None;
    }

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

    return Some(sum);


}

#[test]
fn dot_product_sandbox() {

    let vec1 = vec![0.1,0.2,0.5,124.3];
    let vec2 = vec![0.1,3.2,0.5,12.3];

    let scalar = dot_product(vec1, vec2);

    let value = match scalar {
        Some(value) => value,
        None => panic!(),
    };

    dbg!(value);

    assert_relative_eq!(1529.79, value);

}

#[test]
fn dot_product_error_handling() {

    let vec1 = vec![0.1,0.2,0.5,124.3];
    let vec2 = vec![0.1,3.2,0.5];

    let scalar = dot_product(vec1, vec2);

    let value = scalar.unwrap();

    dbg!(value);

    assert_relative_eq!(1529.79, value);

}




#[test]
fn enum_demo(){
    #[derive(Debug)]
    pub enum EngineeringGrades {
        A(f64),
        B(f64),
        C(f64),
        D(f64),
        E(f64),
        F(f64)
    }
    // let's say on a scale of 0 to 100
    // 75% and above is an A

    fn convert_score_to_grade(score: f64) -> EngineeringGrades {
        if score > 75.0 {
            return EngineeringGrades::A(score);
        } else if score > 65.0 {
            return EngineeringGrades::B(score);
        } else if score > 60.0 {
            return EngineeringGrades::C(score);
        } else if score > 55.0 {
            return EngineeringGrades::D(score);
        } else if score > 50.0 {
            return EngineeringGrades::E(score);
        } else {
            return EngineeringGrades::F(score);
        }
    }
    let my_engineering_grade = EngineeringGrades::B(67.0);

    match my_engineering_grade {
        EngineeringGrades::A(_score) => {
            println!("well done, you are on the Dean's list");
        },
        EngineeringGrades::F(_score) => {
            println!("you need to repeat the course");
        },
        EngineeringGrades::B(_score) => {
            ()
        },
        EngineeringGrades::C(_score) => {
            ()
        },
        EngineeringGrades::D(_score) => {
            ()
        },
        EngineeringGrades::E(_score) => {
            ()
        }

    }
    dbg!(&my_engineering_grade);

    let my_score: f64 = 84.0;

    let my_engineering_grade = convert_score_to_grade(my_score);

    dbg!(&my_engineering_grade);

    // print special message based on your enum
    // if you get A, then print well done!

    match my_engineering_grade {
        EngineeringGrades::A(_) => {
            println!("well done, you are on the Dean's list");
        },
        EngineeringGrades::F(_) => {
            println!("you need to repeat the course");
        },
        // for any score not A and F, do nothing
        _ => (),

    }

}


