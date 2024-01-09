mod newtons_law_of_cooling;
use uom::si::f64::*;
use uom::si::molar_concentration::{mole_per_cubic_decimeter, mole_per_cubic_meter};
use uom::si::molar_volume::{cubic_meter_per_mole, cubic_decimeter_per_mole};
use uom::si::pressure::{atmosphere, bar};
use uom::si::ratio::ratio;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::pressure::pascal;

use crate::newtons_law_of_cooling::implicit_euler_method;

mod ideal_gas;
use crate::ideal_gas::*;

mod pr_eos;
use pr_eos::*;

mod cubic_roots;
mod quadratic_roots;

mod vector_multiplication;

#[macro_use]
extern crate approx;


fn main() {

    println!("Hello, world! I saved my file");
    let x = add_expt_1();
    dbg!(x);

    implicit_euler_method();

    let specific_vol = MolarVolume::new::<cubic_meter_per_mole>(1.0);
    let pressure =  Pressure::new::<atmosphere>(1.0);

    let temp_ideal_gas = get_temperature_ideal_gas(
        pressure, specific_vol);
    
    println!("{:?}",temp_ideal_gas);
    dbg!(temp_ideal_gas);


}

#[test]
fn test_co2_pressure_peng_robinson_subcritial(){

    let co2_pressure = get_pressure_peng_robinson_gas(
        ThermodynamicTemperature::new::<kelvin>(250.0),
        MolarVolume::new::<cubic_decimeter_per_mole>(24.5), 
        Ratio::new::<ratio>(0.239), 
        Pressure::new::<bar>(73.8), 
        ThermodynamicTemperature::new::<kelvin>(304.0));

    let reference_co2_pressure_wolfram = Pressure::new::<pascal>(84183.9);

    assert_relative_eq!(
        co2_pressure.get::<pascal>(),
        reference_co2_pressure_wolfram.get::<pascal>(),
    max_relative = 0.001);

}

fn add_expt_1() -> f32 {
    
    let a = 3;
    let b = 6;

    let c = a + b;

    println!("{}",c);

    //for i in 0..3{
    //    dbg!(&i);
    //}

    1 as f32 +1.5
}
