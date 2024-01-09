
use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::molar_heat_capacity::joule_per_kelvin_mole;

use crate::cubic_roots::find_cubic_roots;
use crate::quadratic_roots::find_real_quadratic_roots;

// if you want a new custom unit 
// https://github.com/iliekturtles/uom/issues/244

///
pub fn get_temperature_peng_robinson_gas(
    pressure: Pressure,
    specific_vol: MolarVolume, 
    acentricity_factor: Ratio, 
    critical_temperature: ThermodynamicTemperature, 
    critical_pressure: Pressure) -> 
Result<ThermodynamicTemperature, CubicEOSError>{

        let b = get_b(critical_pressure, critical_temperature);
        let kappa = get_kappa(acentricity_factor);
        // R
        let molar_gas_constant = 
        MolarHeatCapacity::new::<joule_per_kelvin_mole>(8.314);

        let psi_1: Ratio = 
            pressure * (specific_vol - b) / 
            (molar_gas_constant * critical_temperature);

        // for psi_2

        let psi_2: Ratio = {

            
            let left_term = (specific_vol - b) / 
            (molar_gas_constant * critical_temperature);

            let alpha_coeff = 0.45724 * molar_gas_constant 
                * molar_gas_constant 
                * critical_temperature 
                * critical_temperature 
                / critical_pressure;

            // denominator: (V(V+b) + b(V-b))
            let denominator = specific_vol * (specific_vol + b)
                + b * (specific_vol - b);

            let psi_2 = left_term * alpha_coeff / 
                denominator;

            psi_2


        };
        let kappa_sq: Ratio = kappa * kappa;
        let ratio_a: Ratio = psi_2 * kappa_sq - 
            Ratio::new::<ratio>(1.0);

        let ratio_b: Ratio = -psi_2 * (2.0* kappa_sq - 
            2.0 * kappa);

        let ratio_c: Ratio = psi_1 + psi_2 *(
            Ratio::new::<ratio>(1.0)
            + kappa_sq
            + 2.0 * kappa);

        //let sqrt_reduced_temp_array: [Ratio;2] 
        //    = find_real_quadratic_roots(
        //        ratio_a,
        //        ratio_b,
        //        ratio_c).unwrap();

        //let reduced_temp_arr_result = 
        //    find_real_quadratic_roots(
        //        ratio_a, 
        //        ratio_b, 
        //        ratio_c);
        //
        //let sqrt_reduced_temp_array: [Ratio;2] = 
        //    match reduced_temp_arr_result {
        //        Ok(array) => array,
        //        Err(error) => {
        //            return Err(error);
        //        }
        //    };

        let sqrt_reduced_temp_array: [Ratio;2] 
            = find_real_quadratic_roots(
                ratio_a,
                ratio_b,
                ratio_c)?;

        // todo: which root is the better one to use?
        let reduced_temperature = 
            sqrt_reduced_temp_array[1]
            * sqrt_reduced_temp_array[1];

        let gas_temperature: ThermodynamicTemperature = 
            reduced_temperature.get::<ratio>() * 
            critical_temperature;

        // to do for next time

        return Ok(gas_temperature);
}
#[test]
fn superheated_steam_temperature_test_pr_eos(){
 
    use uom::si::thermodynamic_temperature::{kelvin,degree_celsius};
    use uom::si::pressure::megapascal;
    use uom::si::molar_volume::cubic_meter_per_mole;
    use uom::si::specific_volume::cubic_meter_per_kilogram;
    use uom::si::molar_mass::gram_per_mole;
    // expected behaviour
    let water_crit_temperature = 
        ThermodynamicTemperature::new::<kelvin>(647.1);
    let water_crit_pressure = 
        Pressure::new::<megapascal>(22.06);
    
    let expected_steam_temp = ThermodynamicTemperature::new::<degree_celsius>(600.0);
    let steam_pressure = Pressure::new::<megapascal>(0.40);
    let steam_specific_vol = SpecificVolume::new::
        <cubic_meter_per_kilogram>(1.00558);

    let water_mol_wt = MolarMass::new::<gram_per_mole>(18.0);

    let steam_molar_vol = steam_specific_vol * water_mol_wt;

    //let steam_molar_vol = steam_specific_vol


    // set the input parameters for the function 

    let pressure = steam_pressure;
    let critical_pressure = water_crit_pressure;
    let critical_temperature = water_crit_temperature;
    let acentricity_factor = Ratio::new::<ratio>(0.344);

    // find temperature using pr eos

    let steam_result_temperature: ThermodynamicTemperature =
        get_temperature_peng_robinson_gas(
            pressure, 
            steam_molar_vol, 
            acentricity_factor, 
            critical_temperature, 
            critical_pressure).unwrap();

    // assert approx eq to within 0.1% 

    assert_relative_eq!(steam_result_temperature.get::<kelvin>(),
                        expected_steam_temp.get::<kelvin>(),
                        max_relative=0.001);

}

pub fn get_pressure_peng_robinson_gas(
    temperature: ThermodynamicTemperature,
    molar_volume: MolarVolume,
    acentricity_factor: Ratio,
    critical_pressure: Pressure,
    critical_temperature: ThermodynamicTemperature) -> Pressure{
        
        let b = get_b(critical_pressure, critical_temperature);
        let kappa = get_kappa(acentricity_factor);

        // R
        let molar_gas_constant = 
        MolarHeatCapacity::new::<joule_per_kelvin_mole>(8.314);

        // (RT)/(V-b)
        // V here is the molar volume
        let rt_over_v_minus_b: Pressure 
        = molar_gas_constant* temperature/(molar_volume-b);

        // now to calculate a(T)/(V(V+b) + b(V-b))
        // denominator: (V(V+b) + b(V-b))
        let denominator = molar_volume * (molar_volume + b)
        + b * (molar_volume - b);

        // numerator: a(T) or a_t
        // a_t = 0.45724 r^2*T_c^2/P_c * alpha(t)
        // alpha(t) = [1 + kappa * (1- sqrt(T_r))]^2
        // T_r  is reduced temperature
        let reduced_temperature: Ratio 
        = temperature/critical_temperature;

        let sqrt_alpha = Ratio::new::<ratio>(1.0)
        + kappa * (Ratio::new::<ratio>(1.0)
        - reduced_temperature.sqrt());

        let alpha = sqrt_alpha * sqrt_alpha;

        let numerator = 0.45724 * molar_gas_constant 
        * molar_gas_constant 
        * critical_temperature 
        * critical_temperature 
        / critical_pressure 
        * alpha;

        let acentricity_term: Pressure = 
        numerator/denominator;


        let pressure = rt_over_v_minus_b - acentricity_term;
        return pressure;
}

pub fn get_molar_volume_peng_robinson_eos_gas_and_vapour_like_roots(
    pressure: Pressure,
    temperature: ThermodynamicTemperature,
    acentricity_factor: Ratio,
    critical_pressure: Pressure,
    critical_temperature: ThermodynamicTemperature) -> MolarVolume{
        
        let b = get_b(critical_pressure, critical_temperature);

        // R
        let molar_gas_constant = 
        MolarHeatCapacity::new::<joule_per_kelvin_mole>(8.314);

        // B
        let capital_b: Ratio = b * pressure / molar_gas_constant 
        / temperature;

        // A
        let capital_a: Ratio;

        let kappa = get_kappa(acentricity_factor);
        let reduced_temperature: Ratio 
        = temperature/critical_temperature;
        let sqrt_alpha = Ratio::new::<ratio>(1.0)
        + kappa * (Ratio::new::<ratio>(1.0)
        - reduced_temperature.sqrt());
        let alpha: Ratio = sqrt_alpha * sqrt_alpha;

        let a = 0.45724 * molar_gas_constant 
        * molar_gas_constant 
        * critical_temperature 
        * critical_temperature 
        / critical_pressure 
        * alpha;

        capital_a = a * pressure / molar_gas_constant 
        / molar_gas_constant
        / temperature
        / temperature;

        // Z^3 - (1-B) Z^2 + (A - 3B^2 - 2B) Z 
        // - (AB-B^2 -B^3)
        let one_minus_capital_b: Ratio 
        = Ratio::new::<ratio>(1.0) - capital_b;
        let capital_a_minus_3_cap_b_sq_minus_2_cap_b: Ratio 
        = capital_a - 3.0 * capital_b * capital_b 
        - 2.0 * capital_b;

        let cap_ab_minus_cap_b_sq_minus_b_cubed: Ratio 
        = capital_a * capital_b - 
        capital_b * capital_b 
        - capital_b * capital_b * capital_b;

        // find roots
        let cubic_root_b: f64 = -one_minus_capital_b.get::<ratio>();
        let cubic_root_c: f64 = capital_a_minus_3_cap_b_sq_minus_2_cap_b.
            get::<ratio>();
        let cubic_root_d: f64 = -cap_ab_minus_cap_b_sq_minus_b_cubed.
            get::<ratio>();

        let z_root_values: Vec<f64> = 
            find_cubic_roots(cubic_root_b,cubic_root_c,cubic_root_d);

        // case 1: one real root
        if z_root_values.len() == 1 {
            let z_value: f64 = *z_root_values.first().unwrap();
            
            #[allow(non_snake_case)]
            let Z = Ratio::new::<ratio>(z_value);

            // molar volume 
            let v_tilde: MolarVolume = Z * molar_gas_constant * 
                temperature / pressure;

            return v_tilde;
        }

        // case 2: three real roots
        if z_root_values.len() == 3 {

            // order the z roots first if there are three real roots
            
            let mut soln_vec_clone = z_root_values.clone();
            soln_vec_clone.sort_by(|a, b| a.partial_cmp(b).unwrap());

            // max root of Z is the maximum compressibility factor
            // this corresponds to vapour like roots
            let max_root = *soln_vec_clone.last().unwrap();
            // max root of Z is the maximum compressibility factor
            // this corresponds to liquid like roots
            // let min_root = *soln_vec_clone.first().unwrap();
            //
            //
            #[allow(non_snake_case)]
            let Z = Ratio::new::<ratio>(max_root);

            // molar volume 
            let v_tilde: MolarVolume = Z * molar_gas_constant * 
                temperature / pressure;

            return v_tilde;
        }
        


        panic!("unable to find vapour like volume roots");
}

#[test]
fn crit_point_water_test_pr_eos(){
 
    use uom::si::thermodynamic_temperature::kelvin;
    use uom::si::pressure::megapascal;
    use uom::si::molar_volume::cubic_meter_per_mole;
    // expected behaviour
    let water_crit_temperature = 
        ThermodynamicTemperature::new::<kelvin>(647.1);
    let water_crit_pressure = 
        Pressure::new::<megapascal>(22.06);
    let water_crit_vol = 
        MolarVolume::new::<cubic_meter_per_mole>(0.0560/1000.0);

    // set the input parameters for the function 

    let pressure = water_crit_pressure;
    let critical_pressure = water_crit_pressure;
    let temperature = water_crit_temperature;
    let critical_temperature = water_crit_temperature;
    let acentricity_factor = Ratio::new::<ratio>(0.344);

    // find volume using pr eos
    let water_vap_vol_pr_eos: MolarVolume = 
        get_molar_volume_peng_robinson_eos_gas_and_vapour_like_roots(
            pressure, 
            temperature, 
            acentricity_factor, 
            critical_pressure, 
            critical_temperature);

    // assert approx eq to within 30%

    assert_relative_eq!(water_vap_vol_pr_eos.get::<cubic_meter_per_mole>(),
                        water_crit_vol.get::<cubic_meter_per_mole>(),
                        max_relative=0.30);

}

#[test]
fn superheated_steam_test_pr_eos(){
 
    use uom::si::thermodynamic_temperature::{kelvin,degree_celsius};
    use uom::si::pressure::megapascal;
    use uom::si::molar_volume::cubic_meter_per_mole;
    use uom::si::specific_volume::cubic_meter_per_kilogram;
    use uom::si::molar_mass::gram_per_mole;
    // expected behaviour
    let water_crit_temperature = 
        ThermodynamicTemperature::new::<kelvin>(647.1);
    let water_crit_pressure = 
        Pressure::new::<megapascal>(22.06);
    
    let steam_temp = ThermodynamicTemperature::new::<degree_celsius>(600.0);
    let steam_pressure = Pressure::new::<megapascal>(0.40);
    let steam_specific_vol = SpecificVolume::new::
        <cubic_meter_per_kilogram>(1.00558);

    let water_mol_wt = MolarMass::new::<gram_per_mole>(18.0);

    let expected_steam_molar_vol = steam_specific_vol * water_mol_wt;

    //let steam_molar_vol = steam_specific_vol


    // set the input parameters for the function 

    let pressure = steam_pressure;
    let critical_pressure = water_crit_pressure;
    let temperature = steam_temp;
    let critical_temperature = water_crit_temperature;
    let acentricity_factor = Ratio::new::<ratio>(0.344);

    // find volume using pr eos
    let pr_eos_molar_vol: MolarVolume = 
        get_molar_volume_peng_robinson_eos_gas_and_vapour_like_roots(
            pressure, 
            temperature, 
            acentricity_factor, 
            critical_pressure, 
            critical_temperature);

    // assert approx eq to within 0.1% 

    assert_relative_eq!(pr_eos_molar_vol.get::<cubic_meter_per_mole>(),
                        expected_steam_molar_vol.get::<cubic_meter_per_mole>(),
                        max_relative=0.001);

}

/// this function assumes there are two roots
///
/// TBD: change error type to a proper error
pub fn get_molar_volume_peng_robinson_eos_liquid_like_roots(
    pressure: Pressure,
    temperature: ThermodynamicTemperature,
    acentricity_factor: Ratio,
    critical_pressure: Pressure,
    critical_temperature: ThermodynamicTemperature) -> Result<MolarVolume,CubicEOSError>{
        
        let b = get_b(critical_pressure, critical_temperature);

        // R
        let molar_gas_constant = 
        MolarHeatCapacity::new::<joule_per_kelvin_mole>(8.314);

        // B
        let capital_b: Ratio = b * pressure / molar_gas_constant 
        / temperature;

        // A
        let capital_a: Ratio;

        let kappa = get_kappa(acentricity_factor);
        let reduced_temperature: Ratio 
        = temperature/critical_temperature;
        let sqrt_alpha = Ratio::new::<ratio>(1.0)
        + kappa * (Ratio::new::<ratio>(1.0)
        - reduced_temperature.sqrt());
        let alpha: Ratio = sqrt_alpha * sqrt_alpha;

        let a = 0.45724 * molar_gas_constant 
        * molar_gas_constant 
        * critical_temperature 
        * critical_temperature 
        / critical_pressure 
        * alpha;

        capital_a = a * pressure / molar_gas_constant 
        / molar_gas_constant
        / temperature
        / temperature;

        // Z^3 - (1-B) Z^2 + (A - 3B^2 - 2B) Z 
        // - (AB-B^2 -B^3)
        let one_minus_capital_b: Ratio 
        = Ratio::new::<ratio>(1.0) - capital_b;
        let capital_a_minus_3_cap_b_sq_minus_2_cap_b: Ratio 
        = capital_a - 3.0 * capital_b * capital_b 
        - 2.0 * capital_b;

        let cap_ab_minus_cap_b_sq_minus_b_cubed: Ratio 
        = capital_a * capital_b - 
        capital_b * capital_b 
        - capital_b * capital_b * capital_b;

        // find roots
        let cubic_root_b: f64 = -one_minus_capital_b.get::<ratio>();
        let cubic_root_c: f64 = capital_a_minus_3_cap_b_sq_minus_2_cap_b.
            get::<ratio>();
        let cubic_root_d: f64 = -cap_ab_minus_cap_b_sq_minus_b_cubed.
            get::<ratio>();

        let z_root_values: Vec<f64> = 
            find_cubic_roots(cubic_root_b,cubic_root_c,cubic_root_d);

        // case 1: one real root
        if z_root_values.len() == 1 {
            return Err(CubicEOSError::NoLiquidLikeRoots);
        }

        // case 2: three real roots
        if z_root_values.len() == 3 {

            // order the z roots first if there are three real roots
            
            let mut soln_vec_clone = z_root_values.clone();
            soln_vec_clone.sort_by(|a, b| a.partial_cmp(b).unwrap());

            // max root of Z is the maximum compressibility factor
            // this corresponds to vapour like roots
            // let max_root = *soln_vec_clone.last().unwrap();
            // max root of Z is the maximum compressibility factor
            // this corresponds to liquid like roots
            let min_root = *soln_vec_clone.first().unwrap();
            //
            //
            #[allow(non_snake_case)]
            let Z = Ratio::new::<ratio>(min_root);

            // molar volume 
            let v_tilde: MolarVolume = Z * molar_gas_constant * 
                temperature / pressure;

            return Ok(v_tilde);
        }
        


        panic!("unable to find liquid like volume roots");
}

// more often in Rust codebases, you will see error enums
use thiserror::Error;
#[derive(Error, Debug, PartialEq)]
pub enum CubicEOSError {
    #[error("no liquid like roots")]
    NoLiquidLikeRoots,
    #[error("no real temperature root")]
    NoRealTemperatureRoot,
}

fn get_kappa(acentricity_factor: Ratio) -> Ratio {

    let kappa: Ratio = 
    Ratio::new::<ratio>(0.37464) 
    + 1.54226 * acentricity_factor 
    - 0.26992 * acentricity_factor * acentricity_factor;

    return kappa;

}

fn get_b(critical_pressure: Pressure, 
    critical_temperature: ThermodynamicTemperature) 
    -> MolarVolume {

    // R
    let molar_gas_constant = 
    MolarHeatCapacity::new::<joule_per_kelvin_mole>(8.314);

    let b = 0.07780 * molar_gas_constant* critical_temperature/
    critical_pressure;

    b
}
