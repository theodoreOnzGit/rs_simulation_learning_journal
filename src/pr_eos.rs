use uom::si::{f64::*, thermodynamic_temperature::{degree_celsius, kelvin}, pressure::pascal, molar_volume::cubic_meter_per_mole, molar_heat_capacity::joule_per_kelvin_mole, ratio::ratio};

// if you want a new custom unit 
// https://github.com/iliekturtles/uom/issues/244

///
pub fn get_temperature_peng_robinson_gas(
    pressure: Pressure,
    specific_vol: MolarVolume, 
    accentricity_factor: Ratio, 
    critical_temperature: ThermodynamicTemperature, 
    critical_pressure: Pressure) -> ThermodynamicTemperature{

        let b = get_b(critical_pressure, critical_temperature);
        let kappa = get_kappa(accentricity_factor);

        todo!()
}

pub fn get_pressure_peng_robinson_gas(
    temperature: ThermodynamicTemperature,
    molar_volume: MolarVolume,
    accentricity_factor: Ratio,
    critical_pressure: Pressure,
    critical_temperature: ThermodynamicTemperature) -> Pressure{
        
        let b = get_b(critical_pressure, critical_temperature);
        let kappa = get_kappa(accentricity_factor);

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
        let reduced_temeprature: Ratio 
        = temperature/critical_temperature;

        let sqrt_alpha = Ratio::new::<ratio>(1.0)
        + kappa * (Ratio::new::<ratio>(1.0)
        - reduced_temeprature.sqrt());

        let alpha = sqrt_alpha * sqrt_alpha;

        let numerator = 0.45724 * molar_gas_constant 
        * molar_gas_constant 
        * critical_temperature 
        * critical_temperature 
        / critical_pressure 
        * alpha;

        let accentricity_term: Pressure = 
        numerator/denominator;


        let pressure = rt_over_v_minus_b - accentricity_term;
        return pressure;
}

pub fn get_molar_volume_peng_robinson_gas(
    pressure: Pressure,
    temperature: ThermodynamicTemperature,
    accentricity_factor: Ratio,
    critical_pressure: Pressure,
    critical_temperature: ThermodynamicTemperature) -> MolarVolume{
        
        let b = get_b(critical_pressure, critical_temperature);
        let kappa = get_kappa(accentricity_factor);



        todo!()
}

fn get_kappa(accentricity_factor: Ratio) -> Ratio {

    let kappa: Ratio = 
    Ratio::new::<ratio>(0.37464) 
    + 1.54226 * accentricity_factor 
    - 0.26992 * accentricity_factor * accentricity_factor;

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