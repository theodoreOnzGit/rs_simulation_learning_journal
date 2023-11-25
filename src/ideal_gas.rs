use uom::si::{f64::*, thermodynamic_temperature::{degree_celsius, kelvin}, pressure::pascal, molar_volume::cubic_meter_per_mole, molar_heat_capacity::joule_per_kelvin_mole};


pub fn get_temperature_ideal_gas(
    pressure: Pressure,
    specific_vol: MolarVolume) -> ThermodynamicTemperature{

    // R is in joule/mole/kelvin
    let molar_gas_constant = MolarHeatCapacity::new::<joule_per_kelvin_mole>(8.314);

    let abs_temp_interval: 
    TemperatureInterval = pressure * specific_vol/molar_gas_constant;

    let abs_temp_kelvin_value: f64 = 
    abs_temp_interval.get::<uom::si::temperature_interval::kelvin>();

    let temp = ThermodynamicTemperature::new::<kelvin>(abs_temp_kelvin_value);

    return temp;
}