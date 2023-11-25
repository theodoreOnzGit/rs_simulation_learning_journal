use csv::Writer;
use uom::si::{f64::*, mass::{kilogram, gram}, specific_heat_capacity::joule_per_kilogram_kelvin, heat_transfer::watt_per_square_meter_kelvin, area::square_centimeter, thermodynamic_temperature::degree_celsius, time::second};

pub fn implicit_euler_method(){

    // solve dimensions
    let mass = Mass::new::<gram>(2.0);

    let cp = SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(3000.0);

    let heat_trf_coefficient = HeatTransfer::new::<watt_per_square_meter_kelvin>(20.0);
    let surface_area = Area::new::<square_centimeter>(5.0);

    let temp_surr = ThermodynamicTemperature::new::<degree_celsius>(25.0);

    let timestep = Time::new::<second>(0.1);

    let temp_solid_initial = ThermodynamicTemperature::new::<degree_celsius>(400.0);


    let mut solid_temperature_current_timestep: ThermodynamicTemperature = temp_solid_initial;
    let mut solid_temperature_next_timestep: ThermodynamicTemperature;

    let mut current_timestep = Time::new::<second>(0.0);

    // csv file 
    // use 
    // cargo watch -x run --ignore *.csv

    let mut wtr = Writer::from_path("newtons_law.csv").unwrap();

    let _a = 1;


    for _i in 0..5000{
        // introduce references and borrowing later
        //dbg!(&i);
        //println!("{}",i);

        // start writing csv file
        // write(current_timestep, solid_temperature_current_timestep)
        let current_timestep_string = current_timestep.get::<second>().to_string();
        let solid_temperature_current_timestep_string = 
        solid_temperature_current_timestep.get::<degree_celsius>().to_string();
        wtr.write_record(&[current_timestep_string, solid_temperature_current_timestep_string]).unwrap();

        // start calculations

        let denom_term_1 = mass * cp / timestep;
        let denom_term_2 = heat_trf_coefficient * surface_area;

        let denominator = denom_term_1 + denom_term_2;

        let temp_surr_degc_float = temp_surr.get::<degree_celsius>();

        let temp_surr_interval_deg_c = TemperatureInterval::new::<uom::si::temperature_interval::degree_celsius>(
            temp_surr_degc_float
        );

        let solid_temp_current_timestep_interval_deg_c = TemperatureInterval::new::<uom::si::temperature_interval::degree_celsius>(
            solid_temperature_current_timestep.get::<degree_celsius>()
        );

        let numerator_term_1 = heat_trf_coefficient * surface_area * temp_surr_interval_deg_c;
        let numerator_term_2 = mass * cp * solid_temp_current_timestep_interval_deg_c / timestep;

        // how to add 25 C to 400C?
        // 25 C + 400 C = 425 C

        // (25 + 273) + (400 +273) = 425 + 273 + 273

        let numerator = numerator_term_1 + numerator_term_2;

        solid_temperature_next_timestep = ThermodynamicTemperature::new::<degree_celsius>(
            (numerator/denominator).get::<uom::si::temperature_interval::degree_celsius>()
        );

        current_timestep = current_timestep + timestep;

        //println!("{}", current_timestep);
        //println!("{}",solid_temperature_next_timestep);

        solid_temperature_current_timestep = solid_temperature_next_timestep;

    }


}