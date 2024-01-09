use uom::ConstZero;
use uom::si::f64::*;

use crate::pr_eos::CubicEOSError;

pub fn find_real_quadratic_roots(
    a: Ratio,
    b: Ratio,
    c: Ratio) -> Result<[Ratio;2], CubicEOSError> {

    let bsq_mins_4ac: Ratio = b*b - 4.0 * a * c;

    if bsq_mins_4ac < Ratio::ZERO {
        return Err(CubicEOSError::NoRealTemperatureRoot);
    };

    let x1: Ratio = (-b + bsq_mins_4ac.sqrt())/(2.0 * a);
    let x2: Ratio = (-b - bsq_mins_4ac.sqrt())/(2.0 * a);

    let soln_array = [x1,x2];

    return Ok(soln_array);


}
