//! Primitive adapter for the shared standalone frozen dispersion implementation.
use super::{Data, Value, ast, fields, get, text};
use rsdeseq2::transform::{FrozenLocalDispersion, SavedLocalDispersionFit};

fn numbers(v: &Value) -> Result<&[f64], String> {
    match &v.data {
        Data::Number { values, .. } => Ok(values),
        _ => Err("frozen local inputs must be numeric".into()),
    }
}
fn model(v: &Value) -> Result<FrozenLocalDispersion, String> {
    match text(get(v, "type")?)? {
        "saved" => {
            let mi = get(v, "mi")?;
            let names = mi
                .names
                .as_ref()
                .ok_or("saved local fit requires named profile fields")?;
            let values = numbers(mi)?;
            for (name, expected) in [
                ("p", 3.),
                ("deg0", 2.),
                ("deg", 2.),
                ("d", 1.),
                ("ker", 4.),
                ("kt", 1.),
                ("it", 1.),
                ("ev", 1.),
                ("tg", 68.),
                ("link", 3.),
                ("pc", 1.),
                ("ubas", 0.),
            ] {
                if names.iter().filter(|n| n.as_str() == name).count() != 1
                    || names
                        .iter()
                        .position(|n| n == name)
                        .is_none_or(|i| values[i] != expected)
                {
                    return Err(format!("unsupported saved local fit profile: {name}"));
                }
            }
            if numbers(get(v, "scale")?)? != [1.]
                || numbers(get(v, "style")?)? != [1.]
                || !numbers(get(v, "derivative")?)?.is_empty()
            {
                return Err("unsupported saved local fit scale, style or derivative".into());
            }
            let knots = numbers(get(v, "knots")?)?.to_vec();
            let coefficients = get(v, "coefficients")?;
            if coefficients.dimensions != Some([knots.len(), 11]) {
                return Err(
                    "saved local coefficients must have one eleven-column row per knot".into(),
                );
            }
            let coefficients = numbers(coefficients)?;
            let polynomial = numbers(get(v, "polynomial")?)?;
            let polynomial: [f64; 4] = polynomial
                .get(..4)
                .ok_or("saved local global polynomial is incomplete")?
                .try_into()
                .unwrap();
            SavedLocalDispersionFit::new(
                knots,
                coefficients
                    .as_chunks::<11>()
                    .0
                    .iter()
                    .map(|x| x[0])
                    .collect(),
                coefficients
                    .as_chunks::<11>()
                    .0
                    .iter()
                    .map(|x| x[1])
                    .collect(),
                polynomial,
            )
            .map(FrozenLocalDispersion::Saved)
            .map_err(|e| e.to_string())
        }
        "expression" => {
            let (names, constants) = fields(get(v, "constants")?)?;
            Ok(FrozenLocalDispersion::Expression {
                argument: text(get(v, "argument")?)?.into(),
                expression: ast(get(v, "expression")?, 0)?,
                constants: names
                    .iter()
                    .cloned()
                    .zip(constants.iter().cloned())
                    .collect(),
            })
        }
        _ => Err("unsupported frozen local dispersion representation".into()),
    }
}
pub(super) fn evaluate(request: &Value) -> Result<Value, String> {
    let mut grid = numbers(get(request, "means")?)?.to_vec();
    if grid.len() != 999
        || grid.iter().any(|v| !v.is_finite() || *v < 0.0)
        || grid.windows(2).any(|x| x[0] > x[1])
    {
        return Err("frozen local means require the finite increasing 999-point grid".into());
    }
    let dispersions = model(get(request, "model")?)?
        .evaluate(&grid)
        .map_err(|e| e.to_string())?;
    grid.extend(dispersions);
    Ok(Value::numbers(grid))
}
