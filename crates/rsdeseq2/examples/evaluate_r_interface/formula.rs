use rsdeseq2::design::expression::{
    Argument, Data, EvaluationOptions, Expression, Value, evaluate,
};
use serde::Deserialize;
use serde_json::{Value as Json, json};
fn strings(v: &Json) -> Option<Vec<String>> {
    v.as_array()
        .map(|a| a.iter().map(|v| v.as_str().unwrap().into()).collect())
}
pub fn value(v: &Json) -> Value {
    let list = v["values"].as_array();
    let data = match v["type"].as_str().unwrap() {
        "number" => Data::Number {
            integer: v["integer"].as_bool().unwrap(),
            values: list
                .unwrap()
                .iter()
                .map(|v| f64::from_bits(u64::from_str_radix(v.as_str().unwrap(), 16).unwrap()))
                .collect(),
        },
        "text" => Data::Text(
            list.unwrap()
                .iter()
                .map(|v| v.as_str().map(str::to_owned))
                .collect(),
        ),
        "logical" => Data::Logical(list.unwrap().iter().map(Json::as_bool).collect()),
        "list" => Data::List(list.unwrap().iter().map(value).collect()),
        "null" => Data::Null,
        _ => panic!(),
    };
    let mut result = Value::new(data);
    result.names = strings(&v["names"]);
    result.dimnames = [strings(&v["dimnames"][0]), strings(&v["dimnames"][1])];
    result.dimname_names = [
        strings(&v["dimname_names"][0]),
        strings(&v["dimname_names"][1]),
    ];
    result.dimensions = v["dimensions"].as_array().map(|d| {
        [
            d[0].as_u64().unwrap() as usize,
            d[1].as_u64().unwrap() as usize,
        ]
    });
    if let Some(a) = v["attributes"].as_object() {
        result.attributes = a.iter().map(|(n, v)| (n.clone(), value(v))).collect()
    }
    result
}
pub fn expression(v: &Json) -> Expression {
    match v["type"].as_str().unwrap() {
        "symbol" => Expression::Symbol(v["name"].as_str().unwrap().into()),
        "literal" => Expression::Literal(value(&v["value"])),
        "call" => Expression::Call {
            function: v["function_name"].as_str().unwrap().into(),
            arguments: v["arguments"]
                .as_array()
                .unwrap()
                .iter()
                .map(|a| Argument {
                    name: a["name"].as_str().map(str::to_owned),
                    expression: expression(&a["expression"]),
                })
                .collect(),
        },
        _ => panic!(),
    }
}

pub fn descriptor(v: &Value) -> Json {
    let (kind, integer, values) = match &v.data {
        Data::Null => ("null", false, Json::Null),
        Data::Number { values, integer } => (
            "number",
            *integer,
            json!(
                values
                    .iter()
                    .map(|v| format!("{:016x}", v.to_bits()))
                    .collect::<Vec<_>>()
            ),
        ),
        Data::Logical(values) => ("logical", false, json!(values)),
        Data::Text(values) => ("text", false, json!(values)),
        Data::List(values) => (
            "list",
            false,
            json!(values.iter().map(descriptor).collect::<Vec<_>>()),
        ),
    };
    let attrs = if v.attributes.is_empty() {
        // R's reference encoder retains an empty named attribute container
        // when names/dimensions existed, although those are separate fields.
        if v.names.is_some() || v.dimensions.is_some() || v.dimnames.iter().any(Option::is_some) {
            json!({})
        } else {
            json!([])
        }
    } else {
        Json::Object(
            v.attributes
                .iter()
                .map(|(n, v)| (n.clone(), descriptor(v)))
                .collect(),
        )
    };
    let dimnames = if v.dimensions.is_some() && v.dimnames.iter().all(Option::is_none) {
        json!([])
    } else {
        json!(v.dimnames)
    };
    json!({"type":kind,"integer":integer,"values":values,"dimensions":v.dimensions,"names":v.names,"dimnames":dimnames,"dimname_names":v.dimname_names,"attributes":attrs})
}
#[derive(Deserialize)]
pub struct Input {
    formula_cases: Vec<Json>,
    variables: serde_json::Map<String, Json>,
    collation: String,
}
pub fn panel(input: Input) -> Json {
    let variables = input
        .variables
        .iter()
        .map(|(n, v)| (n.clone(), value(v)))
        .collect::<Vec<_>>();
    let options = EvaluationOptions {
        collation: input.collation,
    };
    json!({"cases":input.formula_cases.iter().map(|case| {
  match evaluate(&expression(&case["expression"]),&variables,&options) {
   Ok(value)=>json!({"name":case["name"],"ok":true,"value":descriptor(&value)}),
   Err(_)=>json!({"name":case["name"],"ok":false}),
  }
 }).collect::<Vec<_>>()})
}
