use rs_edlib::{align::Alignment, config::AlignConfig};

#[test]
fn test_usage() {
    let query: &str = "hello";
    let target: &str = "world!";
    let align_res = Alignment::run(
        AlignConfig::default(),
        query,
        target
    ).unwrap();
    println!("edit_distance('{query}', '{target}') = {:?}", align_res.edit_distance)
}