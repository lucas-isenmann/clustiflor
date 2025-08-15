pub mod biclusters;
pub mod common;
pub mod clustering;

use std::path::Path;
use std::time::{Instant};
use std::{env, fs::File, process::Command};
use std::io::{BufRead, BufReader, Write};

use biclusters::algo::load_wadj_from_matrix;
use biclusters::algo_two_sided::{analyze_ground_biclusters, bicluster_two_sided};
use biclusters::{algo::{bicluster_one_sided, load_wadj_from_csv, print_wadj_stats}, biclustering::Biclustering, weighted_biadj::WeightedBiAdjacency, r_results::load_r_biclusters};
use rand::Rng;
use walkdir::WalkDir;

use crate::clustering::cluster_algo::{run_cluster_solver};

struct Args {
    size: f64,
    split: f64,
    power: usize,
    split_rows: bool,
    verbose: usize,
    matrix: bool
}

impl Default for Args {
    fn default() -> Self {
        Args {
            size: 1.0,
            split: 1.0,
            power: 3,
            split_rows: true,
            verbose: 0,
            matrix: false
        }
    }
}

fn gen_batch_v2(batch_size: usize){
    let base_name = "bigraphs/synth/batch2";
    for i in 0..batch_size {

        let mut rng = rand::thread_rng();

        let n = rng.gen_range(10..=400);
        let m = rng.gen_range(10..=400);
        let c = rng.gen_range(1..(m/2));
        let over = rng.gen_range(0.0..0.3);
        let noise = rng.gen_range(0.0..=0.001);

        println!("n {n} m {m} c {c} over {over} noise {noise}");



        let wadj = WeightedBiAdjacency::rand_v2(n,m, c, over, noise);

        let row_degrees_distrib = wadj.row_degrees_distributon();
        let col_degrees_distrib = wadj.col_degrees_distributon();

        println!("{:?}", row_degrees_distrib);
        println!("{:?}", col_degrees_distrib);
        // wadj.print_matrix();
        
        if let Some(gt) = wadj.get_ground_truth(){
            gt.print();
        }


        // let (biclusters, stats) = bicluster(&mut wadj.clone(), 1.0, 1.0, 3, 0);


        wadj.write_to_file(&format!("{base_name}/{i}.edges"), &format!("# n={n} m={m} c={c} over={over:.3} noise={noise:.3}") );

        let ground_truth = wadj.get_ground_truth();
        let (labels_a, labels_b, nodes_a_map, nodes_b_map) = wadj.get_labels();

        if let Some(ground_truth) = ground_truth {
            ground_truth.write_to_file( &format!("{base_name}/{i}.ground_truth"), Some((labels_a, labels_b)));
        }
    }
}



fn gen_batch(batch_size: usize){
    let base_name = "bigraphs/synth/batch";
    for i in 0..batch_size {

        let mut rng = rand::thread_rng();

        let n = rng.gen_range(10..=100);
        let m = rng.gen_range(10..=100);
        let noise = rng.gen_range(0.0..=0.1);
        let row_overlap = rng.gen_range(1.0..=2.0);
        let row_separation = rng.gen_range(0.0..=1.0);


        let wadj = WeightedBiAdjacency::rand(n,m, noise, row_overlap, row_separation);

        wadj.write_to_file(&format!("{base_name}/{i}.edges"), &format!("# n={n} m={m} noise={noise:.3} row_overlap={row_overlap:.3} row_separation={row_separation:.3}") );

        let ground_truth = wadj.get_ground_truth();
        let (labels_a, labels_b, nodes_a_map, nodes_b_map) = wadj.get_labels();

        if let Some(ground_truth) = ground_truth {
            ground_truth.write_to_file( &format!("{base_name}/{i}.ground_truth"), Some((labels_a, labels_b)));
        }
    }
}

pub fn read_duration_file(
    file_path: &str ) -> f64 {

    let file = File::open(file_path).expect("Failed to open file");
    let reader = BufReader::new(file);

    for line in reader.lines() {
        match line {
            Ok(line) => {
                println!("read: {line}");
                let line = line.trim();
                if let Ok(duration) = line.parse::<f64>() {
                    println!("{duration}");
                    return duration;
                }
            },
            Err(_) => {},
        }
    }
    println!("cant read {file_path}");
    0.
}


fn run_comparison(){
    // Comparison

    let mut file = File::create("comparison.csv").unwrap();

    writeln!(file, "n m real_noise p_noise real_overlap p_overlap p_separation CF_mtc CF_acc CF_time_s BM_mtc BM_acc").unwrap();

    for _ in 0..10000 {
        // let n = 10;
        // let m = 10;
        // let noise = 0.00;
        // let row_overlap = 1.02;
        // let row_separation = 0.8; 

        let mut rng = rand::thread_rng();

        let n = rng.gen_range(10..=140);
        let m = rng.gen_range(10..=140);
        let noise = rng.gen_range(0.0..=0.2);
        let row_overlap = rng.gen_range(1.0..=2.0);
        let row_separation = rng.gen_range(0.0..=1.0);


        let mut wadj = WeightedBiAdjacency::rand(n,m, noise, row_overlap, row_separation);
        // wadj.print();
        wadj.write_to_file("bigraphs/synth/test2.adj", "#");
        wadj.write_to_file("gene.adj", "#");
        let ground_truth = wadj.get_ground_truth();
        let wadj_save = wadj.clone(); // Because clustiflor is deleting edges in wadj
        
        // Clustiflor
        let start_time = Instant::now();
        let (clusti_biclusters, clusti_stats) = bicluster_one_sided(&mut wadj, 1., 1., 3, 0);
        let clustiflor_dur = start_time.elapsed();
        let (labels_a, labels_b, nodes_a_map, nodes_b_map) = wadj.get_labels();
        // biclusters.print_stats(1., 1., 3, &labels_a, &labels_b, Some(&"yo.biclusters"));

        // R Bimax
        Command::new("Rscript")
            .arg("./script.r")
            .status().unwrap();

        let bimax_results = load_r_biclusters( "bimax_results.txt", &nodes_a_map, &nodes_b_map);
        let bimax_duration = read_duration_file("bimax_duration.txt");

        
        // Bibit
        Command::new("python3")
            .arg("./bibit2.py")
            .status().unwrap();

        let bibit_results = load_r_biclusters( "bibit_results.txt", &nodes_a_map, &nodes_b_map);
        let bibit_duration = read_duration_file("bibit_duration.txt");
        
        if let Some(ground_truth) = ground_truth {
            // println!("ground truth");
            // ground_truth.print();
            // println!("biclussters");
            // biclusters.print();
            // println!("R");
            // r.print();
            // println!("{}", ground_truth.matching_score(&biclusters));
            // println!("{}", ground_truth.matching_score(&r));
            
            let real_noise = wadj_save.compute_noise(&ground_truth);
            let clusti_dur = clustiflor_dur.as_secs_f64();
            // let clusti_fscore = ground_truth.f_score(&clusti_biclusters);
            let clusti_accuracy = ground_truth.accuracy(&clusti_biclusters);
            let real_overlap = ground_truth.get_rows_overlapping();
            let bimax_accuracy = ground_truth.accuracy(&bimax_results);

            writeln!(file, "{n} {m} {real_noise:.4} {noise:.4} {real_overlap:.4} {row_overlap:.4} {row_separation:.4} {:.2}  {clusti_dur:.2} {:.2} {bimax_duration:.4} {:.2} {bibit_duration}", 
            ground_truth.matching_score(&clusti_biclusters),
            ground_truth.matching_score(&bimax_results),
            ground_truth.matching_score(&bibit_results)
        ).unwrap();

        }
    }
}



fn run_bicluster_one_sided(){

    let program_args: Vec<String> = env::args().collect();

    let mut args = Args::default();


    for arg in program_args.iter() {
        if arg.starts_with("--size-sensitivity") {
            args.size = arg.split_at(19).1.parse::<f64>().unwrap_or(args.size);
        } else if arg.starts_with("--split-th") {
            args.split = arg.split_at(11).1.parse::<f64>().unwrap_or(args.split);
            println!("split {}", args.split);
        } else if arg.starts_with("--power") {
            args.power = arg.split_at(8).1.parse::<usize>().unwrap_or(args.power);
        } else if arg.starts_with("--verbose") {
            args.verbose = arg.split_at(10).1.parse::<usize>().unwrap_or(args.verbose);
        } else if arg.starts_with("--split-cols") {
            args.split_rows = false;
        } else if arg.starts_with("--matrix") {
            args.matrix = true;
        }
    }


    if program_args.len() < 4 {
        eprintln!("Usage: {} bicluster <file_path> <delimiter> --split-th=<number >= 1.>", program_args[0]);
        std::process::exit(1);
    }

    let file_path = &program_args[2];
    let delimiter = &program_args[3];

    println!("File path: {}", file_path);
    println!("Delimiter: {}", delimiter);

    if args.matrix {
        let wadj =  load_wadj_from_matrix(file_path);
        print_wadj_stats(&wadj, wadj.get_n(), wadj.get_m());
        let mut wadj2 = wadj.clone();
        let (biclusters, algo_stats) = bicluster_one_sided(&mut wadj2, args.size, args.split, args.power, args.verbose);
        let results_path = file_path.to_string() + ".biclusters";
        let mut labels_a = vec![];
        let mut labels_b = vec![];
        for line in 0..wadj.get_n() {
            labels_a.push(line.to_string())
        } 
        for j in 0..wadj.get_m() {
            labels_b.push(format!("c{j}"));
        }
        biclusters.print_stats(args.size, args.split, args.power, &labels_a, &labels_b, Some(&results_path), algo_stats);


    } else {
        match load_wadj_from_csv(file_path, delimiter, args.split_rows) {
            (wadj, n, m, labels_a, labels_b, node_map_a, node_map_b) => {
                print_wadj_stats(&wadj, n, m);
                let mut wadj2 = wadj.clone();
                let (biclusters, algo_stats) = bicluster_one_sided(&mut wadj2, args.size, args.split, args.power, args.verbose);
                let results_path = file_path.to_string() + ".biclusters";
                biclusters.print_stats(args.size, args.split, args.power, &labels_a, &labels_b, Some(&results_path), algo_stats);
    
    
                
                // let r = load_r_biclusters( "results.txt", &node_map_a, &node_map_b);
                // for biclust in r.iter() {
                //     println!("{biclust:?}");
                // }
    
                // println!("{:?} {}", compute_nb_unclustered(&r, n, m), compute_edition_diff(&r, &wadj, n, m));
    
                
            }
        }
    }

    

}





fn run_bicluster_solver(){

    let program_args: Vec<String> = env::args().collect();

    let mut args = Args::default();


    for arg in program_args.iter() {
        if arg.starts_with("--size-sensitivity") {
            args.size = arg.split_at(19).1.parse::<f64>().unwrap_or(args.size);
        } else if arg.starts_with("--split-th") {
            args.split = arg.split_at(11).1.parse::<f64>().unwrap_or(args.split);
            println!("split {}", args.split);
        } else if arg.starts_with("--power") {
            args.power = arg.split_at(8).1.parse::<usize>().unwrap_or(args.power);
        } else if arg.starts_with("--verbose") {
            args.verbose = arg.split_at(10).1.parse::<usize>().unwrap_or(args.verbose);
        }
    }


    if program_args.len() < 2 {
        eprintln!("Usage: {} <data_path> --split-th=<number >= 1.>", program_args[0]);
        std::process::exit(1);
    }

    let data_path = &program_args[1];

    println!("Data path: {}", data_path);


    // Check if path exists
    let path = Path::new(data_path);
    if !path.exists() {
        eprintln!("Error: Path '{}' does not exist", data_path);
        std::process::exit(1);
    }
    
    // Handle directory case
    if path.is_dir() {
        process_directory(data_path, &args);
    } else {
        // Handle single file case
        process_single_file(data_path, &args);
    }


    // let wadj =  load_wadj_from_matrix(data_path);
    // print_wadj_stats(&wadj, wadj.get_n(), wadj.get_m());
    // let mut wadj2 = wadj.clone();
    // let (biclusters, algo_stats) = bicluster_two_sided(&mut wadj2, args.size, args.split, args.power, args.verbose);
    

    // biclusters.print_biclusters( Some(&(data_path.to_string() + ".BiMarkov.results")));


    

}



fn process_directory(dir_path: &str, args: &Args) {
    println!("Processing directory: {}", dir_path);
    
    for entry in WalkDir::new(dir_path)
        .into_iter()
        .filter_map(|e| e.ok())
        .filter(|entry| {
            entry.file_type().is_file() &&
            entry.path().extension().map_or(false, |ext| ext == "data")
        }) {
            
        let file_path = entry.path();
        println!("Processing file: {:?}", file_path);
        
        let wadj = load_wadj_from_matrix(file_path.to_str().unwrap());
        print_wadj_stats(&wadj, wadj.get_n(), wadj.get_m());
        
        if let Some(bc) =wadj.get_ground_truth(){
            println!("{:?}", bc.biclusters());
        }
        
        let mut wadj2 = wadj.clone();
        let (biclusters, algo_stats) = bicluster_two_sided(
            &mut wadj2, 
            1., 
            1., 
            args.power, 
            args.verbose
        );

        if let Some(ground_biclusters) = wadj.get_ground_truth(){
            let ms = ground_biclusters.matching_score(&biclusters);
            let ms2 = biclusters.matching_score(&ground_biclusters);
            println!("matching score: {ms:.3} {ms2:.3} {:.3}", (ms*ms2).sqrt()  );
        }
        
        let result_path = format!("{}.BiMarkov.results", file_path.display());
        biclusters.print_biclusters(Some(&result_path));
    }
}

fn process_single_file(file_path: &str, args: &Args) {
    println!("Processing single file: {}", file_path);
    
    let wadj = load_wadj_from_matrix(file_path);
    print_wadj_stats(&wadj, wadj.get_n(), wadj.get_m());

    if let Some(bc) =wadj.get_ground_truth(){
        println!("{:?}", bc.biclusters());
    }
    
    let mut wadj2 = wadj.clone();
    let (biclusters, algo_stats) = bicluster_two_sided(
        &mut wadj2, 
        1.0, 
        args.split, 
        args.power, 
        args.verbose
    );

    if let Some(ground_biclusters) = wadj.get_ground_truth(){
        let ms = ground_biclusters.matching_score(&biclusters);
        let ms2 = biclusters.matching_score(&ground_biclusters);
        println!("matching score: {ms:.3} {ms2:.3} {:.3}", (ms*ms2).sqrt()  );
    }

    analyze_ground_biclusters(&wadj);
    
    biclusters.print_biclusters(Some(&(file_path.to_string() + ".BiMarkov3.results")));
}








fn main() {

    let program_args: Vec<String> = env::args().collect();

    if program_args.len() < 3 || !(program_args[1] == "cluster" || program_args[1] == "bicluster") {
        eprintln!("Usage: {} cluster|bicluster <data_path>", program_args[0]);
        std::process::exit(1);
    }

    if program_args[1] == "cluster" {
        run_cluster_solver();
    }
    else if program_args[1] == "bicluster" {
        run_bicluster_one_sided();
    }


    return;

    // run_bicluster_solver();
    
    // gen_batch_v2(1);
    
    // run_solver();
    // run_comparison();
    
}
