pub mod biclusters;
pub mod common;

use std::time::{Duration, Instant};
use std::{env, fs::File, process::Command};
use std::io::Write;

use biclusters::{algo::{bicluster, load_wadj_from_csv, print_wadj_stats}, biclust::Biclust, biclustering::Biclustering, weighted_biadj::WeightedBiAdjacency, r_results::load_r_biclusters};
use rand::Rng;

struct Args {
    size: f64,
    split: f64,
    power: usize,
    split_rows: bool,
    verbose: usize
}

impl Default for Args {
    fn default() -> Self {
        Args {
            size: 1.0,
            split: 1.0,
            power: 3,
            split_rows: true,
            verbose: 0
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
        wadj.write_to_file(&format!("{base_name}/{i}.edges") );
        let ground_truth = wadj.get_ground_truth();
        let (labels_a, labels_b, nodes_a_map, nodes_b_map) = wadj.get_labels();

        if let Some(ground_truth) = ground_truth {
            ground_truth.write_to_file( &format!("{base_name}/{i}.ground_truth"), Some((labels_a, labels_b)));
        }
    }
}


fn run_comparison(){
    // Comparison

    let mut file = File::create("comparison.csv").unwrap();

    writeln!(file, "n m real_noise p_noise real_overlap p_overlap p_separation CF_mtc CF_acc CF_time_s BM_mtc BM_acc").unwrap();

    for _ in 0..1000 {
        // let n = 10;
        // let m = 10;
        // let noise = 0.00;
        // let row_overlap = 1.02;
        // let row_separation = 0.8; 

        let mut rng = rand::thread_rng();

        let n = rng.gen_range(10..=100);
        let m = rng.gen_range(10..=100);
        let noise = rng.gen_range(0.0..=0.1);
        let row_overlap = rng.gen_range(1.0..=2.0);
        let row_separation = rng.gen_range(0.0..=1.0);


        let mut wadj = WeightedBiAdjacency::rand(n,m, noise, row_overlap, row_separation);
        // wadj.print();
        wadj.write_to_file("bigraphs/synth/test2.adj");
        wadj.write_to_file("gene.adj");
        let ground_truth = wadj.get_ground_truth();
        let wadj_save = wadj.clone(); // Because clustiflor is deleting edges in wadj
        
        // Clustiflor
        let start_time = Instant::now();
        let clusti_biclusters = bicluster(&mut wadj, 1., 1., 3, 0);
        let clustiflor_dur = start_time.elapsed();
        let (labels_a, labels_b, nodes_a_map, nodes_b_map) = wadj.get_labels();
        // biclusters.print_stats(1., 1., 3, &labels_a, &labels_b, Some(&"yo.biclusters"));

        // R Bimax
        Command::new("Rscript")
            .arg("./script.r")
            .status().unwrap();

        let r = load_r_biclusters( "r_results.txt", &nodes_a_map, &nodes_b_map);
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
            let bimax_accuracy = ground_truth.accuracy(&r);

            writeln!(file, "{n} {m} {real_noise:.4} {noise:.4} {real_overlap:.4} {row_overlap:.4} {row_separation:.4} {:.2} {clusti_accuracy:.2} {clusti_dur:.2} {:.2} {bimax_accuracy:.2}", 
            ground_truth.matching_score(&clusti_biclusters),
            ground_truth.matching_score(&r)).unwrap();

        }
    }
}



fn run_solver(){

    let program_args: Vec<String> = env::args().collect();

    let mut args = Args::default();


    for arg in program_args.iter() {
        if arg.starts_with("--size-sensivity") {
            args.size = arg.split_at(7).1.parse::<f64>().unwrap_or(args.size);
        } else if arg.starts_with("--split-th") {
            args.split = arg.split_at(8).1.parse::<f64>().unwrap_or(args.split);
            println!("split {}", args.split);
        } else if arg.starts_with("--power") {
            args.power = arg.split_at(8).1.parse::<usize>().unwrap_or(args.power);
        } else if arg.starts_with("--verbose") {
            args.verbose = arg.split_at(10).1.parse::<usize>().unwrap_or(args.verbose);
        } else if arg.starts_with("--split-cols") {
            args.split_rows = false;
        }
    }


    if program_args.len() < 3 {
        eprintln!("Usage: {} <file_path> <delimiter> --split-th=<number >= 1.>", program_args[0]);
        std::process::exit(1);
    }

    let file_path = &program_args[1];
    let delimiter = &program_args[2];

    println!("File path: {}", file_path);
    println!("Delimiter: {}", delimiter);

    match load_wadj_from_csv(file_path, delimiter, args.split_rows) {
        (wadj, n, m, labels_a, labels_b, node_map_a, node_map_b) => {
            print_wadj_stats(&wadj, n, m);
            let mut wadj2 = wadj.clone();
            let biclusters = bicluster(&mut wadj2, args.size, args.split, args.power, args.verbose);
            let results_path = file_path.to_string() + ".biclusters";
            biclusters.print_stats(args.size, args.split, args.power, &labels_a, &labels_b, Some(&results_path));


            
            // let r = load_r_biclusters( "results.txt", &node_map_a, &node_map_b);
            // for biclust in r.iter() {
            //     println!("{biclust:?}");
            // }

            // println!("{:?} {}", compute_nb_unclustered(&r, n, m), compute_edition_diff(&r, &wadj, n, m));

            
        }
    }

}

fn main() {

    
    gen_batch(3);
    
    // run_solver();
    // run_comparison();
    

    


    // let (matrix, node_map) = load_adj_list_file("real/kegg.txt", ' ');
    // // let matrix = load_edges_file("real/biodata_Shaw/output_file.txt", ' ');
    // // let matrix = generate_random_matrix(10, 0.6);    
    // // println!("{matrix:?}");

    // let n = matrix.shape()[0];
    // let mut m = 0;
    // for i in 0..n {
    //     for j in i+1..n{
    //         if matrix[[i,j]] > 0. {
    //             m += 1;
    //         }
    //     }
    // }
    // println!("n={n} m={m}");
    // let clusters = cluster_graph(matrix, true, 2, 1, 10.);


    // let mut reversed_map = HashMap::new();
    // for (key, value) in node_map.iter() {
    //     reversed_map.insert(*value, key.clone());
    // }

    // let output_file = File::create("clusters.txt").expect("Failed to open file");
    // let mut writer = BufWriter::new(output_file);

    // for cluster in &clusters {
    //     for v in cluster {
    //         write!(writer, "{} ", reversed_map.get(v).unwrap()).expect("Failed to write to file");
    //     }
    //     writeln!(writer).expect("Failed to write end-of-line to file");
    // }
}
