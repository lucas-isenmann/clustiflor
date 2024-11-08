pub mod biclusters;

use std::env;

use biclusters::{biclust::Biclust, biclustering::Biclustering, r_results::load_r_biclusters};

struct Args {
    cost: f64,
    split: f64,
    power: usize,
    verbose: usize
}

impl Default for Args {
    fn default() -> Self {
        Args {
            cost: 1.0,
            split: 1.0,
            power: 3,
            verbose: 0
        }
    }
}


fn main() {

    match load_wadj_from_csv(file_path, delimiter) {
            (wadj, n, m, labels_a, labels_b, node_map_a, node_map_b) => {

            }

    }

    let a = load_r_biclusters("r_results.txt", node_map_a, node_map_b)
    
    // let program_args: Vec<String> = env::args().collect();

    // let mut args = Args::default();


    // for arg in program_args.iter() {
    //     if arg.starts_with("--cost") {
    //         args.cost = arg.split_at(7).1.parse::<f64>().unwrap_or(args.cost);
    //     } else if arg.starts_with("--split") {
    //         args.split = arg.split_at(8).1.parse::<f64>().unwrap_or(args.split);
    //         println!("split {}", args.split);
    //     } else if arg.starts_with("--power") {
    //         args.power = arg.split_at(8).1.parse::<usize>().unwrap_or(args.power);
    //     } else if arg.starts_with("--verbose") {
    //         args.verbose = arg.split_at(10).1.parse::<usize>().unwrap_or(args.verbose);
    //     }
    // }


    // if program_args.len() < 3 {
    //     eprintln!("Usage: {} <file_path> <delimiter> --split=<>", program_args[0]);
    //     std::process::exit(1);
    // }

    // let file_path = &program_args[1];
    // let delimiter = &program_args[2];

    // println!("File path: {}", file_path);
    // println!("Delimiter: {}", delimiter);

    // match load_wadj_from_csv(file_path, delimiter) {
    //     (wadj, n, m, labels_a, labels_b, node_map_a, node_map_b) => {
    //         print_wadj_stats(&wadj, n, m);
    //         let mut wadj2 = wadj.clone();
    //         let biclusters = bicluster(&mut wadj2, n, m, args.cost, args.split, args.power, args.verbose);
    //         let results_path = file_path.to_string() + ".biclusters";
    //         print_biclusters_stats(&biclusters, args.cost, args.split, args.power, n, &labels_a, &labels_b, Some(&results_path));


            
    //         // let r = load_r_biclusters( "results.txt", &node_map_a, &node_map_b);
    //         // for biclust in r.iter() {
    //         //     println!("{biclust:?}");
    //         // }

    //         // println!("{:?} {}", compute_nb_unclustered(&r, n, m), compute_edition_diff(&r, &wadj, n, m));

            
    //     }
    // }



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
