extern crate clap;
extern crate shrust;
use clap::{arg, Command};
use std::collections::{HashMap, HashSet};
use bio::io::{fasta,  fastq};
use generalized_suffix_tree::suffix_tree::KGST;
use indicatif::{ProgressBar, ProgressStyle};
use serde::{Serialize, Deserialize};
use std::{fs, fmt};
use std::io::Write;
use error_chain::error_chain;
use std::time::{Duration, Instant};
use glob::{glob_with, MatchOptions};
use rayon::prelude::*;

error_chain! {
    foreign_links {
        Glob(glob::GlobError);
        Pattern(glob::PatternError);
    }
}

#[derive(Eq, Hash, PartialEq, Copy, Clone, Debug, Serialize, Deserialize)]
enum SeqElement {
    A, G, T, C, E
}

impl fmt::Display for SeqElement {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            SeqElement::A => write!(f, "A"),
            SeqElement::G => write!(f, "G"),
            SeqElement::T => write!(f, "T"),
            SeqElement::C => write!(f, "C"),
            SeqElement::E => write!(f, "$"),
        }
    }
}

fn get_files(dir:&str) -> Vec<String>{
    let mut files: Vec<String> = Vec::new();
    let options = MatchOptions {
        case_sensitive: false,
        ..Default::default()
    };
    for entry in glob_with(format!("{}/**/single*.fastq", dir).as_str(), options).unwrap() {
        files.push(format!("{}", entry.expect("file").display()));
    }
    files
}

fn get_files_all(dir:&str) -> Vec<String>{
    let mut files: Vec<String> = Vec::new();
    let options = MatchOptions {
        case_sensitive: false,
        ..Default::default()
    };
    for entry in glob_with(format!("{}/**/*.fastq", dir).as_str(), options).unwrap() {
        files.push(format!("{}", entry.expect("file").display()));
    }
    files
}


fn build_tree(file:&str, max_depth:usize, num_seq: u32)->KGST<SeqElement, String>{
    println!("Building tree from {}", file);
    let reader = fasta::Reader::from_file(file).unwrap();

    let total_size = reader.records().count();

    let pb = ProgressBar::new(total_size as u64);
    pb.set_style(ProgressStyle::with_template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {bytes}/{total_bytes} ({eta})")
        .unwrap()
        .progress_chars("#>-"));
    
    let mut tree: KGST<SeqElement, String> = KGST::new(SeqElement::E);

    let reader = fasta::Reader::from_file(file).unwrap();
    let mut strings:HashMap<String, Vec<SeqElement>> = HashMap::new();

    let mut count = 0;
    
    for result in reader.records() {

        let result_data = result.unwrap();

        let x: Vec<char> = result_data.seq()
        .to_vec()
        .iter()
        .map(|x| *x as char)
        .collect();
        
        let seq: Vec<SeqElement> = x.iter()
            .map(|x|{
                match x{
                    'A' => SeqElement::A,
                    'G' => SeqElement::G,
                    'T' => SeqElement::T,
                    'C' => SeqElement::C,
                    _ => SeqElement::E,
                }
            })
            .collect();
        let string_len = seq.len();
    
        if string_len>max_depth.try_into().unwrap(){
            let num_iter = string_len+1-(max_depth);
            for (n, depth) in (0..num_iter).enumerate(){
                tree.add_string(seq[depth..depth+(max_depth)].to_vec(), format!("{}___{}", result_data.id(), n));
            }
            
        }
        strings.insert(format!("{}", result_data.id()), seq);
        pb.inc(1);   
        count+=1;
        if count==num_seq {
            break;
        }
    }
    tree.set_strings(strings);
    tree
}

fn match_prob(q_seq_match_vec: Vec<bool>, quality_score_vec: Vec<u8>) -> f32{
    return 0.0;
}

fn complement(q_seq: Vec<SeqElement>)->Vec<SeqElement>{
    q_seq.iter()
        .map(|x|{
            match x{
                SeqElement::T => SeqElement::A,
                SeqElement::C => SeqElement::G,
                SeqElement::A => SeqElement::T,
                SeqElement::G => SeqElement::C,
                _ => SeqElement::E,
            }
        })
        .collect()
}

fn preprocess_read(read: &[u8])->Option<Vec<SeqElement>>{
    let x: Vec<char> = read
        .to_vec()
        .iter()
        .map(|x| *x as char)
        .collect();
    
        if x.contains(&'N'){
            return None;
        }
        else{
            let seq: Vec<SeqElement> = x.iter()
            .map(|x|{
                match x{
                    'A' => SeqElement::A,
                    'G' => SeqElement::G,
                    'T' => SeqElement::T,
                    'C' => SeqElement::C,
                    _ => SeqElement::E,
                }
            })
            .collect();
            return Some(seq);
        }
}

fn hamming_distance(ref_seq: &Vec<SeqElement>, read_seq: &Vec<SeqElement>)->(usize, Vec<bool>){
    let count = ref_seq.iter().zip(read_seq.iter()).filter(|(x, y)| x!=y).count();
    let match_vec = ref_seq.iter().zip(read_seq.iter()).map(|(x, y)| {
        if x!=y{
            return false;
        }
        else{
            return true;
        }
    }).collect();
    (count, match_vec)
}

fn score(match_seq: Vec<bool>) -> usize{
    let mut total: usize = 0;
    let mut count:usize = 0;
    for el in match_seq.iter(){
        match el{
            &true => count +=1,
            &false => {
                if count>=22{
                    total += (count-15).pow(2);
                }
                count = 0;
            }
        }
    }
    if count>=15{
        total += (count-15).pow(2);
    }
    total
}

fn check_pos(ref_seq_len: &usize, depth: &usize, ref_seq_pos: &usize, q_len: &usize)->bool{
    depth<=ref_seq_pos && (q_len-depth)<=(ref_seq_len-ref_seq_pos)
}

fn query_tree(tree:&KGST<SeqElement, String>, q_seq:Vec<SeqElement>, q_seq_id:String, percent_match: f32)->Vec<(String, String, usize, usize)>{
    let mut match_set:Vec<(String, String, usize, usize)> = Vec::new();    //ref seq id, read id, score, start_pos
    let string_len: usize = q_seq.len();
    let mismatch_lim: usize = (q_seq.len() as f32 * percent_match).floor() as usize;
    let chunk_size: usize = string_len/(mismatch_lim+1);
    if string_len>=chunk_size{
        // (0..string_len+1-(chunk_size)).into_par_iter().for_each_with(match_set.clone(), |x, depth| {
        //     let mut lmer_matches: Vec<(String, String, usize, usize)> = tree.find(q_seq[depth..depth+(chunk_size)].to_vec()).iter()
        //                                                              .map(|(ref_id, ref_pos)| {
        //                                                                  (ref_id.split("___").collect::<Vec<&str>>()[0].to_string(), ref_id.split("___").collect::<Vec<&str>>()[1].parse().unwrap())
        //                                                              })
        //                                                              .filter(|(ref_id, ref_seq_pos)| check_pos(&tree.get_string(&(ref_id)).len(), &depth, ref_seq_pos, &string_len))
        //                                                              .map(|(ref_id, ref_seq_pos)| {
        //                                                                 (ref_id.clone(), q_seq_id.clone(), score(hamming_distance(&tree.get_string(&ref_id)[ref_seq_pos-depth..ref_seq_pos-depth+string_len].to_vec(), &q_seq).1), ref_seq_pos)
        //                                                              })
        //                                                              .collect();
        //     x.append(&mut lmer_matches);
        // });
        for depth in (0..string_len+1-(chunk_size)) {
            match_set.par_extend(tree.find(q_seq[depth..depth+(chunk_size)].to_vec()).par_iter()
            .map(|(ref_id, ref_pos)| {
                (ref_id.split("___").collect::<Vec<&str>>()[0].to_string(), ref_id.split("___").collect::<Vec<&str>>()[1].parse().unwrap())
            })
            .filter(|(ref_id, ref_seq_pos)| check_pos(&tree.get_string(&(ref_id)).len(), &depth, ref_seq_pos, &string_len))
            .map(|(ref_id, ref_seq_pos)| {
               (ref_id.clone(), q_seq_id.clone(), score(hamming_distance(&tree.get_string(&ref_id)[ref_seq_pos-depth..ref_seq_pos-depth+string_len].to_vec(), &q_seq).1), ref_seq_pos)
            }));
            let mut lmer_matches: Vec<(String, String, usize, usize)> = tree.find(q_seq[depth..depth+(chunk_size)].to_vec()).iter()
                                                                     .map(|(ref_id, ref_pos)| {
                                                                         (ref_id.split("___").collect::<Vec<&str>>()[0].to_string(), ref_id.split("___").collect::<Vec<&str>>()[1].parse().unwrap())
                                                                     })
                                                                     .filter(|(ref_id, ref_seq_pos)| check_pos(&tree.get_string(&(ref_id)).len(), &depth, ref_seq_pos, &string_len))
                                                                     .map(|(ref_id, ref_seq_pos)| {
                                                                        (ref_id.clone(), q_seq_id.clone(), score(hamming_distance(&tree.get_string(&ref_id)[ref_seq_pos-depth..ref_seq_pos-depth+string_len].to_vec(), &q_seq).1), ref_seq_pos)
                                                                     })
                                                                     .collect();
            match_set.append(&mut lmer_matches);

        }

    }
    match_set
}

fn save_tree(tree: &mut KGST<SeqElement, String>, output_path: String){
    println!("Saving tree to {}.", &output_path);
    std::fs::write(
        output_path,
        serde_json::to_string_pretty(tree).unwrap(),
    ).unwrap();
    println!("Saved");
}

fn load_tree(fname:&String) -> KGST<SeqElement, String>{
    println!("Loading tree from {}", fname);
    let json_str:String = fs::read_to_string(fname).unwrap();
    let tree: KGST<SeqElement, String> = serde_json::from_str(&json_str).unwrap();
    tree
}

fn search_fastq(tree:&KGST<SeqElement, String>, fastq_file:&str, result_file:&str, percent_match:f32){
    println!("Classifying read file: {}", &fastq_file);
    let reader = fastq::Reader::from_file(fastq_file).unwrap();

    let mut match_set:Vec<(String, String, usize, usize)> = Vec::new();

    let mut id_flag:bool = false;

    let total_size = reader.records().count();

    let pb = ProgressBar::new(total_size as u64);
    pb.set_style(ProgressStyle::with_template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {bytes}/{total_bytes} ({eta})")
        .unwrap()
        .progress_chars("#>-"));
    

    let reader = fastq::Reader::from_file(fastq_file).unwrap();

    if std::path::Path::new(result_file).exists() {
        fs::remove_file(result_file).unwrap();
      }

    let mut file_ref = fs::OpenOptions::new().create_new(true).append(true).open(result_file).expect("Unable to open file");
    let result_header:String = "readID\tseqID\tscore\thit position\n".to_string();
    file_ref.write_all(result_header.as_bytes()).expect("write failed");

    for read in reader.records() {

        let read_data = read.unwrap();
        
        let read_id: String = read_data.id().to_string();
        // let read_num: usize = read_data.desc().unwrap().split(' ').collect::<Vec<&str>>()[0].parse().unwrap();
        let _read_len: usize = read_data.desc().unwrap().split(' ').collect::<Vec<&str>>()[1].split('=').collect::<Vec<&str>>()[1].parse().unwrap();
        let read_qual: Vec<u8> = read_data.qual().to_vec().iter().map(|x| x-33).collect();
        

        let x: Vec<char> = read_data.seq()
        .to_vec()
        .iter()
        .map(|x| *x as char)
        .collect();
        
        let seq: Vec<SeqElement> = x.iter()
            .map(|x|{
                match x{
                    'A' => SeqElement::A,
                    'G' => SeqElement::G,
                    'T' => SeqElement::T,
                    'C' => SeqElement::C,
                    _ => SeqElement::E,
                }
            })
            .collect();
        let mut matches = query_tree(tree, seq, read_id, percent_match);
        if !matches.is_empty(){
            if !id_flag{
                id_flag = !id_flag;
                pb.println("Positive ID");
            }
            match_set.append(&mut matches);
        }
        for (seq_id, read_id, match_score, hit_pos) in (matches).iter(){
            let out_string:String = format!("{}\t{}\t{}\t{}\n", seq_id, read_id, match_score, hit_pos);
            file_ref.write_all(out_string.as_bytes()).expect("write failed");
        }
        pb.inc(1);
        

    }
}

fn main() {
    std::env::set_var("RUST_BACKTRACE", "1");

    let matches = Command::new("Metagenomic Classification")
        .version("1.0")
        .author("Sriram Vijendran <vijendran.sriram@gmail.com>")
        .subcommand(Command::new("build")
            .about("Build suffix tree index from reference fasta file")
            .arg(arg!(-s --source <SRC_FILE> "Source file with sequences(fasta)")
                .required(true)
                )
            .arg(arg!(-m --max <MAX_DEPTH> "Max depth of the tree")
                .required(true)
                .value_parser(clap::value_parser!(usize))
                )
            .arg(arg!(-o --out <SAVE_FILE> "save file")
                .required(true)
                )
            .arg(arg!(-n --num <NUM_SEQ> "Number of seq. (0==all)")
                .required(true)
                .value_parser(clap::value_parser!(u32))
                )
        )
        .subcommand(Command::new("query")
            .about("Classify reads from fastq file")
            .arg(arg!(-r --reads <READS>"Queries tree from read_file")
                .required(true)
                )
            .arg(arg!(-p --percent_match <PERCENT_MISMATCH>"Percent mismatch to reference sequences")
                .required(true)
                .value_parser(clap::value_parser!(usize))
                )
            .arg(arg!(-t --tree <TREE_FILE>"Queries tree")
                .required(true)
                )
            )
        .subcommand(Command::new("query_dir")
            .about("Classify reads from fastq dir")
            .arg(arg!(-r --read_dir <READS>"directory with barcodes")
                .required(true)
                )
            .arg(arg!(-p --percent_match <PERCENT_MISMATCH>"Percent mismatch to reference sequences")
                .required(true)
                .value_parser(clap::value_parser!(usize))
                )
            .arg(arg!(-t --tree <TREE_FILE>"Queries tree")
                .required(true)
                )
            .arg(arg!(-o --out <OUT_FILE>"output dir")
                .required(true)
                )
            )
        .subcommand(Command::new("query_dir_all")
            .about("Classify reads from fastq dir")
            .arg(arg!(-r --read_dir <READS>"directory with barcodes")
                .required(true)
                )
            .arg(arg!(-p --percent_match <PERCENT_MISMATCH>"Percent mismatch to reference sequences")
                .required(true)
                .value_parser(clap::value_parser!(usize))
                )
            .arg(arg!(-t --tree <TREE_FILE>"Queries tree")
                .required(true)
                )
            .arg(arg!(-o --out <OUT_FILE>"output dir")
                .required(true)
                )
            )
        .subcommand(Command::new("quick_build")
            .about("Build suffix tree index from reference fasta file")
            .arg(arg!(-s --source <SRC_FILE> "Source file with sequences(fasta)")
                .required(true)
                )
            .arg(arg!(-n --num <NUM_SEQ> "Number of seq. (0==all)")
                .required(true)
                .value_parser(clap::value_parser!(u32))
                )
            .arg(arg!(-p --percent_match <PERCENT_MISMATCH>"Percent mismatch to reference sequences")
                .required(true)
                .value_parser(clap::value_parser!(usize))
                )
            .arg(arg!(-r --reads <READS>"Queries tree from read_file")
                .required(true)
                )
            .arg(arg!(-o --out <OUT_FILE>"Output file")
                .required(true)
                )
        )
        .about("Metagenomic classifier using Suffix trees")
        .get_matches();
    
    match matches.subcommand(){
        Some(("build",  sub_m)) => {
            let mut tree: KGST<SeqElement, String> = build_tree(sub_m.get_one::<String>("source").expect("required").as_str(), *sub_m.get_one::<usize>("max").expect("required"), *sub_m.get_one::<u32>("num").expect("required"));
            save_tree(&mut tree, sub_m.get_one::<String>("out").expect("required").to_string());
        },
        Some(("query",  sub_m)) => {
            let mut tree: KGST<SeqElement, String> = load_tree(&sub_m.get_one::<String>("tree").expect("required").to_string());
            let percent_mismatch: f32 = (*sub_m.get_one::<usize>("percent_match").expect("required") as f32)/100.0;
            search_fastq(&mut tree, sub_m.get_one::<String>("reads").expect("required").as_str(), format!("{}.txt", sub_m.get_one::<String>("reads").expect("required")).as_str(), percent_mismatch);

        },
        Some(("quick_build",  sub_m)) => {
            let mut tree: KGST<SeqElement, String> = build_tree(sub_m.get_one::<String>("source").expect("required").as_str(), *sub_m.get_one::<usize>("max").expect("required"), *sub_m.get_one::<u32>("num").expect("required"));
            let percent_mismatch: f32 = (*sub_m.get_one::<usize>("percent_match").expect("required") as f32)/100.0;
            search_fastq(&mut tree, sub_m.get_one::<String>("reads").expect("required").as_str(), sub_m.get_one::<String>("out").expect("required").as_str(), percent_mismatch);
        },
        Some(("query_dir",  sub_m)) => {
            let mut tree: KGST<SeqElement, String> = load_tree(&sub_m.get_one::<String>("tree").expect("required").to_string());
            let percent_mismatch: f32 = (*sub_m.get_one::<usize>("percent_match").expect("required") as f32)/100.0;
            let files = get_files(sub_m.get_one::<String>("read_dir").expect("required").as_str());
            for file in files{
                let save_file = format!("{}/{}.txt", sub_m.get_one::<String>("out").expect("required").as_str(), file.split("/").collect::<Vec<&str>>().last().unwrap());
                search_fastq(&mut tree, &file, &save_file, percent_mismatch);
            }
        },
        Some(("query_dir_all",  sub_m)) => {
            let mut tree: KGST<SeqElement, String> = load_tree(&sub_m.get_one::<String>("tree").expect("required").to_string());
            let percent_mismatch: f32 = (*sub_m.get_one::<usize>("percent_match").expect("required") as f32)/100.0;
            let files = get_files_all(sub_m.get_one::<String>("read_dir").expect("required").as_str());
            for file in files{
                let save_file = format!("{}/{}.txt", sub_m.get_one::<String>("out").expect("required").as_str(), file.split("/").collect::<Vec<&str>>().last().unwrap());
                search_fastq(&mut tree, &file, &save_file, percent_mismatch);
            }
        },
        _ => {
            println!("Either build a tree or query an existing tree. Refer help page (-h flag)");
        }
    }
}
