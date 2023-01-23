extern crate clap;
use clap::{arg, Command};
use std::collections::{HashMap, HashSet};
use bio::io::{fasta,  fastq};
use generalized_suffix_tree::suffix_tree::KGST;
use indicatif::{ProgressBar, ProgressStyle};
use serde::{Serialize, Deserialize};
use std::{fs, string};
use std::io::Write;


#[derive(Eq, Hash, PartialEq, Copy, Clone, Debug, Serialize, Deserialize)]
enum SeqElement {
    A, G, T, C, E
}

fn build_tree(file:&str, max_depth:i32)->KGST<SeqElement, String>{
    println!("Building tree from {}", file);
    let reader = fasta::Reader::from_file(file).unwrap();

    let total_size = reader.records().count();

    let pb = ProgressBar::new(total_size as u64);
    pb.set_style(ProgressStyle::with_template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {bytes}/{total_bytes} ({eta})")
        .unwrap()
        // .with_key("eta", |state: &ProgressState, w: &mut dyn Write| write!(w, "{:.1}s", state.eta().as_secs_f64()).unwrap())
        .progress_chars("#>-"));
    
    let mut tree: KGST<SeqElement, String> = KGST::new(SeqElement::E);

    let reader = fasta::Reader::from_file(file).unwrap();

    let mut count = 20;
    
    for result in reader.records() {

        let result_data = result.unwrap();

        let x: Vec<char> = result_data.seq()
        .to_vec()
        .iter()
        .map(|x| *x as char)
        .collect();
        
        if x.contains(&'N'){
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
            let string_len = seq.len();
        
            if string_len>max_depth.try_into().unwrap(){
                let num_iter = string_len+1-(max_depth as usize);
                for (n, depth) in (0..num_iter).enumerate(){
                    tree.add_string(seq[depth..depth+(max_depth as usize)].to_vec(), format!("{}_{}", result_data.id(), n));
                }
                
            }
            // pb.println(result_data.id());
            pb.inc(1);   
            count+=1;
            if(count%20==0){
                break;
            }
        }
    }

    let reader = fasta::Reader::from_file(file).unwrap();
    let mut strings:HashMap<String, Vec<SeqElement>> = HashMap::new();
    
    for result in reader.records() {

        let result_data = result.unwrap();

        let x: Vec<char> = result_data.seq()
        .to_vec()
        .iter()
        .map(|x| *x as char)
        .collect();
        
        if x.contains(&'N'){
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
            
            strings.insert(format!("{}", result_data.id()), seq);
        }
    }
    tree.set_strings(strings);
    tree
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

fn hamming_distance(ref_seq: &Vec<SeqElement>, read_seq: &Vec<SeqElement>)->usize{
    let mut count:usize = 0;
    let it = ref_seq.iter().zip(read_seq.iter());
    for (x,y) in it{
        if x!=y{count +=1;}
    }
    count
}

fn query_tree(tree:&KGST<SeqElement, String>, q_seq:Vec<SeqElement>, q_seq_id:String, max_depth:i32)->HashSet<(String, String, usize)>{
    let mut match_set:HashSet<(String, String, usize)> = HashSet::new();
    let string_len = q_seq.len();
    if string_len>=max_depth.try_into().unwrap(){
        let num_iter = string_len+1-(max_depth as usize);
        for (n,depth) in (0..num_iter).enumerate(){
            let sub_seq = q_seq[depth..depth+(max_depth as usize)].to_vec();
            // println!("{:?}", sub_seq);
            let matches: Vec<(&String, &i32)> = tree.find(sub_seq);
            // println!("{:?}", matches);
            if !matches.is_empty(){
                for (hit_id, hit_idx) in matches.iter(){
                    let hit_pos: usize = hit_id.split('_').collect::<Vec<&str>>()[1].parse().unwrap();
                    let ref_seq:&Vec<SeqElement> = tree.get_string(&(**hit_id).split('_').collect::<Vec<&str>>()[0].to_string());
                    if &n<=&hit_pos && (&hit_pos+&string_len)<=ref_seq.len(){
                        let ref_sice: Vec<SeqElement> = ref_seq[hit_pos..hit_pos+string_len].to_vec();
                        // println!("{}", &(**hit_id).split('_').collect::<Vec<&str>>()[0].to_string());
                        let dist:usize = hamming_distance(&ref_sice, &q_seq);
                        match_set.insert(((**hit_id).split('_').collect::<Vec<&str>>()[0].to_string(), q_seq_id.clone(), dist));
                    }
                }
            }
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

fn search_fastq(tree:&KGST<SeqElement, String>, fastq_file:&str, max_depth:i32, result_file:&str){
    println!("Classifying read file: {}", &fastq_file);
    let reader = fastq::Reader::from_file(fastq_file).unwrap();

    let mut match_set:HashSet<(String, String, usize)> = HashSet::new();

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

    for read in reader.records() {

        let read_data = read.unwrap();
        
        let read_id: String = read_data.id().to_string();
        let read_num: i32 = read_data.desc().unwrap().split(' ').collect::<Vec<&str>>()[0].parse().unwrap();
        let read_len: i32 = read_data.desc().unwrap().split(' ').collect::<Vec<&str>>()[1].split('=').collect::<Vec<&str>>()[1].parse().unwrap();
        let read_qual: Vec<u8> = read_data.qual().to_vec();
        

        let x: Vec<char> = read_data.seq()
        .to_vec()
        .iter()
        .map(|x| *x as char)
        .collect();
        
        if x.contains(&'N'){
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
            let mut matches = query_tree(tree, seq, read_id, max_depth);
            if !matches.is_empty(){
                if !id_flag{
                    id_flag = !id_flag;
                    pb.println("Positive ID");
                }
                match_set.extend(matches.clone());
            }
            for (seq_id, read_id, seq_idx) in matches.iter(){
                let out_string:String = format!("{}\t{}\t{}\n", seq_id, read_id, seq_idx);
                pb.println("wrote output");
                file_ref.write_all(out_string.as_bytes()).expect("write failed");
            }
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
                .value_parser(clap::value_parser!(i32))
                )
            .arg(arg!(-o --out <SAVE_FILE> "save file")
                .required(true)
                )
        )
        .subcommand(Command::new("query")
            .about("Classify reads from fastq file")
            .arg(arg!(-r --reads <READS>"Queries tree from read_file")
                .required(true)
                )
            .arg(arg!(-m --max <MAX_DEPTH> "Max depth of the tree")
                .required(true)
                .value_parser(clap::value_parser!(i32))
                )
            .arg(arg!(-t --tree <TREE_FILE>"Queries tree")
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
            let mut tree: KGST<SeqElement, String> = build_tree(sub_m.get_one::<String>("source").expect("required").as_str(), *sub_m.get_one::<i32>("max").expect("required"));
            save_tree(&mut tree, sub_m.get_one::<String>("out").expect("required").to_string());
        },
        Some(("query",  sub_m)) => {
            let mut tree: KGST<SeqElement, String> = load_tree(&sub_m.get_one::<String>("tree").expect("required").to_string());
            search_fastq(&mut tree, sub_m.get_one::<String>("reads").expect("required").as_str(), *sub_m.get_one::<i32>("max").expect("required"), sub_m.get_one::<String>("out").expect("required").as_str());

        },
        _ => {
            println!("Either build a tree or query an existing tree");
        }
    }
    // let mut tree: KGST<SeqElement, String> = build_tree(matches.get_one::<String>("file").expect("required").as_str(), *matches.get_one::<i32>("max_depth").expect("required"));

    // save_tree(&mut tree, matches.get_one::<String>("save_file").expect("required").to_string());

    // println!("Searching {} for matches", matches.get_one::<String>("search_file").expect("required").as_str());
    // search_fastq(&mut tree, matches.get_one::<String>("search_file").expect("required").as_str(), *matches.get_one::<i32>("max_depth").expect("required"), matches.get_one::<String>("result_file").expect("required").as_str());
}
