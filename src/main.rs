extern crate clap;
use clap::{arg, Command};
use std::collections::{HashMap, HashSet};
use bio::io::{fasta,  fastq};
use generalized_suffix_tree::suffix_tree::KGST;
use indicatif::{ProgressBar, ProgressStyle};
use serde::{Serialize, Deserialize};
use std::fs;
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
            
            strings.insert(format!("{}\n", result_data.id()), seq);
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

fn query_tree(tree:&mut KGST<SeqElement, String>, q_seq:Vec<SeqElement>, max_depth:i32)->HashSet<(String, i32)>{
    let mut match_set:HashSet<(String, i32)> = HashSet::new();
    let string_len = q_seq.len();
    if string_len>=max_depth.try_into().unwrap(){
        let num_iter = string_len+1-(max_depth as usize);
        for  depth in 0..num_iter{
            let sub_seq = q_seq[depth..depth+(max_depth as usize)].to_vec();
            let matches: Vec<(&String, &i32)> = tree.find(sub_seq);
            if !matches.is_empty(){
                for (hit_id, hit_idx) in matches.iter(){
                    match_set.insert(((**hit_id).clone(), **hit_idx));
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
}

fn load_tree(fname:&String) -> KGST<SeqElement, String>{
    println!("Loading tree from {}", fname);
    let json_str:String = fs::read_to_string(fname).unwrap();
    let tree: KGST<SeqElement, String> = serde_json::from_str(&json_str).unwrap();
    tree
}

fn search_fastq(tree:&mut KGST<SeqElement, String>, fastq_file:&str, max_depth:i32, result_file:&str){
    let reader = fastq::Reader::from_file(fastq_file).unwrap();

    let mut match_set:HashSet<(String, i32)> = HashSet::new();

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
            let mut matches = query_tree(tree, seq, max_depth);
            if !matches.is_empty(){
                if !id_flag{
                    id_flag = !id_flag;
                    pb.println("Positive ID");
                }
                match_set.extend(matches.clone());
            }
            for (id, idx) in matches.iter(){
                let out_string:String = format!("{}\t{}\t{}\n", result_data.id(), id, idx);
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
        .arg(arg!(
            -b --build <FILE> "Builds tree from file"
        )
        .required(false)
        )
        .arg(arg!(
            -m --max-depth <MAX_DEPTH> "Max depth of the tree"
        )
        .required(false)
        )
        .arg(arg!(
            -s --search <READ_FILE>"Queries tree from read_file"
        )
        .required(false)
        )
        .arg(arg!(
            -t --tree-file <TREE_FILE>"Queries tree"
        )
        .required(false)
        )
        .arg(arg!(
            -o --out-file <OUT_FILE>"Output file"
        )
        .required(false)
        )
        .get_matches();

    // if let Some(fname) = matches.get_one::<String>("build") {
    //         println!("Value for name: {}", fname);
    //     }
    
    // let mut tree: KGST<SeqElement, String> = build_tree(matches.get_one::<String>("file").expect("required").as_str(), *matches.get_one::<i32>("max_depth").expect("required"));

    // save_tree(&mut tree, matches.get_one::<String>("save_file").expect("required").to_string());

    // println!("Searching {} for matches", matches.get_one::<String>("search_file").expect("required").as_str());
    // search_fastq(&mut tree, matches.get_one::<String>("search_file").expect("required").as_str(), *matches.get_one::<i32>("max_depth").expect("required"), matches.get_one::<String>("result_file").expect("required").as_str());
}
