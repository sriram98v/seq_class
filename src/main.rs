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
    let reader = fasta::Reader::from_file(file).unwrap();

    let total_size = reader.records().count();

    let pb = ProgressBar::new(total_size as u64);
    pb.set_style(ProgressStyle::with_template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {bytes}/{total_bytes} ({eta})")
        .unwrap()
        // .with_key("eta", |state: &ProgressState, w: &mut dyn Write| write!(w, "{:.1}s", state.eta().as_secs_f64()).unwrap())
        .progress_chars("#>-"));
    
    let mut tree: KGST<SeqElement, String> = KGST::new(SeqElement::E);

    let reader = fasta::Reader::from_file(file).unwrap();
    let mut count:i32 = 0;
    
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
            pb.println(format!("{}", result_data.id()));
            pb.inc(1);   
            count+=1;
            if(count%20==0){
                break;
            }
        }
    }

    let reader = fasta::Reader::from_file(file).unwrap();
    let mut strings:HashMap<String, Vec<SeqElement>> = HashMap::new();
    let mut count:u32 = 0;
    
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
    // tree.set_strings(strings);
    tree
}


fn search_fastq(mut tree:KGST<SeqElement, String>, fastq_file:&str, max_depth:i32, result_file:&str){
    let reader = fastq::Reader::from_file(fastq_file).unwrap();

    let mut match_set:HashSet<String> = HashSet::new();

    let mut id_flag:bool = false;

    let total_size = reader.records().count();

    let pb = ProgressBar::new(total_size as u64);
    pb.set_style(ProgressStyle::with_template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {bytes}/{total_bytes} ({eta})")
        .unwrap()
        .progress_chars("#>-"));
    

    let reader = fastq::Reader::from_file(fastq_file).unwrap();
    let mut count:u32 = 0;

    let mut fileRef = fs::OpenOptions::new().create_new(true).append(true).open(result_file).expect("Unable to open file");   



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
            if string_len>=max_depth.try_into().unwrap(){
                let num_iter = string_len+1-(max_depth as usize);
                for (n, depth) in (0..num_iter).enumerate(){
                    let sub_seq = seq[depth..depth+(max_depth as usize)].to_vec();
                    let matches: Vec<(&String, &i32)> = tree.find(sub_seq);
                    if !matches.is_empty(){
                        if !id_flag{
                            id_flag = !id_flag;
                            pb.println(format!("Positive ID"));
                        }
                        for (hit_ID, idx) in matches.iter(){
                            let mut ID:Vec<&str> = hit_ID.split('_').collect();
                            if !match_set.contains(&format!("{:?}", ID[0])){
                                match_set.insert(format!("{:?}", ID[0]));
                                // fileRef.write_all(format!("{:?}\n", ID[0]).as_bytes()).expect("write failed");

                                // fs::write(result_file, format!("{:?}", ID[0])).expect("Unable to write file");
                            }
                            
                        }
                        
                    }
                }
                for id in match_set.iter(){
                    let mut out_string:String = format!("{}\t{}\n", result_data.id(), id);
                    fileRef.write_all(out_string.as_bytes()).expect("write failed");
                }
                
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
        .arg(arg!([file] "fasta file with reference sequences").required(true))
        .arg(arg!([max_depth] "max depth of suffix tree").required(true).value_parser(clap::value_parser!(i32)))
        .arg(arg!([search_file] "Fastq file").required(true))
        .arg(arg!([result_file] "result file path").required(true))
        .get_matches();

    println!("Building tree from {}", matches.get_one::<String>("file").expect("required").as_str());
    let mut tree: KGST<SeqElement, String> = build_tree(matches.get_one::<String>("file").expect("required").as_str(), *matches.get_one::<i32>("max_depth").expect("required"));

    println!("Searching {} for matches", matches.get_one::<String>("search_file").expect("required").as_str());
    search_fastq(tree, matches.get_one::<String>("search_file").expect("required").as_str(), *matches.get_one::<i32>("max_depth").expect("required"), matches.get_one::<String>("result_file").expect("required").as_str());

}
