use glob::{glob_with, MatchOptions};
use crate::alphabet::*;

pub fn get_files(dir:&str) -> Vec<String>{
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

pub fn get_files_all(dir:&str) -> Vec<String>{
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

pub fn match_prob(_q_seq_match_vec: Vec<bool>, _quality_score_vec: Vec<u8>) -> f32{
    0.0
}

pub fn complement(q_seq: Vec<SeqElement>)->Vec<SeqElement>{
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