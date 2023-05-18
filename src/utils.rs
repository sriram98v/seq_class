use glob::{glob_with, MatchOptions};

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