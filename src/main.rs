use bio::io::fasta::Reader;
use generalized_suffix_tree::suffix_tree::KGST;
use indicatif::{ProgressBar, ProgressStyle};


#[derive(Eq, Hash, PartialEq, Copy, Clone, Debug)]
enum SeqElement {
    A, G, T, C, E
}


fn main() {
    std::env::set_var("RUST_BACKTRACE", "1");

    let filename = "/home/sriramv/Datasets/PDesign/database/PA.fasta";
    let max_depth = 5;


    let reader = Reader::from_file(filename).unwrap();

    let total_size = reader.records().count();

    let pb = ProgressBar::new(total_size as u64);
    pb.set_style(ProgressStyle::with_template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {bytes}/{total_bytes} ({eta})")
        .unwrap()
        // .with_key("eta", |state: &ProgressState, w: &mut dyn Write| write!(w, "{:.1}s", state.eta().as_secs_f64()).unwrap())
        .progress_chars("#>-"));
    
    let mut tree: KGST<SeqElement, String> = KGST::new(max_depth, SeqElement::E);

    let reader = Reader::from_file(filename).unwrap();
    
    for result in reader.records() {

        let result_data = result.unwrap();

        let x: Vec<char> = result_data.seq()
        .to_vec()
        .iter()
        .map(|x| *x as char)
        .collect();
        
        if x.contains(&'N'){
            // println!("{}", result_data.id());
        }
        else{
            // println!("{:?}", x);
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
        
            // println!("{}, {}, {}",result_data.id(), s, string_len);
            if string_len>max_depth.try_into().unwrap(){
                let num_iter = string_len+1-(max_depth as usize);
                for (n, depth) in (0..num_iter).enumerate(){
                    tree.add_string(seq[depth..depth+(max_depth as usize)].to_vec(), format!("{}_{}", result_data.id(), n));
                    pb.println(format!("{}_{}", result_data.id(), n));
                }
                
            }
        pb.inc(1);        
        }
                                    
    }
}
