/*
author: Noah Legall
email: noahaus@uga.edu
purpose: given a multi-sequence alignment, calculate the Tajima's D statistic in a sliding window.
*/
fn main() {
use seq_io::fasta::{Reader,Record}; //https://docs.rs/seq_io/0.3.1/seq_io/
use std::str;

// read in the initial fasta file
let mut reader = Reader::from_path("example.fa").unwrap();


// is there a way to create specific windows of a fasta?
// well, can we get the fasta in a form we want first?
// A vec of Strings is my go to as a novice Rust programmer.
let mut window = Vec::<String>::new();
while let Some(record) = reader.next() {
    let record = record.expect("Error reading record");
    for s in record.seq_lines() {
        // this is where the magic should happen
        // window takes the value of the FASTA record 
        // it comes out as a array of utf numbers. this can be converted into characters.
        window.push(str::from_utf8(&s).unwrap().to_string());
    }

}
// we now successfully have a list of sequences. 
// uncomment the print statement to check.
// println!("{:?}",window);

//TODO: I'll need to figure out how to take the vec of strings and perform a sliding window analysis on it.



}
