/*
author: Noah Legall
email: noahaus@uga.edu
purpose: given a multi-sequence alignment, calculate the Tajima's D statistic in a sliding window.
*/
fn main() {
use seq_io::fasta::{Reader,Record}; //https://docs.rs/seq_io/0.3.1/seq_io/

let mut reader = Reader::from_path("example.fa").unwrap();

let mut n = 0;
let mut sum = 0;
while let Some(record) = reader.next() {
    let record = record.expect("Error reading record");
    for s in record.seq_lines() {
        sum += s.len();
    }
    n += 1;
}
println!("mean sequence length of {} records: {:.1} bp", n, sum as f32 / n as f32);   


}
