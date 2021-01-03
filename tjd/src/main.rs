/*
author: Noah Legall
email: noahaus@uga.edu
purpose: given a multi-sequence alignment, calculate the Tajima's D statistic in a sliding window.
*/
fn main() {
use seq_io::fasta::{Reader,Record}; //https://docs.rs/seq_io/0.3.1/seq_io/
use std::str;

// read in the initial fasta file
let mut reader = Reader::from_path("example2.fa").unwrap();


// is there a way to create specific windows of a fasta?
// well, can we get the fasta in a form we want first?
// A vec of Strings is my go to as a novice Rust programmer.
let mut alignment = Vec::<String>::new();
while let Some(record) = reader.next() {
    let record = record.expect("Error reading record");
    for s in record.seq_lines() {
        // this is where the magic should happen
        // window takes the value of the FASTA record 
        // it comes out as a array of utf numbers. this can be converted into characters.
        alignment.push(str::from_utf8(&s).unwrap().to_string());
    }

}

let temp = calculate_pi(&alignment);
println!("{}",temp);

// we now successfully have a list of sequences. 
// uncomment the print statement to check.
// println!("{:?}",window);

// let's create the window size first.
/*let mut window = Vec::<String>::new();
let mut start: usize = 0;
let mut end: usize = 3;
while end <= alignment[0].len(){
    for row in &alignment {
        //the 'row' value needs to be referenced here to be sliced,
        //println!("{:?}",&row[0..3]);
        window.push(row[start..end].to_string());
}

// move the window up one for this example.
// println!("{:?}",window);
let temp = calculate_pi(&window);
println!("{}",temp);
start = start + 1;
end = end + 1;
window.clear();
}*/

}

// given a vector of strings, calculate pi, which is calculated according to this website:
// https://arundurvasula.wordpress.com/2015/02/18/interpreting-tajimas-d/
fn calculate_pi(alignment: &Vec<String>) -> f64 {
    let align_length = alignment.len();
    let mut distances = Vec::new();
    for i in 0..align_length {
        for j in 0..align_length {
            if i < j {
                distances.push(hamming_distance(&alignment[i],&alignment[j]));
            } else {
                continue;
            }
        }
    }
    let align_length = align_length as i32;
    let dist_sum = distances.iter().sum::<i32>();
    let pi = (dist_sum*2) as f64 /(align_length*(align_length-1)) as f64;
    return pi;
    
}

// given two strings of equal size, calculate the number of mismatching characters.
fn hamming_distance(seq1: &String,seq2: &String) -> i32 {
    let mut distance = 0;
    for i in 0..seq1.len() {
        if seq1.chars().nth(i) == seq2.chars().nth(i) {            
            continue;            
        } else {
            distance += 1;
        }
    }
    return distance;
}

// TODO: create 'calculate_theta' method and add sliding window support