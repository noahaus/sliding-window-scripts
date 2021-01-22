/*
author: Noah Legall
email: noahaus@uga.edu
purpose: given a multi-sequence alignment, calculate the Tajima's D statistic in a sliding window.
arguments: input FASTA, output CSV, window size, window step size
*/

/****
TODOS
*****/

// practice using rust wrapper with python https://developers.redhat.com/blog/2017/11/16/speed-python-using-rust/


/******************** 
 main logic of code 
 ********************/

 use seq_io::fasta::Reader; //https://docs.rs/seq_io/0.3.1/seq_io/ we use this external crate to make FASTA parsing a bit easier.
 use std::str;
 use std::fs::File;
 use std::io::Write;
 use std::env;
 
fn main() -> std::io::Result<()> {

let args: Vec<String> = env::args().collect();

let input_path = &args[1];
let output_path = &args[2];
let window_size = &args[3].parse::<usize>().unwrap();
let window_step = &args[4].parse::<usize>().unwrap();

let mut output = File::create(output_path)?;

//// 1 read in fasta file into Rust
// 1.0 read in the fasta given a proper path
let mut reader = Reader::from_path(input_path).unwrap();

// 1.1 create a Vec data structure that will hold the sequence strings
let mut alignment = Vec::<String>::new();
println!("start,end,tjd"); // this will be important for CSV file headers
writeln!(output,"start,end,tjd")?;

// 1.2 loop through the file and add values to the 'alignment' variable.
// this is pretty boiler-plate from the seq_io website.
println!("reading in alignment");
while let Some(record) = reader.next() {
    let record = record.expect("Error reading record");
    for s in record.seq_lines() {
        alignment.push(str::from_utf8(&s).unwrap().to_string());
    }
}
println!("alignment stored");

//// 2 generate a sliding window of the alignment data
// 2.0 initialize the variables we need to create a sliding window.
// currently we need a 'step_size' value, but this can be integrated trivially
let mut window = Vec::<String>::new(); // vector to store window
let mut start: usize = 0; // index to start the window
let mut end: usize = start + window_size; // index to end the window
let mut data = String::new();

// 2.1 while the size of 'end' is smaller than total size of the sequence, update the window

while end <= alignment[0].len(){
    for row in &alignment {
        window.push(row[start..end].to_string());
        //println!("window length pushed. end value equals {}", end);
}

// 2.2 calculate the tajima's d statistic and print
let temp = tajimas_d(&window);
println!("tajima's d calculated");
println!("{},{},{}",start,end,temp);
data = format!("{},{},{}",start,end,temp);
writeln!(output,"{}",data)?;

// 2.3 update the step of the variables. reset the 'window' variable to be empty
start += window_step;
end += window_step;
window.clear();
}

Ok(())

}

/************************
 supplemental functions 
 ************************/

// 1 given an alignment, will calculate the tajima's d statistic.
fn tajimas_d(align: &Vec<String>) -> f64 {
    return calculate_pi(&align) - calculate_theta(&align);
}

// 1 given an alignment, will calculate pairwise variation between the sequences (pi)
fn calculate_pi(alignment: &Vec<String>) -> f64 {
    //// 1.0 initialize the variables needed for pairwise comparison
    let align_length = alignment.len(); // number of sequences being analyzed
    let mut distances = Vec::new(); // vector that holds all the values of the pairwise distances

    // 1.1 'distances' vector is updated as the hamming distance is calculated between two sequences
    println!("prepare for pairwise comparisons");
    for i in 0..align_length {
        for j in 0..align_length {
            if i < j {
                distances.push(hamming_distance(&alignment[i],&alignment[j]));
            } else {
                continue;
            }
        }
    }
    println!("pairwise distances computed");

    // 1.2 sum the pairwise distances and then divide by the n(n-1)
    let align_length = align_length as i32;
    let dist_sum = distances.iter().sum::<i32>();
    let pi = (dist_sum*2) as f64 /(align_length*(align_length-1)) as f64;
    println!("pi computed");
    return pi;

    //// 2 given two string sequences, output the hamming distance between them
    // assuming that these sequences are the same size, then we can just calculate the number of times the strings are different
    fn hamming_distance(seq1: &String,seq2: &String) -> i32 {
        //Remember this 
        let distance = seq1.as_bytes().iter().zip(seq2.as_bytes()).fold(0,|a,(b,c)| a + (*b != *c) as i32);
        return distance;
    }
    
}

//// 1 given vector of sequences, calculate the value theta, which is an estimate of equilibrium.
fn calculate_theta(alignment: &Vec<String>) -> f64 {
    // 1.0 initialize the variables
    let mut seg_sites = 0; // the number of segregating site, columns that are polymorphic
    let mut align_column = String::new(); // string to represent an alignment column
    let mut nuc_set = Vec::<char>::new(); // the nucleotides that make up a column. will be used to calculate a number to compare to '1'
    let mut a: f64 = 0.0; // this quantity is used to help calculate theta
    
    // 1.1 loop through the alignment and extract the columns
    // do calculation based on these columns
    for i in 0..alignment[0].len() {

        // 1.2 for a given column, get all the nucleotides
        for j in 0..alignment.len() {
            align_column.push(alignment[j].chars().nth(i).unwrap());
        }

        // 1.3 remove duplicate nucleotides from the column 
        // vectors can remove duplicates, but the vector needs to be sorted first
        nuc_set = align_column.chars().collect::<Vec<char>>();
        nuc_set.sort();
        nuc_set.dedup();

        // 1.4 based on the number of unique characters in the column, update value of segregating sites
        // after the calculation, reset the column variable
        if nuc_set.len() > 1 {
            seg_sites += 1;
        } else {
            continue;
        }
        align_column.clear();
    }

    // 1.5 update the value of 'a'
    for i in 1..alignment.len() {
        a += 1.0/(i as f64); 
    }

    return (seg_sites as f64)/a;

}