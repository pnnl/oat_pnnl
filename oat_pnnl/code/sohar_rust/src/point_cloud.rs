//! Functions to generate point clouds


use std::{f64::consts::PI, ops::Range};

use rand::Rng;

/// Return `m` points evenly spaced on the unit circle
pub fn unit_circle( m: usize, noise_range: Option<Range<f64>> ) -> Vec< Vec< f64 > > {
    let mut rng = rand::thread_rng();

    let circpoint = |k: usize | { 
        let theta = k as f64 * 2.0 * PI / m as f64; 
        let mut vec = vec![ theta.cos(), theta.sin() ]; 
        return vec          
    };
    let mut vec: Vec<_> = (0..m).map(|k| circpoint(k) ).collect();
    if let Some( r ) = noise_range {
        for row in vec.iter_mut() {
            for t in row.iter_mut() { *t += rng.gen_range( r.clone()) }
        }
    }           
    return  vec
}
