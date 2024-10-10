//! Elementary statistics



//  ---------------------------------------------------------------------------
//  COUNTING
//  ---------------------------------------------------------------------------

/// Given an object that implements `Iterator< Item=usize >`, count the number
/// of occurences of each integer.
pub fn  histogram 
        < I: Iterator< Item = usize > > 
        ( iter: I ) 
        -> 
        Vec< usize > 
{
    let mut hist = Vec::new();
    for i in iter.into_iter() {
        while i + 1 > hist.len() { hist.push(0); }
        hist[ i ] +=1;
    }
    hist
}