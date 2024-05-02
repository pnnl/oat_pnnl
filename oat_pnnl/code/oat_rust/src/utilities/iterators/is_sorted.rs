//! Determines whether an iterator is sorted.
//! 
//! Code for this module is inspired by the crate [`is_sorted`](https://docs.rs/is_sorted/0.1.1/is_sorted/trait.IsSorted.html).
//! See below for the copyright and license.  The oat_rust developers originally attempted to use the `is_sorted` crate as a dependency but encountered build issues.
//! 
//! 
//! 
//! > Copyright (c) 2018 Gonzalo Brito Gadeschi Copyright (c) 2017 The Rust Project Developers
//! >
//! > Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//! > 
//! > The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//! > 
//! > THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.



pub trait IsSortedBy: Iterator {

    /// Returns `true` if the iterator returns item in ascneding order according to `comparator`.
    /// 
    /// # Examples 
    /// 
    /// ```
    /// use oat_rust::utilities::iterators::is_sorted::IsSortedBy;
    /// 
    /// let v = vec![ 1, 2, 3];
    /// assert!( v.iter().is_sorted_by(| a, b | a < b) );
    /// 
    /// let u = vec![ 1, 2, 1];
    /// assert!( ! u.iter().is_sorted_by(| a, b | a < b) );
    /// ```
    fn is_sorted_by<F: FnMut(& Self::Item, & Self::Item) -> bool > (&mut self, mut comparator: F) -> bool {
        match self.next() {
            Some( mut item_second_from_last ) => {
                for item_last in self {
                    if ! comparator( & item_second_from_last, & item_last ){ 
                        return false 
                    } else {
                        item_second_from_last = item_last;
                    }
                }
                return true
            },
            None => { return true }
        }
    }
}

impl < I: Iterator > IsSortedBy for I
{}





//  =========================================================================================================
//  TESTS
//  =========================================================================================================


//  ---------------------------------------------------------------------
//  Doc-test drafts
//  ---------------------------------------------------------------------

//  We use the following module to draft doc tests, which are easier to debug here than in doc strings.

#[cfg(test)]
mod doc_test_drafts {


    #[test]
    fn test_skip_until() {
        use crate::utilities::iterators::is_sorted::IsSortedBy;
        
        let v = vec![ 1, 2, 3];
        assert!( v.iter().is_sorted_by(| a, b | a < b) );
        
        let u = vec![ 1, 2, 1];
        assert!( ! u.iter().is_sorted_by(| a, b | a < b) );    
    }
}