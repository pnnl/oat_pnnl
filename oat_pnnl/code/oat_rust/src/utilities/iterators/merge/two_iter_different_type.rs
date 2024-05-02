use crate::utilities::partial_order::StrictlyLess;
use crate::utilities::iterators::general::PeekUnqualified;



/// Merges two iterators; if both iterators return entries in ascencing order (according to the order comaparator)
/// then so does the merged iterator.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::iterators::merge::two_iter_different_type::MergeTwoItersByOrderComparator;        
/// use oat_rust::utilities::partial_order::OrderComparatorAutoAnyType;
/// 
/// let vec1 = vec![ 0, 2, 4 ];
/// let vec2 = vec![ 1, 3, 5 ];
/// let merged 
///     = MergeTwoItersByOrderComparator::new( 
///             vec1.iter().peekable(), 
///             vec2.iter().peekable(), 
///             OrderComparatorAutoAnyType
///         );
/// itertools::assert_equal( merged, vec![ 0,1,2,3,4,5 ].iter() );
/// ```
pub struct MergeTwoItersByOrderComparator< Iterator1, Iterator2, OrderComparator >
    where 
        Iterator1:          Iterator + PeekUnqualified,
        Iterator2:          Iterator< Item = Iterator1::Item > + PeekUnqualified,
        OrderComparator:    StrictlyLess<  Iterator1::Item >
{
    iter1:              Iterator1,
    iter2:              Iterator2,
    order_comparator:   OrderComparator,
}

impl < Iterator1, Iterator2, OrderComparator >

    MergeTwoItersByOrderComparator
        < Iterator1, Iterator2, OrderComparator >

    where 
        Iterator1:          Iterator + PeekUnqualified,
        Iterator2:          Iterator< Item = Iterator1::Item > + PeekUnqualified,
        OrderComparator:    StrictlyLess<  Iterator1::Item >        

{
    pub fn new( iter1: Iterator1, iter2: Iterator2, order_comparator: OrderComparator )
            ->
            MergeTwoItersByOrderComparator
                < Iterator1, Iterator2, OrderComparator >                
    {
        MergeTwoItersByOrderComparator{ iter1, iter2, order_comparator }
    }                
}        




impl < Iterator1, Iterator2, OrderComparator >

    Iterator for

    MergeTwoItersByOrderComparator
        < Iterator1, Iterator2, OrderComparator >

    where 
        Iterator1:          Iterator + PeekUnqualified,
        Iterator2:          Iterator< Item = Iterator1::Item > + PeekUnqualified,
        OrderComparator:    StrictlyLess<  Iterator1::Item >        
{
    type Item = Iterator1::Item;
    
    fn next( &mut self ) -> Option< Self::Item > {
        match self.iter1.peek_unqualified() {
            Some( item1 ) => {
                match self.iter2.peek_unqualified() {
                    Some( item2 ) => {
                        match self.order_comparator.strictly_less( item1, item2 ) {
                            true    => { return self.iter1.next() },
                            false   => { return self.iter2.next() },
                        }
                    },
                    None => { return self.iter1.next() }
                }
            },
            None => { return self.iter2.next() }            
        }
    }
}      









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
    fn test_merge_by_order_comparator() {
        use crate::utilities::iterators::merge::two_iter_different_type::MergeTwoItersByOrderComparator;  
        use crate::utilities::partial_order::OrderComparatorAutoAnyType;      

        let vec1 = vec![ 0, 2, 4 ];
        let vec2 = vec![ 1, 3, 5 ];
        let merged 
            = MergeTwoItersByOrderComparator::new( 
                    vec1.iter().peekable(), 
                    vec2.iter().peekable(), 
                    OrderComparatorAutoAnyType
                );
        itertools::assert_equal( merged, vec![ 0,1,2,3,4,5 ].iter() );
    }    
}