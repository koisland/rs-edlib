use std::ops::Index;

use anyhow::{bail, Context};

use crate::align::MAX_UCHAR;

#[derive(Debug, Clone)]
/// Defines two given characters as equal.
pub struct EqualityPair {
    /// First character.
    pub first: char,
    /// Second character.
    pub second: char,
}

/// Defines equality relation on alphabet characters.
#[derive(Debug, Clone)]
pub struct EqualityDefinition {
    pub(crate) alphabet: String,
    matrix: Vec<bool>,
}

impl EqualityDefinition {
    /// Initialize a new `EqualityDefinition`.
    ///
    /// # Arguments
    /// * `alphabet`: All possible characters.
    /// * `added_equalities`: Additional equalities to add. Characters must exist within alphabet.
    ///
    /// # Returns
    /// * New instance of `EqualityDefinition`.
    ///
    /// # Examples
    /// No added equalities.
    /// ```
    /// use rs_edlib::equal::EqualityDefinition;
    ///
    /// let eq_def = EqualityDefinition::new("ATGC", None);
    /// ```
    /// Make `A` also equal `T`.
    /// ```
    /// use rs_edlib::equal::{EqualityDefinition, EqualityPair};
    ///
    /// let eq_def = EqualityDefinition::new(
    ///     "ATGC",
    ///     &[EqualityPair  { first: 'A', second: 'T' }]
    /// );
    /// ```
    pub fn new(alphabet: &str, added_equalities: Option<&[EqualityPair]>) -> Self {
        let mut matrix = vec![false; MAX_UCHAR * MAX_UCHAR];

        // (1,2) = 7
        /*
         * * *
         * * *
         * * *
         */
        for x in 0..alphabet.len() {
            for y in 0..alphabet.len() {
                matrix[(y * alphabet.len()) + x] = x == y
            }
        }

        if let Some(added_equalities) = added_equalities {
            for added_equality in added_equalities.iter() {
                let first_transformed = alphabet.find(added_equality.first);
                let second_transformed = alphabet.find(added_equality.second);
                if let (Some(first_transform_pos), Some(second_transform_pos)) =
                    (first_transformed, second_transformed)
                {
                    // SAFETY: First transformed char and second transformed character guaranteed to exist in alphabet.
                    // Matrix constructed from alphabet.
                    unsafe {
                        *matrix.get_unchecked_mut(
                            first_transform_pos + (alphabet.len() * second_transform_pos),
                        ) = true
                    }
                }
            }
        }
        EqualityDefinition {
            alphabet: alphabet.to_owned(),
            matrix,
        }
    }

    /// Check if characters are equivalent.
    ///
    /// # Arguments
    /// * `a`: Character.
    /// * `b`: Character
    ///
    /// # Examples
    ///
    /// Normal usage.
    /// ```
    /// use rs_edlib::equal::EqualityDefinition;
    ///
    /// let eq_def = EqualityDefinition::new("ATGC", None);
    /// assert!(eq_def.are_equal('T', 'T'));
    /// ```
    ///
    /// Invalid character.
    /// ```should_panic
    /// use rs_edlib::equal::EqualityDefinition;
    /// ///
    /// let eq_def = EqualityDefinition::new("ATGC", None)
    /// // Z not in alphabet.
    /// eq_def.are_equal('Z', 'Z');
    /// ```
    pub fn are_equal(&self, a: char, b: char) -> anyhow::Result<bool> {
        let (Some(pos_x), Some(pos_y)) = (self.alphabet.find(a), self.alphabet.find(b)) else {
            bail!(
                "One or more characters ({a}, {b}) do not exist in alphabet {}",
                self.alphabet
            );
        };
        Ok(self[(pos_x, pos_y)])
    }

    /// Get symbol in alphabet by index.
    ///
    /// # Arguments
    /// * `index`: Index of character in alphabet.
    ///
    /// # Returns
    /// * Character in alphabet, if exists.
    pub fn symbol(&self, index: usize) -> Option<char> {
        self.alphabet.chars().nth(index)
    }
}

/// Index into `EqualityDefinition` matrix by (row, col).
impl Index<(usize, usize)> for EqualityDefinition {
    type Output = bool;

    fn index(&self, index: (usize, usize)) -> &bool {
        self.matrix
            .get(index.0 + self.alphabet.len() * index.1)
            .with_context(|| format!("Invalid index {index:?}."))
            .unwrap()
    }
}

mod test {
    use super::*;

    #[test]
    fn test_init_equality_definition() {
        let eq_def = EqualityDefinition::new("ATGC", None);
        assert_eq!(eq_def.matrix.len(), 16);
        assert_eq!(eq_def.matrix.iter().filter(|elem| **elem).count(), 4)
    }

    #[test]
    fn test_init_equality_definition_extra_def() {
        let eq_def = EqualityDefinition::new(
            "ATGC",
            Some(&[EqualityPair {
                first: 'A',
                second: 'T',
            }]),
        );
        assert_eq!(eq_def.matrix.len(), 16);
        // First A and second T now count as equal.
        assert_eq!(eq_def.matrix.iter().filter(|elem| **elem).count(), 5)
    }

    #[test]
    fn test_equality_definition_check_equal() {
        let eq_def = EqualityDefinition::new("ATGC", None);
        assert!(eq_def.are_equal('T', 'T').unwrap());
        assert!(!eq_def.are_equal('T', 'G').unwrap());
    }

    #[test]
    fn test_equality_definition_check_invalid_equal() {
        let eq_def = EqualityDefinition::new("ATGC", None);
        assert!(eq_def.are_equal('X', 'X').is_err());
    }

    #[test]
    fn test_equality_definition_index() {
        let eq_def = EqualityDefinition::new("ATGC", None);
        // T == T. T is at position 1.
        assert!(eq_def[(1, 1)]);
        // T != G. T is at position 1. G is at position 2.
        assert!(!eq_def[(1, 2)]);
    }

    #[test]
    #[should_panic]
    fn test_equality_invalid_definition_index() {
        let eq_def = EqualityDefinition::new("ATGC", None);
        // 7 is not a valid position in matrix.
        assert!(eq_def[(7, 7)]);
    }
}
