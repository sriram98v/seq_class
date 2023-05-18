use serde::{Serialize, Deserialize};
use std::{fmt};


#[derive(Eq, Hash, PartialEq, Copy, Clone, Debug, Serialize, Deserialize)]
pub enum SeqElement {
    A, G, T, C, E
}

impl fmt::Display for SeqElement {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            SeqElement::A => write!(f, "A"),
            SeqElement::G => write!(f, "G"),
            SeqElement::T => write!(f, "T"),
            SeqElement::C => write!(f, "C"),
            SeqElement::E => write!(f, "$"),
        }
    }
}